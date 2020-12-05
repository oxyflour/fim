#include "fit.h"
#include "generated/inc_chunk_cu.h"
#include "utils/assert.h"

#include <regex>
#include <filesystem>
using namespace std;

fit::Port::Port(Grid &grid, cst::port_type &port, float epsi) {
    src = port.src;
    dst = port.dst;
    pos = grid.ParsePort(src, dst, epsi);
    auto c = pos.size() / 2;
    auto d = pos[c] - pos[c - 1];
    idx = grid.GetFlatIndex(pos[c - 1], 0);
    dir = d.x ? 0 : d.y ? 1 : 2;
};

fit::Coefficient::Coefficient(Matrix &mat, float dt) {
    grid = mat.grid;
    auto nvar = grid->xs.size() * grid->ys.size() * grid->zs.size() * 3;
    le = new float[nvar];
    re = new float[nvar];
    lh = new float[nvar];
    rh = new float[nvar];

    float kap = 0, rho = 0;
    for (int i = 0; i < nvar; i ++) {
        auto ep = mat.eps[i], mu = mat.mue[i];
        le[i] = (1 - kap * ep * dt / 2) / (1 + kap * ep * dt / 2);
        re[i] = dt * ep / (1 + kap * ep * dt / 2);
        lh[i] = (1 - rho * mu * dt / 2) / (1 + rho * mu * dt / 2);
        rh[i] = dt * mu / (1 + rho * mu * dt / 2);
    }
}

fit::Coefficient::~Coefficient() {
    delete le;
    delete re;
    delete lh;
    delete rh;
}

void fit::Coefficient::UpdateFromPort(fit::Port &port) {
    auto &pos = port.pos;
    for (int i = 0, len = pos.size(), c = len / 2; i < len - 1; i ++) {
        auto d = pos[i + 1] - pos[i];
        auto g = grid->GetFlatIndex(pos[i], d.x ? 0 : d.y ? 1 : 2);
        if (i != c) {
            //le[g] = 1;
            //re[g] = 0;
        }
    }
}

static auto Compile(Grid &grid, fit::Port &port) {
    auto source = string(SRC_CHUNK_CU);
    source = regex_replace(source, regex("_\\$i(\\W)"), "_0$1");
    source = regex_replace(source, regex("\\$nx(\\W)"), to_string(grid.xs.size()) + "$1");
    source = regex_replace(source, regex("\\$ny(\\W)"), to_string(grid.ys.size()) + "$1");
    source = regex_replace(source, regex("\\$nz(\\W)"), to_string(grid.zs.size()) + "$1");
    source = regex_replace(source, regex("\\$sg(\\W)"), to_string(port.idx) + "$1");
    source = regex_replace(source, regex("\\$sd(\\W)"), to_string(port.dir) + "$1");

    auto tmp = filesystem::temp_directory_path() / ("fit-compile-" + utils::random(8));
    auto inc = (tmp / "include").u8string(),
        src = (tmp / "tpl.cu").u8string(),
        dll = (tmp / "tpl.dll").u8string(),
        log = (tmp / "tpl.log").u8string(),
        cmd = "nvcc -I" + inc + " " + src + " --shared -o " + dll +" > " + log + " 2>&1";

    filesystem::create_directories(tmp);
    utils::writeFile(src, source);
    for (int i = 0, n = sizeof(INC_FILES) / sizeof(char *) - 1; i < n; i +=2) {
        auto path = tmp / INC_FILES[i];
        filesystem::create_directories(path.parent_path());
        utils::writeFile(path.u8string(), INC_FILES[i + 1]);
    }

    auto code = system(cmd.c_str());
    ASSERT(code == 0, "compile " + dll + " failed,\ncmd: " + cmd + "\nmessage: " + utils::readFile(log));

    return new utils::DLL(dll);
}

fit::Solver::Solver(Matrix &mats, float dt, vector<Port> &ports) {
    auto coe = Coefficient(mats, dt);
    for (auto &port : ports) {
        coe.UpdateFromPort(port);
    }

    auto &grid = *mats.grid;
    ASSERT(ports.size() == 1, "FIXME: only one port supported");
    dll = Compile(grid, ports[0]);
    FnInit = (decltype(FnInit)) dll->getProc("init_0");
    FnStep = (decltype(FnStep)) dll->getProc("step_0");
    FnQuit = (decltype(FnQuit)) dll->getProc("step_0");
    ASSERT(FnInit && FnStep && FnQuit, "get function failed");

    ASSERT(FnInit(coe.le, coe.re, coe.lh, coe.rh) == 0, "solver init failed");
}

fit::Solver::~Solver() {
    FnQuit();
    auto path = dll->path;
    delete dll;
    filesystem::remove_all(filesystem::path(path).parent_path());
}

float fit::Solver::Step(float s) {
    return FnStep(s);
}
