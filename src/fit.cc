#include <regex>
#include "fit.h"
#include "utils/assert.h"
using namespace std;

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
            le[g] = 1;
            re[g] = 0;
        }
    }
}

static auto Compile(Grid &grid, fit::Port &port) {
    auto root = utils::dirname(utils::dirname(__FILE__)),
        source = utils::readFile(root + "/src/chunk.cu");
    source = regex_replace(source, regex("_\\$i(\\W)"), "_0$1");
    source = regex_replace(source, regex("\\$nx(\\W)"), to_string(grid.xs.size()) + "$1");
    source = regex_replace(source, regex("\\$ny(\\W)"), to_string(grid.ys.size()) + "$1");
    source = regex_replace(source, regex("\\$nz(\\W)"), to_string(grid.zs.size()) + "$1");
    source = regex_replace(source, regex("\\$sg(\\W)"), to_string(port.idx) + "$1");
    source = regex_replace(source, regex("\\$sd(\\W)"), to_string(port.dir) + "$1");

    auto inc = root + "/include",
        src = root + "/build/tpl.cu",
        dll = root + "/build/tpl.dll",
        log = root + "/build/tpl.log",
        cmd = "nvcc -I" + inc + " " + src + " --shared -o " + dll +" > " + log + " 2>&1";
    utils::writeFile(src, source);
    auto code = system(cmd.c_str());
    ASSERT(code == 0, "compile " + dll + " failed,\ncmd: " + cmd + "\nmessage: " + utils::readFile(log));
    return dll;
}

fit::Solver::Solver(Matrix &mats, float dt, vector<Port> &ports) {
    ASSERT(ports.size() == 1, "FIXME: only one port fupported");
    auto coe = Coefficient(mats, dt);
    for (auto &port : ports) {
        coe.UpdateFromPort(port);
    }

    auto &grid = *mats.grid;
    dll = new utils::DLL(Compile(grid, ports[0]));
    FnInit = (InitPTR) dll->getProc("init_0");
    FnStep = (StepPTR) dll->getProc("step_0");
    FnQuit = (QuitPTR) dll->getProc("step_0");
    ASSERT(FnInit && FnStep && FnQuit, "get function failed");

    ASSERT(FnInit(coe.le, coe.re, coe.lh, coe.rh) == 0, "solver init failed");
}

fit::Solver::~Solver() {
    FnQuit();
    delete dll;
}

float fit::Solver::Step(float s) {
    return FnStep(s);
}
