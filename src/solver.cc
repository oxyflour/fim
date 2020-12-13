#include <regex>
#include <iostream>
#include <fstream>
#include <filesystem>

#include "generated/inc_chunk_cu.h"
#include "solver.h"
#include "check.h"

using namespace std;
using namespace fit;
using namespace grid;
using namespace solver;

Solver::Solver(Matrix &mats, float dt, vector<Port> &ports) {
    coe = new Coefficient(mats, dt);
    for (auto &port : ports) {
        coe->Add(port);
    }
}

Solver::~Solver() {
    delete coe;
}

float Solver::Step(float s) {
    CHECK(FALSE, "step method of base solver should not be called");
    return s;
}

CPUSolver::CPUSolver(Matrix &mats, float dt, vector<Port> &ports) : Solver(mats, dt, ports) {
    grid = mats.grid;

    CHECK(ports.size() == 1, "only one port supported");
    idx = ports[0].idx;
    dir = ports[0].dir;

    auto nvar = grid->nvar, nxyz = grid->nxyz;
    E = new float[nvar]; fill(E, E + nvar, 0.f);
    H = new float[nvar]; fill(H, H + nvar, 0.f);
    Ex = E; Ey = Ex + nxyz; Ez = Ey + nxyz;
    Hx = H; Hy = Hx + nxyz; Hz = Hy + nxyz;
    LEx = coe->le; LEy = LEx + nxyz; LEz = LEy + nxyz; REx = coe->re; REy = REx + nxyz; REz = REy + nxyz;
    LHx = coe->lh; LHy = LHx + nxyz; LHz = LHy + nxyz; RHx = coe->rh; RHy = RHx + nxyz; RHz = RHy + nxyz;
}

float CPUSolver::Step(float s) {
    auto get_idx = [&](int i, int j, int k){ return grid->GetIndex(i, j, k, 0); };
    auto SG = idx, SD = dir;

    for (int i = 1; i < grid->nx - 1; i ++) {
        for (int j = 1; j < grid->ny - 1; j ++) {
            for (int k = 1; k < grid->nz - 1; k ++) {
                auto g = get_idx(i, j, k);
                Hx[g] = LHx[g] * Hx[g] - RHx[g] * (Ey[get_idx(i+1, j, k)] - Ey[get_idx(i+1, j, k+1)] - Ez[get_idx(i+1, j, k)] + Ez[get_idx(i+1, j+1, k)]);
                Hy[g] = LHy[g] * Hy[g] - RHy[g] * (Ez[get_idx(i, j+1, k)] - Ez[get_idx(i+1, j+1, k)] - Ex[get_idx(i, j+1, k)] + Ex[get_idx(i, j+1, k+1)]);
                Hz[g] = LHz[g] * Hz[g] - RHz[g] * (Ex[get_idx(i, j, k+1)] - Ex[get_idx(i, j+1, k+1)] - Ey[get_idx(i, j, k+1)] + Ey[get_idx(i+1, j, k+1)]);
            }
        }
    }

    float p = 0;
    for (int i = 1; i < grid->nx - 1; i ++) {
        for (int j = 1; j < grid->ny - 1; j ++) {
            for (int k = 1; k < grid->nz - 1; k ++) {
                auto g = get_idx(i, j, k);
                float sx = 0, sy = 0, sz = 0;
                if (g == SG) {
                    SD == 0 ? (sx = s) : SD == 1 ? (sy = s) : (sz = s);
                }
                Ex[g] = LEx[g] * Ex[g] + REx[g] * (Hy[get_idx(i, j-1, k-1)] - Hy[get_idx(i, j-1, k)] - Hz[get_idx(i, j-1, k-1)] + Hz[get_idx(i, j, k-1)] + sx);
                Ey[g] = LEy[g] * Ey[g] + REy[g] * (Hz[get_idx(i-1, j, k-1)] - Hz[get_idx(i, j, k-1)] - Hx[get_idx(i-1, j, k-1)] + Hx[get_idx(i-1, j, k)] + sy);
                Ez[g] = LEz[g] * Ez[g] + REz[g] * (Hx[get_idx(i-1, j-1, k)] - Hx[get_idx(i-1, j, k)] - Hy[get_idx(i-1, j-1, k)] + Hy[get_idx(i, j-1, k)] + sz);
                if (g == SG) {
                    p = SD == 0 ? Ex[g] : SD == 1 ? Ey[g] : Ez[g];
                }
            }
        }
    }
    return p;
}

CPUSolver::~CPUSolver() {
    delete E;
    delete H;
}

static auto Compile(Grid &grid, Port &port) {
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
    CHECK(code == 0, "compile " + dll + " failed,\ncmd: " + cmd + "\nmessage: " + utils::readFile(log));

    return new utils::DLL(dll);
}

GPUSolver::GPUSolver(Matrix &mats, float dt, vector<Port> &ports) : Solver(mats, dt, ports) {
    auto &grid = *mats.grid;
    CHECK(ports.size() == 1, "FIXME: only one port supported");
    dll = Compile(grid, ports[0]);
    FnInit = (decltype(FnInit)) dll->getProc("init_0");
    FnStep = (decltype(FnStep)) dll->getProc("step_0");
    FnQuit = (decltype(FnQuit)) dll->getProc("step_0");
    CHECK(FnInit && FnStep && FnQuit, "get function failed");
    CHECK(FnInit(coe) == 0, "solver init failed");
}

GPUSolver::~GPUSolver() {
    CHECK(FnQuit() == 0, "solver cleanup failed");
    auto path = dll->path;
    delete dll;
    filesystem::remove_all(filesystem::path(path).parent_path());
}

float GPUSolver::Step(float s) {
    return FnStep(s);
}
