#include <ctpl_stl.h>

#include "utils/check.h"
#include "utils/plt.h"

#include "cst.h"
#include "solver.h"
#include "occ.h"
#include "stl.h"

using namespace std;

auto solve() {
    auto proj = cst::Project("C:\\Projects\\fim-cylinder.cst", "2019", true, true);
    auto grid = grid::Grid(proj.xs, proj.ys, proj.zs);
    cout << "INFO: grid = " << grid.xs.size() << "x" << grid.ys.size() << "x" << grid.zs.size() << endl;

    auto mats = map<string, stl::Fragments>();
    mutex lock;

    ctpl::thread_pool pool;
    for (auto &solid : proj.solids) {
        pool.push([&] (int id) {
            auto mesh = stl::load(solid.stl);
            auto fragments = stl::Spliter(grid, mesh).fragments;
            {
                lock_guard guard(lock);
                if (mats.count(solid.material)) {
                    mats[solid.material] += fragments;
                } else {
                    mats[solid.material] = fragments;
                }
            }
        });
    }
    pool.stop(true);

/*

    CHECK(proj.dt > 0, "cannot get dt from project");
    cout << "INFO: dt = " << proj.dt * 1e9 << " ns" << endl;

    auto eps = proj.GetMatrix(100), mue = proj.GetMatrix(101);
    auto mats = fit::Matrix(grid, eps, mue);
    auto port = fit::Port(grid, proj.ports[0]);
    auto cpuSolver = solver::CPUSolver(mats, proj.dt, vector({ port }));
    auto gpuSolver = solver::GPUSolver(mats, proj.dt, vector({ port }));

    int steps = 3074;
    vector<float> sigt(steps), sigs(steps), sigy1(steps), sigy2(steps);
    for (auto c : utils::range(steps)) {
        auto t = sigt[c] = c * proj.dt;
        sigs[c] = utils::interp1(proj.excitation.x, proj.excitation.y, t) * port.power;
    }

    for (auto c : utils::range(steps)) {
        sigy1[c] = cpuSolver.Step(sigs[c]);
        sigy2[c] = gpuSolver.Step(sigs[c]);
        utils::outputProgress((c + 1) * 1. / steps);
    }

    plt::plot(sigt, sigy1, {{ "label", "cpu" }, { "ls", "-" }});
    plt::plot(sigt, sigy2, {{ "label", "gpu" }, { "ls", "--" }});
    plt::xlabel("time (s)");
    plt::ylabel("voltage (V)");
    plt::legend();
    plt::show();
 */
}

int main() {
    try {
        solve();
        return 0;
    } catch (exception &e) {
        cerr << "FATAL: " << e.what() << endl;
        return -1;
    }
}
