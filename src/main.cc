#include "cst.h"
#include "solver.h"
#include "occ.h"
#include "plt.h"
#include "check.h"

using namespace std;

auto mesh() {
    auto grid = grid::Grid(utils::range(-5., 5., 1.), utils::range(-5., 5., 1.), utils::range(-5., 5., 1.));
    auto shape = occ::Builder::sphere(make_float3(0, 0, 0), 3.5);
    //auto shape = occ::Builder::box(make_float3(-3, -3, -3), make_float3(3, 3, 3));
    auto mesher = occ::Mesher(grid, shape);
    mesher.Save("E:\\mesh-sphere.stp");
}

auto solve() {
    auto proj = cst::Project("C:\\Projects\\fim-example.cst", "2019", false, true);
    CHECK(proj.ports.size() == 1, "Only one port supported, got " + to_string(proj.ports.size()));

    auto grid = grid::Grid(proj.xs, proj.ys, proj.zs);
    cout << "INFO: grid = " << grid.xs.size() << "x" << grid.ys.size() << "x" << grid.zs.size() << endl;
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
        //mesh();
        return 0;
    } catch (exception &e) {
        cerr << "FATAL: " << e.what() << endl;
        return -1;
    }
}
