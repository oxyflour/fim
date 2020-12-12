#include "cst.h"
#include "solver.h"
#include "occ.h"
#include "plt.h"
#include "check.h"

using namespace std;

auto mesh() {
    auto grid = grid::Grid { utils::range<double>(-5, 5, 1), utils::range<double>(-5, 5, 1), utils::range<double>(-5, 5, 1) };
    occ::Step::save(string("E:\\test-sphere.stp"), occ::Builder::sphere(make_float3(0, 0, 0), 3.5));
    //occ::Step::save(string("E:\\test-sphere.stp"), occ::Builder::box(float3(-3, -3, -3), float3(3, 3, 3)));
    auto mesher = occ::Mesher(grid, string("E:\\test-sphere.stp"));
}

auto solve() {
    auto proj = cst::Project(string("E:\\Projects\\cst-demo\\dipole-test.cst"), string("2019"), true);
    CHECK(proj.ports.size() == 1, "Only one port supported, got " + to_string(proj.ports.size()));

    auto grid = proj.GetHexGrid();
    cout << "INFO: grid = " << grid.xs.size() << "x" << grid.ys.size() << "x" << grid.zs.size() << endl;

    CHECK(proj.dt > 0, "cannot get dt from project");
    cout << "INFO: dt = " << proj.dt * 1e9 << " ns" << endl;

    auto eps = proj.GetMatrix(100), mue = proj.GetMatrix(101);
    auto mats = fit::Matrix(grid, eps, mue);
    auto port = fit::Port(grid, proj.ports[0]);
    auto solver = solver::Solver(mats, proj.dt, vector<fit::Port>({ port }));

    int steps = 3074;
    vector<float> sigt(steps), sigs(steps), sigy(steps);
    for (auto c : utils::range(steps)) {
        sigs[c] = utils::interp1(proj.excitation.x, proj.excitation.y, sigt[c] = c * proj.dt);
    }
    for (auto c : utils::range(steps)) {
        sigy[c] = solver.Step(sigs[c]);
    }

    plt::plot(sigt, sigy);
    plt::show();
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
