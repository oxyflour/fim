#include "cst.h"
#include "fit.h"
#include "occ.h"

// https://github.com/lava/matplotlib-cpp
#define WITHOUT_NUMPY
// https://bugs.python.org/issue38728
#ifdef _WIN32 && _DEBUG
#undef _DEBUG
#include "matplotlibcpp.h"
#define _DEBUG
#else
#include "matplotlibcpp.h"
#endif

using namespace std;
namespace plt = matplotlibcpp;

auto solve() {
    //auto grid = Grid { range<double>(-5, 5, 1), range<double>(-5, 5, 1), range<double>(-5, 5, 1) };
    //occ::Step::save(string("E:\\test-sphere.stp"), occ::Builder::sphere(float3 { 0, 0, 0 }, 3.5));
    //occ::Step::save(string("E:\\test-sphere.stp"), occ::Builder::box(float3 { -3, -3, -3 }, float3 { 3, 3, 3 }));
    //auto mesher = occ::Mesher(grid, string("E:\\test-sphere.stp"));

    auto proj = cst::Project(string("E:\\Projects\\cst-demo\\dipole-test.cst"), string("2019"), true);
    ASSERT(proj.ports.size() == 1, "Only one port supported, got " + to_string(proj.ports.size()));
    ASSERT(proj.dt > 0, "cannot get dt from project");

    auto grid = proj.GetHexGrid();
    auto eps = proj.GetMatrix(100).data(), mue = proj.GetMatrix(101).data();
    auto mats = fit::Matrix(grid, eps, mue);
    auto port = fit::Port(grid, proj.ports[0]);
    auto solver = fit::Solver(mats, proj.dt, vector<fit::Port>({ port }));

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
        return 0;
    } catch (exception &e) {
        cerr << "FATAL: " << e.what() << endl;
        return -1;
    }
}
