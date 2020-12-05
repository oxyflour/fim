#include "cst.h"
#include "fit.h"
#include "occ.h"

using namespace std;

template <typename T> auto range(T from, T to, T delta) {
    auto ret = vector<T>();
    for (auto val = from; val < to; val += delta) {
        ret.push_back(val);
    }
    return ret;
}

int wmain() {
    //auto grid = Grid { range<double>(-5, 5, 1), range<double>(-5, 5, 1), range<double>(-5, 5, 1) };
    //occ::Step::save(string("E:\\test-sphere.stp"), occ::Builder::sphere(float3 { 0, 0, 0 }, 3.5));
    //occ::Step::save(string("E:\\test-sphere.stp"), occ::Builder::box(float3 { -3, -3, -3 }, float3 { 3, 3, 3 }));
    //auto mesher = occ::Mesher(grid, string("E:\\test-sphere.stp"));

    auto proj = cst::Project(string("E:\\Projects\\cst-demo\\dipole-test.cst"), string("2019"));
    ASSERT(proj.ports.size() == 1, "Only one port supported");
    ASSERT(proj.dt > 0, "cannot get dt from project");

    auto grid = proj.GetHexGrid();
    auto eps = proj.GetMatrix(100).data(), mue = proj.GetMatrix(101).data();
    auto mats = fit::Matrix(grid, eps, mue);
    auto pos = proj.ports[0];
    auto port = fit::Port(grid, pos.src, pos.dst);

    // TODO
    auto solver = fit::Solver(mats, proj.dt, vector<fit::Port>({ port }));

    int steps = 1074;
    vector<float> sigt(steps), sigs(steps), sigy(steps);

    double t = 0;
    for (int c = 0; c < steps; c ++, t += proj.dt) {
        sigt[c] = t;
        sigs[c] = utils::interp1(proj.excitation.x, proj.excitation.y, t);
    }
    for (int c = 0; c < steps; c ++) {
        sigy[c] = solver.Step(sigs[c]);
    }
    for (int c = 0; c < steps; c ++) {
        if (c % 10 == 0) {
            cout << sigt[c] << ", " << sigs[c] << ", " << sigy[c] << endl;
        }
    }

    return 0;
}

int main() {
    try {
        return wmain();
    } catch (exception &e) {
        cerr << e.what() << endl;
        return -1;
    }
}
