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

int main() {
    auto grid = Grid { range<double>(-5, 5, 1), range<double>(-5, 5, 1), range<double>(-5, 5, 1) };
    //occ::Step::save(string("E:\\test-sphere.stp"), occ::Builder::sphere(float3 { 0, 0, 0 }, 3.5));
    occ::Step::save(string("E:\\test-sphere.stp"), occ::Builder::box(float3 { -3, -3, -3 }, float3 { 3, 3, 3 }));
    auto mesher = occ::Mesher(grid, string("E:\\test-sphere.stp"));
    /*

    auto proj = cst::Project(string("E:\\Projects\\cst-demo\\dipole-test.cst"), string("2019"));
    auto grid = proj.GetHexGrid();
    auto eps = proj.GetMatrix(100).data(), mue = proj.GetMatrix(101).data();
    auto mats = fit::Matrix(grid, eps, mue);

    ASSERT(proj.ports.size() == 1, "Only one port supported");
    auto pos = proj.ports[0];
    auto port = fit::Port(grid, pos.src, pos.dst);

    // TODO
    float t = 0, steps = 1e5;
    auto solver = fit::Solver(mats, proj.dt, vector<fit::Port>({ port }));

    vector<float> sigx(steps), sigy(steps);
    for (int c = 0; c < steps; c ++, t += proj.dt) {
        auto s = utils::interp1(proj.excitation.x, proj.excitation.y, t);
        sigx[c] = t;
        sigy[c] = solver.Step(s);
    }

    cout << "ok" << endl;
     */
    return 0;
}
