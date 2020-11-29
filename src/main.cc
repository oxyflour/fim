#include "cst.h"
#include "fit.h"
#include "occ.h"

using namespace std;

int main() {
    auto proj = cst::Project(string("E:\\Projects\\cst-demo\\dipole-test.cst"), string("2019"));
    auto grid = proj.GetHexGrid();
    auto eps = proj.GetMatrix(100).data(), mue = proj.GetMatrix(101).data();
    auto mats = fit::Matrix(grid, eps, mue);

    //auto mesher = occ::Mesher(grid, proj.units.geometry, string("E:\\Projects\\cst-demo\\dipole-test-sphere.stp"));
    //return 0;

    ASSERT(proj.ports.size() == 1, "Only one port supported");
    auto pos = proj.ports[0];
    auto port = fit::Port(grid, pos.src, pos.dst);

    // TODO
    float t = 0, dt = 1e-12, steps = 1e5;
    auto solver = fit::Solver(mats, dt, vector<fit::Port>({ port }));

    vector<float> sigx(steps), sigy(steps);
    for (int c = 0; c < steps; c ++, t += dt) {
        auto s = utils::interp1(proj.excitation.x, proj.excitation.y, t);
        sigx[c] = t;
        sigy[c] = solver.Step(s);
    }

    cout << "ok" << endl;
    return 0;
}
