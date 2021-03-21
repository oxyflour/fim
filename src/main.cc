#include <ctpl_stl.h>
#include <regex>

#include "utils/check.h"
#include "utils/plt.h"

#include "cst.h"
#include "solver.h"
#include "occ.h"
#include "stl.h"

using namespace std;

auto is_dielec(string name) {
    return name != "PEC";
}

auto mat_priority = map<string, int>();

auto export_svg(string file, map<int, stl::Shape> &shapes, function<bool(int)> test) {
    auto list = vector<stl::Shape *>();
    for (auto &[n, s] : shapes) {
        if (test(n)) {
            list.push_back(&s);
        }
    }
    if (list.size()) {
        ofstream fn(file);
        stl::bg::svg_mapper<stl::Point> map(fn, 500, 500);
        for (auto ptr : list) {
            map.add(ptr->polys);
        }
        for (auto ptr : list) {
            map.map(ptr->polys, "fill:blue;stroke:black;stroke-width:0.1");
        }
    }
}

constexpr double TOL = 1e-6, EXT = 1e-2;

auto solve() {
    auto proj = cst::Project("C:\\Projects\\fim-example.cst", "2019", true, true);
    auto grid = grid::Grid(proj.xs, proj.ys, proj.zs);
    cout << "INFO: grid = " << grid.xs.size() << "x" << grid.ys.size() << "x" << grid.zs.size() << endl;

    auto start = utils::clockNow();
    auto mats = map<string, stl::Fragments>();
    for (auto &solid : proj.solids) {
        auto mesh = stl::load(solid.stl, TOL);
        // used when merging all metals
        mesh.order = mat_priority[solid.material];

        auto fragments = stl::Spliter(grid, mesh, TOL, EXT).fragments;
        if (mats.count(solid.material)) {
            mats[solid.material] += fragments;
        } else {
            mats[solid.material] = fragments;
        }
    }
    cout << "PERF: meshed in " << utils::secondsSince(start) << " s" << endl;

    auto dielecs = vector<string>();
    auto allMetal = stl::Fragments();
    for (auto &[mat, frags] : mats) {
        if (is_dielec(mat)) {
            dielecs.push_back(mat);
        } else {
            allMetal += frags;
        }
    }

    sort(dielecs.begin(), dielecs.end(), [](string a, string b) {
        return mat_priority[a] - mat_priority[b];
    });

    for (int i = 0; i < dielecs.size(); i ++) {
        auto &a = mats[dielecs[i]];
        for (int j = i + 1; j < dielecs.size(); j ++) {
            a -= mats[dielecs[j]];
        }
        a -= allMetal;
    }

    auto bound = stl::extract_boundary(allMetal.z, grid, 2, EXT);
    for (int k = 0; k < grid.nz; k ++) {
        export_svg("ALL_METAL.z-" + to_string(k) + "-" + to_string(grid.zs[k]) + ".svg",
            bound, [&, k](int n) { return grid.GetIndex(n).z == k; });
    }

    // debug
    mats["ALL_METAL"] = allMetal;
    for (auto &[mat, frags] : mats) {
        for (int i = 0; i < grid.nx; i ++) {
            export_svg("mat-" + mat + ".x-" + to_string(i) + "-" + to_string(grid.xs[i]) + ".svg",
                frags.x, [&, i](int n) { return grid.GetIndex(n).x == i; });
        }
        for (int j = 0; j < grid.ny; j ++) {
            export_svg("mat-" + mat + ".y-" + to_string(j) + "-" + to_string(grid.ys[j]) + ".svg",
                frags.y, [&, j](int n) { return grid.GetIndex(n).y == j; });
        }
        for (int k = 0; k < grid.nz; k ++) {
            export_svg("mat-" + mat + ".z-" + to_string(k) + "-" + to_string(grid.zs[k]) + ".svg",
                frags.z, [&, k](int n) { return grid.GetIndex(n).z == k; });
        }
    }

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
