#include <ctpl.h>

#include "cst.h"
#include "occ.h"
#include "stl.h"
#include "utils.h"

constexpr double TOL = 1e-6, EXT = 1e-2;

using namespace std;

template <class T> auto export_shape(string file, T polys) {
    ofstream fn(file);
    stl::bg::svg_mapper<stl::Point> map(fn, 500, 500);
    map.add(polys);
    map.map(polys, "fill:blue;stroke:black;stroke-width:0.1");
}

struct Joint {
    int i, j, k;
    stl::Point p;
    bool removed;
};
struct Ending {
    Joint *a, *b;
};
auto slice(stl::Polygon &poly, int dir, double pos, int side) {
    vector<vector<stl::Point>> rings(1 + poly.inners().size());
    auto &outer = poly.outer();
    for (int i = 0; i < outer.size() - 1; i ++) {
        rings[0].push_back(outer[i]);
    }
    auto &inners = poly.inners();
    for (int k = 0; k < inners.size(); k ++) {
        for (int i = 0; i < inners[k].size() - 1; i ++) {
            rings[k + 1].push_back(inners[k][i]);
        }
    }

    vector<Joint *> joints;
    vector<vector<Ending>> ends;
    vector<int> holes;
    for (auto &ring : rings) {
        ends.push_back(vector<Ending>(ring.size()));
    }
    for (int k = 0; k < rings.size(); k ++) {
        auto &r = rings[k];
        auto &e = ends[k];
        auto is_sliced = false;
        for (int i = 0, n = r.size(); i < n; i ++) {
            auto j = (i + 1) % n;
            auto &a = r[i], &b = r[j];
            auto test = dir == 0 ? (a.x() - pos) * (b.x() - pos)   : (a.y() - pos) * (b.y() - pos),
                fac =   dir == 0 ? (a.x() - pos) / (a.x() - b.x()) : (a.y() - pos) / (a.y() - b.x());
            if (test < 0) {
                auto p = a * (1 - fac) + b * fac;
                auto joint = new Joint { i, j, k, p, false };
                joints.push_back(joint);
                e[i].a == NULL ? (e[i].a = joint) : (e[i].b = joint);
                e[j].a == NULL ? (e[j].a = joint) : (e[j].b = joint);
                is_sliced = true;
            }
        }
        if (!is_sliced) {
            holes.push_back(k);
        }
    }
    sort(joints.begin(), joints.end(), [dir](Joint *a, Joint *b) {
        return dir == 0 ? a->p.y() < b->p.y() : a->p.x() < b->p.x();
    });

    stl::MultiPolygon ret;
    auto iter = find_if(joints.begin(), joints.end(), [](Joint *a) { return !a->removed; });
    while (iter != joints.end()) {
        auto idx = distance(joints.begin(), iter);
        auto a = joints[idx], b = joints[idx + 1];
        stl::Polygon poly;
        auto &ring = poly.outer();

        while (a && b) {
            ring.push_back(a->p);
            ring.push_back(b->p);
            a->removed = b->removed = true;

            auto &r = rings[b->k];
            auto &e = ends[b->k];
            auto inversed = ((dir == 0) ? r[b->i].x() : r[b->i].y()) * side > 0;
            auto i = inversed ? b->i : b->j;
            auto d = inversed ? -1 : 1;
            while (true) {
                ring.push_back(r[i]);
                auto joint = e[i].a == b ? e[i].b : e[i].a;
                if (joint && joint != b) {
                    if (joint->removed) {
                        a = b = NULL;
                    } else {
                        a = joint;
                        b = *(find(joints.begin(), joints.end(), joint) + 1);
                    }
                    break;
                } else {
                    i = (i + d + r.size()) % r.size();
                }
            }
        }
        stl::bg::correct(poly);
        ret.push_back(poly);

        iter = find_if(joints.begin(), joints.end(), [](Joint *a) { return !a->removed; });
    }

    for (auto k : holes) {
        auto p = stl::Polygon();
        auto &outer = p.outer();
        outer.insert(outer.begin(), inners[k - 1].begin(), inners[k - 1].end());
        stl::bg::correct(p);
        ret = ret - p;
    }
    return ret;
}

auto clip(stl::MultiPolygon &c, int dir, double min, double max) {
    stl::MultiPolygon e;
    for (auto &p : c) {
        e = e + slice(p, dir, min, -1);
    }
    stl::MultiPolygon f;
    for (auto &p : e) {
        f = f + slice(p, dir, max, 1);
    }
    return f;
}

int main() {
    /*
    auto shape = occ::Builder::box(float3 { -1, .3, -1 }, float3 { 1, .6, 1 });
    auto geometry = occ::Mesh::triangulate(shape);
    auto mesh = stl::load(geometry.verts, geometry.faces);

    auto grid = grid::Grid(utils::range(-2., 3., 1.), utils::range(-2., 3., 1.), utils::range(-2., 3., 1.));
    ctpl::thread_pool pool(std::thread::hardware_concurrency());
    auto splited = stl::Spliter(grid, mesh, pool);
    splited.fragments.Dump("a-test", grid);

    auto bound = splited.fragments.GetBoundary(grid, TOL, 0.1);
    bound.Dump("a-bound", grid);
     */

    stl::Polygon a, b;
    stl::bg::read_wkt("POLYGON("
        "(100 100, 100 200, 200 200, 200 100, 100 100)"
    ")", a);
    stl::bg::read_wkt("POLYGON("
        "(120 120, 120 180, 180 180, 180 120, 120 120)"
    ")", b);
    auto c = stl::MultiPolygon() + a - b;

    {
        auto start = utils::clockNow();
        for (int i = 0; i < 10000; i ++) {
            auto d = clip(c, 0, 110, 190);
        }
        printf("t1: %f s\n", utils::secondsSince(start));
    }
    {
        auto b = stl::Box({ 110, 0 }, { 190, 500 });
        auto start = utils::clockNow();
        for (int i = 0; i < 10000; i ++) {
            auto d = c * b;
        }
        printf("t2: %f s\n", utils::secondsSince(start));
    }

    //export_shape("out.svg", clip(c, 0, 110, 190));
    return 0;
}
