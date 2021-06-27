#include <ctpl.h>

#include "cst.h"
#include "occ.h"
#include "stl.h"
#include "utils.h"

constexpr double TOL = 1e-6, EXT = 1e-2;

using namespace std;

auto test_clip() {
    stl::Polygon a, b;
    stl::bg::read_wkt("POLYGON("
        "(100 100, 100 200, 250 250, 200 100, 100 100)"
    ")", a);
    stl::bg::read_wkt("POLYGON("
        "(120 120, 120 180, 220 220, 180 120, 120 120)"
    ")", b);
    auto c = stl::MultiPolygon() + a - b;

    double v0 = 130.123456789, v1 = 190.123456789;
    stl::MultiPolygon u, v;
    {
        auto start = utils::clockNow();
        for (int i = 0; i < 10000; i ++) {
            u = stl::clip(c, 1, v0, v1);
        }
        printf("t1: %f s\n", utils::secondsSince(start));
    }
    {
        auto b = stl::MultiPolygon() + stl::Box({ 0, v0 }, { 500, v1 });
        auto start = utils::clockNow();
        for (int i = 0; i < 10000; i ++) {
            v = c * b;
        }
        printf("t2: %f s\n", utils::secondsSince(start));
    }
    std::cout << std::setprecision(12) << stl::bg::wkt(u) << endl;
    std::cout << std::setprecision(12) << stl::bg::wkt(v) << endl;
    stl::export_shape("d1.svg", u);
    stl::export_shape("d2.svg", v);
}

int main() {
    auto shape = occ::Builder::box(float3 { -1, -1, -1 }, float3 { 1, 1, 1 });
    auto geometry = occ::Mesh::triangulate(shape);
    auto mesh = stl::load(geometry.verts, geometry.faces);

    auto grid = grid::Grid(utils::range(-2., 3., 0.1), utils::range(-2., 3., 0.1), utils::range(-2., 3., 0.1));
    auto start = utils::clockNow();
    ctpl::thread_pool pool(std::thread::hardware_concurrency());
    auto splited = stl::Spliter(grid, mesh, pool);
    cout << utils::secondsSince(start) << endl;
    splited.fragments.Dump("a-test", grid);

    auto bound = splited.fragments.GetBoundary(grid, TOL, 0.1);
    bound.Dump("a-bound", grid);
    return 0;
}
