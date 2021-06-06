#include "cst.h"
#include "occ.h"
#include "stl.h"

constexpr double TOL = 1e-6, EXT = 1e-2;

int main() {
    auto shape = occ::Builder::box(float3 { -1, .3, -1 }, float3 { 1, .6, 1 });
    auto geometry = occ::Mesh::triangulate(shape);
    auto mesh = stl::load(geometry.verts, geometry.faces);

    auto grid = grid::Grid(utils::range(-2., 3., 1.), utils::range(-2., 3., 1.), utils::range(-2., 3., 1.));
    auto splited = stl::Spliter(grid, mesh);
    splited.fragments.Dump("a-test", grid);

    auto bound = splited.fragments.GetBoundary(grid, TOL, 0.1);
    bound.Dump("a-bound", grid);
}
