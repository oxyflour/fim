#include "cst.h"
#include "occ.h"
#include "stl.h"

constexpr double TOL = 1e-6, EXT = 1e-2;

int main() {
    auto center = float3 { 1.5, 1.5, 1.5 };
    auto shape = occ::Bool::fuse(
        occ::Bool::cut(
            occ::Builder::sphere(center, 1.5),
            occ::Builder::sphere(center, 1)
        ),
        occ::Builder::sphere(center, 0.5)
    );
    auto grid = grid::Grid(utils::range(-1., 4., .2), utils::range(-1., 4., .2), utils::range(-1., 5., 1.));
    auto geometry = occ::Mesh::triangulate(shape);

    auto mesh = stl::load(geometry.verts, geometry.faces);
    stl::save("a.stl", mesh);

    auto splited = stl::Spliter(grid, mesh);
    splited.fragments.Dump("a-test", grid);

    auto bound = splited.fragments.GetBoundary(grid, TOL, 0.1);
    bound.Dump("a-bound", grid);
}
