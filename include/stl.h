#include <string>
#include <vector>
#include <mutex>

#include "utils/cuda.h"
#include "grid.h"

#ifndef STL_H
#define STL_H

#undef __CUDACC__
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

namespace stl {
    namespace bg = boost::geometry;
    typedef bg::model::d2::point_xy<double> Point;
    typedef bg::model::polygon<Point> Polygon;
    typedef bg::model::multi_polygon<Polygon> MultiPolygon;
    typedef bg::model::box<Point> Box;

    struct Bound {
        double3 min = { 1e9, 1e9, 1e9 }, max = { -1e9, -1e9, -1e9 };
    };

    struct Joint {
        int2 e;
        double u;
    };

    struct Mesh {
        std::vector<double3> vertices;
        std::vector<int3> faces;
        std::vector<double3> normals;
        std::vector<Bound> bounds;
        double3 min = { 1e9, 1e9, 1e9 }, max = { -1e9, -1e9, -1e9 };
        double3 get(Joint &item);
    };

    void save(std::string file, Mesh &mesh);
    Mesh load(std::string file, double tol = 1e-6);

    struct Mesher {
        Mesher(grid::Grid &grid, Mesh &mesh);
        ~Mesher();
    private:
        void SplitX(Mesh &mesh, int i);
        std::vector<double> xs, ys, zs;
        int nx, ny, nz;
        double3 min, max;

        double3 *vertices, *normals;
        int3 *faces;
        Bound *bounds;
        int faceNum, vertNum;

        std::mutex lock;
        std::map<int, MultiPolygon> sx, sy, sz;

        double tol;
    };
}

#endif
