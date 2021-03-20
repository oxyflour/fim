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
#include <boost/geometry/geometries/linestring.hpp>
#include <boost/geometry/geometries/multi_linestring.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace stl {
    namespace bg = boost::geometry;
    typedef bg::model::d2::point_xy<double> Point;
    typedef bg::model::polygon<Point> Polygon;
    typedef bg::model::multi_polygon<Polygon> MultiPolygon;
    typedef bg::model::linestring<Point> Line;
    typedef bg::model::multi_linestring<Line> MultiLine;
    typedef bg::model::segment<Point> Segment;
    typedef bg::model::box<Point> Box;

    namespace bgi = boost::geometry::index;
    typedef std::pair<Segment, double3> RTValue;
    typedef bgi::rtree<RTValue, bgi::quadratic<8, 4>> RTree;

    struct Shape {
        MultiPolygon polys;
        // Note: check performance issue with lines even you don't use it
        // MultiLine lines;
        RTree tree;
        Box bound;
    };

    struct Bound {
        double3 min = { 1e9, 1e9, 1e9 }, max = { -1e9, -1e9, -1e9 };
    };

    struct Mesh {
        std::vector<double3> vertices;
        std::vector<int3> faces;
        std::vector<double3> normals;
        std::vector<Bound> bounds;
        double3 min = { 1e9, 1e9, 1e9 }, max = { -1e9, -1e9, -1e9 };
    };

    void save(std::string file, Mesh &mesh);
    Mesh load(std::string file, double tol = 1e-6);

    struct Fragments {
        std::map<int, Shape> x, y, z;
    };

    struct Locks {
        std::mutex x, y, z;
    };

    struct Spliter {
        Spliter(grid::Grid &grid, Mesh &mesh, double tol = 1e-6);
        ~Spliter();
        Fragments fragments;
    private:
        std::vector<double> xs, ys, zs;
        int nx, ny, nz;
        double3 min, max;

        double3 *vertices, *normals;
        int3 *faces;
        Bound *bounds;
        int faceNum, vertNum;

        Locks locks;
        double tol, htol;

        Shape Slice(Mesh &mesh, double pos, int dir);
        void SliceX(grid::Grid &grid, Mesh &mesh, int i);
        void SliceY(grid::Grid &grid, Mesh &mesh, int j);
        void SliceZ(grid::Grid &grid, Mesh &mesh, int k);
    };
}

stl::Fragments& operator+=(stl::Fragments &a, stl::Fragments &b);
stl::Fragments& operator-=(stl::Fragments &a, stl::Fragments &b);

#endif
