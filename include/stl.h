#include <string>
#include <vector>
#include <mutex>
#include <ctpl.h>

#include "utils/cuda.h"
#include "utils/clipper.h"
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
    using namespace std;
    using namespace clipper;

    namespace bg = boost::geometry;
    typedef bg::model::d2::point_xy<double> Point;
    typedef bg::model::polygon<Point> Polygon;
    typedef bg::model::multi_polygon<Polygon> MultiPolygon;
    typedef bg::model::linestring<Point> Line;
    typedef bg::model::multi_linestring<Line> MultiLine;
    typedef bg::model::segment<Point> Segment;
    typedef bg::model::box<Point> Box;

    namespace bgi = boost::geometry::index;
    struct RTValue {
        Segment segment;
        double3 normal;
    };
    typedef pair<Box, RTValue> RTPair;
    typedef bgi::rtree<RTPair, bgi::quadratic<8, 4>> RTree;

    struct Shape {
        int order;
        MultiPolygon polys;
        RTree tree;
        Box bound;
    };

    struct Bound {
        double3 min = { 1e9, 1e9, 1e9 }, max = { -1e9, -1e9, -1e9 };
    };

    struct Mesh {
        int order;
        vector<double3> vertices;
        vector<int3> faces;
        vector<double3> normals;
        vector<Bound> bounds;
        double3 min = { 1e9, 1e9, 1e9 }, max = { -1e9, -1e9, -1e9 };
    };

    void save(string file, Mesh &mesh);
    Mesh load(string file, double tol = 1e-4);
    Mesh load(vector<double3> verts, vector<int3> faces, double tol = 1e-4);

    struct Fragments {
        map<int, Shape> x, y, z;
        Fragments GetBoundary(grid::Grid &grid, double tol, double len);
        void Dump(string prefix, grid::Grid &grid);
    };

    struct Locks {
        mutex x, y, z;
    };

    struct Spliter {
        Spliter(grid::Grid &grid, Mesh &mesh, ctpl::thread_pool &pool, double tol = 1e-4, double ext = 1e-2);
        ~Spliter();
        Fragments fragments;
    private:
        vector<double> xs, ys, zs;
        int nx, ny, nz;
        double3 min, max;

        double3 *vertices, *normals;
        int3 *faces;
        Bound *bounds;
        int faceNum, vertNum;

        Locks locks;
        double tol, ext;

        Shape Slice(Mesh &mesh, double pos, int dir);
        void SliceX(grid::Grid &grid, Mesh &mesh, int i);
        void SliceY(grid::Grid &grid, Mesh &mesh, int j);
        void SliceZ(grid::Grid &grid, Mesh &mesh, int k);
    };

    // TODO
    map<int, Shape> extract_boundary(map<int, Shape> &a, grid::Grid &grid, int dir, double tol, double len = 1.);
    void export_svg(string file, map<int, Shape> &shapes, function<bool(int)> test);
    MultiPolygon clip(MultiPolygon &c, int dir, double min, double max);

    template <class T> auto export_shape(string file, T polys) {
        ofstream fn(file);
        bg::svg_mapper<Point> map(fn, 400, 400);
        map.add(polys);
        map.map(polys, "fill:blue;stroke:black;stroke-width:0.1");
    }
}

stl::Fragments& operator+= (stl::Fragments &a, stl::Fragments &b);
stl::Fragments& operator-= (stl::Fragments &a, stl::Fragments &b);

stl::Point operator+ (stl::Point &a, stl::Point &b);
stl::Point operator/ (stl::Point &a, double f);
stl::Point operator- (stl::Point &a, stl::Point &b);
stl::Point operator* (stl::Point &a, double f);

template <typename T> auto operator+(stl::MultiPolygon &shape, T &poly) {
    stl::MultiPolygon output;
    stl::bg::union_(shape, poly, output);
    return output;
}
template <typename T> auto operator-(stl::MultiPolygon &shape, T &poly) {
    stl::MultiPolygon output;
    stl::bg::difference(shape, poly, output);
    return output;
}
template <typename T> auto operator-(stl::Box &shape, T &poly) {
    stl::MultiPolygon output;
    stl::bg::difference(shape, poly, output);
    return output;
}
template <typename T> auto operator*(stl::MultiPolygon &shape, T &poly) {
    stl::MultiPolygon output;
    stl::bg::intersection(shape, poly, output);
    return output;
}

#endif
