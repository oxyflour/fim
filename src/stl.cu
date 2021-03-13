#include "stl.h"

#include <fstream>
#include <iostream>
#include <map>

#undef __CUDACC__
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

using namespace std;
using namespace stl;

namespace bg = boost::geometry;
typedef bg::model::d2::point_xy<double> Point;
typedef bg::model::polygon<Point> Polygon;
typedef bg::model::multi_polygon<Polygon> MultiPolygon;
typedef bg::model::box<Point> Box;

inline __host__ __device__ double fmin(double a, double b) {
    return a < b ? a : b;
}

inline __host__ __device__ double fmax(double a, double b) {
    return a > b ? a : b;
}

inline __host__ __device__ double3 operator+(double3 a, double3 b) {
    return double3 { a.x + b.x, a.y + b.y, a.z + b.z };
}

inline __host__ __device__ double3 operator-(double3 a, double3 b) {
    return double3 { a.x - b.x, a.y - b.y, a.z - b.z };
}

inline __host__ __device__ double3 operator*(double3 a, double b) {
    return double3 { a.x * b, a.y * b, a.z * b };
}

inline __host__ __device__ double3 operator/(double3 a, double b) {
    return double3 { a.x / b, a.y / b, a.z / b };
}

inline __host__ __device__ double3 lerp(double3 a, double3 b, double f) {
    return a * (1 - f) + b * f;
}

inline __host__ __device__ double dot(double3 a, double3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline __host__ __device__ double length(double3 a) {
    return sqrt(dot(a, a));
}

inline __host__ __device__ double3 cross(double3 a, double3 b) { 
    return double3 { a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x }; 
}

inline __host__ __device__ double3 normalize(double3 v) {
    return v / length(v);
}

inline __host__ __device__ double3 fmin(double3 a, double3 b) { 
    return double3 { fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z) }; 
}

inline __host__ __device__ double3 fmax(double3 a, double3 b) { 
    return double3 { fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z) }; 
}

inline __host__ __device__ double3 round_by(double3 a, double tol) {
    return double3 { round(a.x / tol) * tol, round(a.y / tol) * tol, round(a.z / tol) * tol };
}

inline bool operator<(double3 a, double3 b) {
    return (a.x < b.x || (a.x == b.x && a.y < b.y) || (a.x == b.x && a.y == b.y && a.z < b.z));
}

constexpr int
    DIR_X = 0,
    DIR_Y = 1,
    DIR_Z = 2;

struct Splited {
    int f;
    Joint from, to;
};

struct Edge {
    int i, f;
};

struct Conn {
    Edge a, b;
};

struct Ring {
    vector<double3> pts;
    int2 f;
};

struct Loop {
    vector<double3> pts;
    bool hole;
};

__global__ void kernel_round(double3 *vertices, int vertNum, double tol) {
    for (int i = cuIdx(x); i < vertNum; i += cuDim(x)) {
        vertices[i] = round_by(vertices[i], tol);
    }
}

__global__ void kernel_prepare(
    double3 *vertices, int3 *faces, int faceNum,
    double3 *normals, Bound *bounds) {
    for (int i = cuIdx(x); i < faceNum; i += cuDim(x)) {
        auto &f = faces[i];
        auto &a = vertices[f.x], &b = vertices[f.y], &c = vertices[f.z];
        normals[i] = normalize(cross(a - b, b - c));
        bounds[i].min = fmin(fmin(a, b), c);
        bounds[i].max = fmax(fmax(a, b), c);
    }
}

__global__ void kernel_split(
    double3 *vertices, int3 *faces, int* bucket, int binLen, double val, int dir,
    Splited *out, int *idx) {
    for (int i = cuIdx(x); i < binLen; i += cuDim(x)) {
        auto k = bucket[i];
        auto &f = faces[k];
        auto &a = vertices[f.x], &b = vertices[f.y], &c = vertices[f.z];
        double3 m;
        if (dir == DIR_X) {
            m.x = a.x; m.y = b.x; m.z = c.x;
        } else if (dir == DIR_Y) {
            m.x = a.y; m.y = b.y; m.z = c.y;
        } else if (dir == DIR_Z) {
            m.x = a.z; m.y = b.z; m.z = c.z;
        }
        double3 u;
        u.x = (m.x - val) / (m.x - m.y);
        u.y = (m.y - val) / (m.y - m.z);
        u.z = (m.z - val) / (m.z - m.x);
        out[atomicAdd(idx, 1)] =
            u.x > 0 && u.x < 1 && u.y > 0 && u.y < 1 ?
                Splited { k, f.x, f.y, u.x, f.y, f.z, u.y } :
            u.y > 0 && u.y < 1 && u.z > 0 && u.z < 1 ?
                Splited { k, f.y, f.z, u.y, f.z, f.x, u.z } :
                Splited { k, f.z, f.x, u.z, f.x, f.y, u.x };
    }
}

__global__ void kernel_group(
    double3 *vertices, int3 *faces, int faceNum, Bound *bounds, double val, int dir,
    int *out, int *idx) {
    for (int i = cuIdx(x); i < faceNum; i += cuDim(x)) {
        auto &b = bounds[i];
        auto min = dir == DIR_X ? b.min.x : dir == DIR_Y ? b.min.y : b.min.z,
            max = dir == DIR_X ? b.max.x : dir == DIR_Y ? b.max.y : b.max.z;
        if (min < val && val < max) {
            out[atomicAdd(idx, 1)] = i;
        }
    }
}

void stl::save(string file, Mesh &mesh) {
    ofstream fn(file);
    fn << "solid component" << endl;
    for (int i = 0; i < mesh.faces.size(); i ++) {
        auto f = mesh.faces[i];
        auto a = mesh.vertices[f.x], b = mesh.vertices[f.y], c = mesh.vertices[f.z];
        if (mesh.normals.size() > i) {
            auto n = mesh.normals[i];
            fn << "  facet normal " << n.x << " " << n.y << " " << n.z << endl;
        } else {
            fn << "  facet" << endl;
        }
        fn << "    outer loop" << endl;
        fn << "      vertex " << a.x << " " << a.y << " " << a.z << endl;
        fn << "      vertex " << b.x << " " << b.y << " " << b.z << endl;
        fn << "      vertex " << c.x << " " << c.y << " " << c.z << endl;
        fn << "    endloop" << endl;
        fn << "  endfacet" << endl;
    }
    fn << "endsolid component" << endl;
    fn.close();
}

void save_loop(string file, Loop &loop, double3 offset) {
    auto mesh = Mesh();
    auto num = (int) loop.pts.size();
    if (num > 2) {
        auto idx = [&](int i) { return (i + num * 2) % (num * 2); };
        for (auto &pt : loop.pts) {
            auto start = (int) mesh.faces.size();
            mesh.vertices.push_back(pt);
            mesh.vertices.push_back(pt + offset);
            mesh.faces.push_back(int3 { start, idx(start + 1), idx(start - 1) });
            mesh.faces.push_back(int3 { start, idx(start - 1), idx(start - 2) });
        }
    } else {
        cout << "WARN: only " << num << " points in loop, ignoring" << endl;
    }
    save(file, mesh);
}

Mesh stl::load(string file, double tol) {
    Mesh mesh;
    map<double3, int> indices;

    ifstream fn(file);
    string line;
    while (getline(fn, line)) {
        if (line.find("face") == string::npos || line.find("endfacet") != string::npos) {
            continue;
        }
        double3 norm;
        auto start = line.c_str() + line.find_first_of("face normal");
        sscanf(start, "face normal %lf %lf %lf", &norm.x, &norm.y, &norm.z);

        // skip "outer loop"
        getline(fn, line);

        int idx[3];
        for (int i = 0; i < 3; i ++) {
            double3 pt;
            getline(fn, line);
            auto start = line.c_str() + line.find_first_of("vertex");
            sscanf(start, "vertex %lf %lf %lf", &pt.x, &pt.y, &pt.z);
            auto rounded = round_by(pt, tol);
            if (!indices.count(rounded)) {
                indices[rounded] = mesh.vertices.size();
                mesh.vertices.push_back(rounded);
            }
            idx[i] = indices[rounded];
            mesh.min = fmin(mesh.min, pt);
            mesh.max = fmax(mesh.max, pt);
        }

        if (idx[0] != idx[1] && idx[1] != idx[2] && idx[2] != idx[0]) {
            mesh.faces.push_back(int3 { idx[0], idx[1], idx[2] });
            mesh.normals.push_back(norm);
        }
    }
    fn.close();

    mesh.normals.resize(mesh.faces.size());
    mesh.bounds.resize(mesh.faces.size());

    return mesh;
}

double3 stl::Mesh::get(Joint &joint) {
    auto a = vertices[joint.e.x], b = vertices[joint.e.y];
    return lerp(a, b, joint.u);
}

auto idx_of(int2 &i, int n) {
    return i.x < i.y ? i.x * n + i.y : i.y * n + i.x;
}

auto get_rings(Mesh &mesh, Splited *splited, int num) {
    auto maxVert = mesh.vertices.size() + 1;
    auto conns = map<int, Conn>();
    auto coord = map<int, double3>();
    for (int i = 0; i < num; i ++) {
        auto &edge = splited[i];
        auto f = edge.f;
        auto from = idx_of(edge.from.e, maxVert),
            to = idx_of(edge.to.e, maxVert);
        conns.count(from) ? (conns[from].b = { to, f }) : (conns[from].a = { to, f });
        conns.count(to) ? (conns[to].b = { from, f }) : (conns[to].a = { from, f });
        coord[from] = mesh.get(edge.from);
        coord[to] = mesh.get(edge.to);
    }

    auto loops = vector<Ring>();
    while (conns.size()) {
        auto loop = Ring();
        auto &begin = *conns.begin();
        auto current = begin.first, prev = begin.second.a.i;
        loop.f = { begin.second.a.f, begin.second.b.f };
        while (conns.count(current)) {
            loop.pts.push_back(coord[current]);
            auto &pair = conns[current];
            auto next = pair.a.i == prev ? pair.b.i : pair.a.i;
            conns.erase(current);
            prev = current;
            current = next;
        }
        loops.push_back(loop);
    }

    return loops;
}

auto is_reversed(vector<double3> &pts, int dir) {
    auto sum = 0;
    for (int i = 0; i < pts.size() - 1; i ++) {
        auto &a = pts[i], &b = pts[i + 1];
        if (dir == DIR_X) {
            sum += (b.y - a.y) * (b.z + a.z);
        }
    }
    return sum < 0;
}

auto get_loops(Mesh &mesh, Splited *splited, int num, int dir) {
    auto loops = get_rings(mesh, splited, num);
    auto norm =
        dir == 0 ? double3 { 1, 0, 0 } :
        dir == 1 ? double3 { 0, 1, 0 } : double3 { 0, 0, 1 };
    auto ret = vector<Loop>();
    for (auto &loop : loops) {
        auto &pts = loop.pts;
        auto face = loop.f.x;

        auto reversed = is_reversed(pts, DIR_X);
        if (reversed) {
            reverse(pts.begin(), pts.end());
            face = loop.f.y;
        }

        auto hole = false;
        if (pts.size() >= 2) {
            auto &a = pts[0], &b = pts[1];
            hole = dot(cross(a - b, mesh.normals[face]), norm) > 0;
        }

        ret.push_back({ loop.pts, hole });
    }
    return ret;
}

auto make_poly(Loop &loop) {
    Polygon poly;
    for (auto &pt : loop.pts) {
        bg::append(poly, Point { pt.y * 100, pt.z * 100 });
    }
    bg::correct(poly);
    return poly;
}

template <typename T> auto operator+(MultiPolygon &shape, T &poly) {
    MultiPolygon input;
    bg::union_(shape, poly, input);
    return input;
}

template <typename T> auto operator-(MultiPolygon &shape, T &poly) {
    MultiPolygon input;
    bg::difference(shape, poly, input);
    return input;
}

template <typename T> auto operator*(MultiPolygon &shape, T &poly) {
    MultiPolygon input;
    bg::intersection(shape, poly, input);
    return input;
}

void stl::Mesher::SplitX(Mesh &mesh, int i, double x) {
    auto bucket = malloc_device<int>(faceNum);
    auto idx = malloc_device<int>(1);
    auto splited = malloc_device<Splited>(faceNum);

    int bucketLen = 0;
    to_device(&bucketLen, 1, idx);
    kernel_group CU_DIM(256, 128) (vertices, faces, faceNum, bounds, x, DIR_X, bucket, idx);
    CU_ASSERT(cudaGetLastError());
    from_device(idx, 1, &bucketLen);
    if (bucketLen == 0) {
        return;
    }

    int splitedLen = 0;
    to_device(&splitedLen, 1, idx);
    kernel_split CU_DIM(256, 128) (vertices, faces, bucket, bucketLen, x, DIR_X, splited, idx);
    CU_ASSERT(cudaGetLastError());

    auto splitedArr = new Splited[faceNum];
    from_device(splited, bucketLen, splitedArr);
    auto loops = get_loops(mesh, splitedArr, bucketLen, DIR_X);
    /*
    for (int m = 0; m < loops.size(); m ++) {
        save_loop("loop" + to_string(i) + "-n" + to_string(m) + ".stl", loops[m], double3 { 0.01, 0, 0 });
    }
     */
    MultiPolygon shape;
    for (auto &loop : loops) {
        if (!loop.hole) {
            shape = shape + make_poly(loop);
        }
    }
    for (auto &loop : loops) {
        if (loop.hole) {
            shape = shape - make_poly(loop);
        }
    }

    //shape = shape * Box(Point { 0, 0 }, Point { 300, 300 });

    ofstream fn("svg" + to_string(i) + ".svg");
    bg::svg_mapper<Point> mapper(fn, 500, 500);
    mapper.add(shape);
    mapper.map(shape, "fill:black;stroke:blue");

    cout << "loop " << i << " at " << x << " loop " << loops.size() << endl;
}

stl::Mesher::Mesher(grid::Grid &grid, Mesh &mesh) {
    xs = grid.xs; ys = grid.ys; zs = grid.zs;

    faces = to_device(mesh.faces.data(), mesh.faces.size());
    faceNum = mesh.faces.size();
    vertices = to_device(mesh.vertices.data(), mesh.vertices.size());
    vertNum = mesh.vertices.size();
    normals = malloc_device<double3>(faceNum);
    bounds = malloc_device<Bound>(faceNum);

    tol = 1e-5;
    kernel_round CU_DIM(256, 128) (vertices, vertNum, tol);
    CU_ASSERT(cudaGetLastError());

    kernel_prepare CU_DIM(256, 128) (vertices, faces, faceNum, normals, bounds);
    CU_ASSERT(cudaGetLastError());

    from_device(vertices, vertNum, mesh.vertices.data());
    from_device(normals, faceNum, mesh.normals.data());
    from_device(bounds, faceNum, mesh.bounds.data());

    for (int i = 0; i < grid.xs.size(); i ++) {
        auto x = (round(grid.xs[i] / tol) + 0.5) * tol;
        SplitX(mesh, i, x);
    }
}
