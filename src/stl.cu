#include "stl.h"
#include "ctpl_stl.h"

#include <fstream>
#include <iostream>
#include <map>

using namespace std;
using namespace stl;

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

inline __host__ __device__ double round_by(double a, double tol) {
    return round(a / tol) * tol;
}

inline __host__ __device__ double3 round_by(double3 a, double tol) {
    return double3 { round_by(a.x, tol), round_by(a.y, tol), round_by(a.z, tol) };
}

inline bool operator<(double3 a, double3 b) {
    return (a.x < b.x || (a.x == b.x && a.y < b.y) || (a.x == b.x && a.y == b.y && a.z < b.z));
}

auto swap_xy(int2 &a) {
    auto t = a.x;
    a.x = a.y;
    a.y = t;
    return a;
}

constexpr int
    DIR_X = 0,
    DIR_Y = 1,
    DIR_Z = 2;

constexpr auto
    NORM_X = double3 { 1, 0, 0 },
    NORM_Y = double3 { 0, 1, 0 },
    NORM_Z = double3 { 0, 0, 1 };

struct Joint {
    int2 e;
    double3 p;
};

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

struct Point3D {
    double3 p;
    int2 f;
};

struct Loop {
    vector<Point3D> pts;
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
                Splited { k, f.x, f.y, lerp(a, b, u.x), f.y, f.z, lerp(b, c, u.y) } :
            u.y > 0 && u.y < 1 && u.z > 0 && u.z < 1 ?
                Splited { k, f.y, f.z, lerp(b, c, u.y), f.z, f.x, lerp(c, a, u.z) } :
                Splited { k, f.z, f.x, lerp(c, a, u.z), f.x, f.y, lerp(a, b, u.x) };
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
            mesh.vertices.push_back(pt.p);
            mesh.vertices.push_back(pt.p + offset);
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
        coord[from] = edge.from.p;
        coord[to] = edge.to.p;
    }

    auto rings = vector<vector<Point3D>>();
    while (conns.size()) {
        auto ring = vector<Point3D>();
        auto &begin = *conns.begin();
        auto current = begin.first, prev = begin.second.a.i;
        while (conns.count(current)) {
            auto &pair = conns[current];
            auto face = int2 { pair.a.f, pair.b.f };
            ring.push_back(Point3D { coord[current], pair.a.i == prev ? swap_xy(face) : face });
            auto next = pair.a.i == prev ? pair.b.i : pair.a.i;
            conns.erase(current);
            prev = current;
            current = next;
        }
        rings.push_back(ring);
    }

    return rings;
}

auto is_reversed(vector<Point3D> &pts, int dir) {
    double sum = 0;
    for (int i = 0, n = pts.size(); i < n; i ++) {
        auto &a = pts[i].p, &b = pts[(i + 1) % n].p;
        if (dir == DIR_X) {
            sum += (b.y - a.y) * (b.z + a.z);
        } else if (dir == DIR_Y) {
            sum += (b.z - a.z) * (b.x + a.x);
        } else {
            sum += (b.x - a.x) * (b.y + a.y);
        }
    }
    return sum < 0;
}

auto get_loops(Mesh &mesh, Splited *splited, int num, int dir) {
    auto rings = get_rings(mesh, splited, num);
    auto norm = dir == 0 ? NORM_X : dir == 1 ? NORM_Y : NORM_Z;
    auto ret = vector<Loop>();
    for (auto &ring : rings) {
        auto &pts = ring;

        auto reversed = is_reversed(pts, dir);
        if (reversed) {
            reverse(pts.begin(), pts.end());
        }

        auto sum = 0.;
        for (int i = 0, s = pts.size(); i < s; i ++) {
            auto &a = pts[i], &b = pts[(i + 1) % s];
            if (reversed) swap_xy(a.f);
            auto &n = mesh.normals[a.f.x];
            sum += dot(cross(a.p - b.p, n), norm);
        }

        ret.push_back({ pts, sum > 0 });
    }
    return ret;
}

inline auto pt_2d(Point3D pt, int dir) {
    auto &p = pt.p;
    return dir == DIR_X ? Point { p.y, p.z } :
           dir == DIR_Y ? Point { p.z, p.x } :
                          Point { p.x, p.y };
}

auto make_poly(Loop &loop, int dir) {
    Polygon poly;
    for (auto &pt : loop.pts) {
        bg::append(poly, pt_2d(pt, dir));
    }
    bg::append(poly, pt_2d(loop.pts[0], dir));
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

// Note: structured binding not working in nvcc
map<int, MultiPolygon>& operator+=(map<int, MultiPolygon> &a, map<int, MultiPolygon> &b) {
    for (auto &pair : b) {
        auto i = pair.first;
        a[i] = a.count(i) ? a[i] + b[i] : b[i];
    }
    return a;
}

map<int, MultiPolygon>& operator-=(map<int, MultiPolygon> &a, map<int, MultiPolygon> &b) {
    auto removed = vector<int>();
    for (auto &pair : b) {
        auto i = pair.first;
        if (a.count(i)) {
            a[i] = a[i] - b[i];
            if (!a[i].size()) {
                removed.push_back(i);
            }
        }
    }
    for (auto i : removed) {
        a.erase(i);
    }
    return a;
}

Fragments& operator+=(Fragments &a, Fragments &b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

Fragments& operator-=(Fragments &a, Fragments &b) {
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

auto get_lines(MultiPolygon &polys) {
    MultiLine lines;
    lines.resize(polys.size());
    for (auto &poly : polys) {
        Line edge;
        for (auto &pt : poly.outer()) {
            bg::append(edge, pt);
        }
        for (auto &loop : poly.inners()) {
            Line edge;
            for (auto &pt : loop) {
                bg::append(edge, pt);
            }
        }
        bg::append(lines, edge);
    }
    return lines;
}

inline auto operator+(Shape &a, Shape &b) {
    Shape shape;
    shape.polys = a.polys + b.polys;
    bg::envelope(shape.polys, shape.bound);
    shape.lines = get_lines(shape.polys);
    shape.tree = a.tree;
    for (auto &item : b.tree) {
        shape.tree.insert(item);
    }
    return shape;
}

inline auto operator*(Shape &a, Box &b) {
    Shape shape;
    if (a.polys.size() && bg::intersects(a.bound, b)) {
        shape.polys = a.polys * b;
        a.tree.query(bgi::intersects(b), bgi::inserter(shape.tree));
        bg::envelope(shape.polys, shape.bound);
        MultiLine lines;
        bg::intersection(shape.lines, b, lines);
        shape.lines = lines;
    }
    return shape;
}

Shape stl::Spliter::Slice(Mesh &mesh, double pos, int dir) {
    auto bucket = malloc_device<int>(faceNum);
    auto idx = malloc_device<int>(1);
    auto splited = malloc_device<Splited>(faceNum);

    int bucketLen = 0;
    to_device(&bucketLen, 1, idx);
    kernel_group CU_DIM(256, 128) (vertices, faces, faceNum, bounds, pos, dir, bucket, idx);
    CU_ASSERT(cudaGetLastError());
    from_device(idx, 1, &bucketLen);

    Shape shape;
    if (bucketLen == 0) {
        cudaFree(bucket);
        cudaFree(idx);
        cudaFree(splited);
        return shape;
    }

    int splitedLen = 0;
    to_device(&splitedLen, 1, idx);
    kernel_split CU_DIM(256, 128) (vertices, faces, bucket, bucketLen, pos, dir, splited, idx);
    CU_ASSERT(cudaGetLastError());
    from_device(idx, 1, &splitedLen);

    auto splitedArr = new Splited[faceNum];
    from_device(splited, splitedLen, splitedArr);
    auto loops = get_loops(mesh, splitedArr, splitedLen, dir);

    delete splitedArr;
    cudaFree(bucket);
    cudaFree(idx);
    cudaFree(splited);

    for (auto &loop : loops) {
        if (!loop.hole) {
            shape.polys = shape.polys + make_poly(loop, dir);
        }
    }
    for (auto &loop : loops) {
        if (loop.hole) {
            shape.polys = shape.polys - make_poly(loop, dir);
        }
    }

    for (auto &loop : loops) {
        for (int i = 0, n = loop.pts.size(); i < n; i ++) {
            auto &a = loop.pts[i], &b = loop.pts[(i + 1) % n];
            auto s = Segment(pt_2d(a, dir), pt_2d(b, dir));
            shape.tree.insert({ s, a.f.x });
        }
    }

    shape.lines = get_lines(shape.polys);
    return shape;
}

auto check_edges(Shape &shape, Mesh &mesh) {
    printf("begin multiline\n");
    for (auto &line : shape.lines) {
        printf("begin line\n");
        for (auto &pt : line) {
            vector<RTValue> data;
            shape.tree.query(bgi::intersects(pt), back_inserter(data));
            printf("begin pt %zd\n", data.size());
            for (auto &s : data) {
                auto f = s.second;
                if (f >= 0 && f < mesh.faces.size()) {
                    auto n = mesh.normals[f];
                    printf("%f %f: %f %f %f (%d)\n", pt.x(), pt.y(), n.x, n.y, n.z, f);
                }
            }
        }
        printf("end line\n");
    }
    printf("end multiline\n");
}

void stl::Spliter::SliceX(grid::Grid &grid, Mesh &mesh, int i) {
    auto a = Slice(mesh, xs[i] - htol, DIR_X),
        b = Slice(mesh, xs[i] + htol, DIR_X),
        shape = a + b;
    for (int j = 0; j < ny - 1; j ++) {
        auto stride = shape * Box({ ys[j] - htol, min.z }, { ys[j + 1] + htol, max.z });
        for (int k = 0; k < nz - 1; k ++) {
            auto patch = stride * Box({ min.y, zs[k] - htol }, { max.y, zs[k + 1] + htol });
            if (patch.polys.size()) {
                lock_guard<mutex> guard(locks.x);
                fragments.x[grid.GetIndex(i, j, k)] = patch.polys;
            }
        }
    }
}

void stl::Spliter::SliceY(grid::Grid &grid, Mesh &mesh, int j) {
    auto a = Slice(mesh, ys[j] - htol, DIR_Y),
        b = Slice(mesh, ys[j] + htol, DIR_Y),
        shape = a + b;
    for (int k = 0; k < nz - 1; k ++) {
        auto stride = shape * Box({ zs[k] - htol, min.x }, { zs[k + 1] + htol, max.x });
        for (int i = 0; i < nx - 1; i ++) {
            auto patch = stride * Box({ min.z, xs[i] - htol }, { max.z, xs[i + 1] + htol });
            if (patch.polys.size()) {
                lock_guard<mutex> guard(locks.y);
                fragments.y[grid.GetIndex(i, j, k)] = patch.polys;
            }
        }
    }
}

void stl::Spliter::SliceZ(grid::Grid &grid, Mesh &mesh, int k) {
    auto a = Slice(mesh, zs[k] - htol, DIR_Z),
        b = Slice(mesh, zs[k] + htol, DIR_Z),
        shape = a + b;
    for (int i = 0; i < nx - 1; i ++) {
        auto stride = shape * Box({ xs[i] - htol, min.y }, { xs[i + 1] + htol, max.y });
        for (int j = 0; j < ny - 1; j ++) {
            auto patch = stride * Box({ min.x, ys[j] - htol }, { max.x, ys[j + 1] + htol });
            if (patch.polys.size()) {
                lock_guard<mutex> guard(locks.z);
                fragments.z[grid.GetIndex(i, j, k)] = patch.polys;
            }
            if (patch.lines.size() && k == 27) {
                check_edges(patch, mesh);
            }
        }
    }
}

stl::Spliter::Spliter(grid::Grid &grid, Mesh &mesh, double tol) {
    this->tol = tol;
    this->htol = tol / 2;
    min = mesh.min; max = mesh.max;
    xs = grid.xs; ys = grid.ys; zs = grid.zs;
    nx = grid.nx; ny = grid.ny; nz = grid.nz;
    for (auto &v : xs) v = round_by(v, tol);
    for (auto &v : ys) v = round_by(v, tol);
    for (auto &v : zs) v = round_by(v, tol);

    faceNum = mesh.faces.size();
    faces = to_device(mesh.faces.data(), faceNum);
    vertNum = mesh.vertices.size();
    vertices = to_device(mesh.vertices.data(), vertNum);
    normals = malloc_device<double3>(faceNum);
    bounds = malloc_device<Bound>(faceNum);

    kernel_round CU_DIM(256, 128) (vertices, vertNum, tol);
    CU_ASSERT(cudaGetLastError());

    kernel_prepare CU_DIM(256, 128) (vertices, faces, faceNum, normals, bounds);
    CU_ASSERT(cudaGetLastError());

    from_device(vertices, vertNum, mesh.vertices.data());
    from_device(normals, faceNum, mesh.normals.data());
    from_device(bounds, faceNum, mesh.bounds.data());

    ctpl::thread_pool pool(thread::hardware_concurrency());
    for (int i = 0; i < nx; i ++) {
        pool.push([&, i](int id) { SliceX(grid, mesh, i); });
    }
    for (int j = 0; j < ny; j ++) {
        pool.push([&, j](int id) { SliceY(grid, mesh, j); });
    }
    for (int k = 0; k < nz; k ++) {
        pool.push([&, k](int id) { SliceZ(grid, mesh, k); });
    }
    pool.stop(true);
}

stl::Spliter::~Spliter() {
    cudaFree(faces);
    cudaFree(vertices);
    cudaFree(normals);
    cudaFree(bounds);
}
