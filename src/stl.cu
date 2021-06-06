#include "stl.h"
#include "ctpl_stl.h"
#include "utils.h"

#include <fstream>
#include <iostream>
#include <map>
#include <set>

using namespace std;
using namespace stl;

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
    double area;
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

Mesh stl::load(vector<double3> verts, vector<int3> faces, double tol) {
    Mesh mesh;
    map<double3, int> indices;
    map<int, int> pts;

    for (int i = 0; i < verts.size(); i ++) {
        auto pt = verts[i];
        auto rounded = round_by(pt, tol);
        if (!indices.count(rounded)) {
            indices[rounded] = mesh.vertices.size();
            mesh.vertices.push_back(rounded);
        }
        pts[i] = indices[rounded];
        mesh.min = fmin(mesh.min, pt);
        mesh.max = fmax(mesh.max, pt);
    }
    for (auto fv : faces) {
        auto face = int3 { pts[fv.x], pts[fv.y], pts[fv.z] };
        if (face.x != face.y && face.y != face.z && face.z != face.x) {
            mesh.faces.push_back(face);
        }
    }
    mesh.normals.resize(mesh.faces.size());
    mesh.bounds.resize(mesh.faces.size());

    return mesh;
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

        int3 face;
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
            auto idx = indices[rounded];
            i == 0 ? (face.x = idx) :
            i == 1 ? (face.y = idx) :
                     (face.z = idx);
            mesh.min = fmin(mesh.min, pt);
            mesh.max = fmax(mesh.max, pt);
        }

        if (face.x != face.y && face.y != face.z && face.z != face.x) {
            mesh.faces.push_back(face);
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

auto get_area(vector<Point3D> &pts, int dir) {
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
    return sum;
}

auto get_loops(Mesh &mesh, Splited *splited, int num, int dir) {
    auto rings = get_rings(mesh, splited, num);
    auto norm = dir == 0 ? NORM_X : dir == 1 ? NORM_Y : NORM_Z;
    auto ret = vector<Loop>();
    for (auto &ring : rings) {
        auto &pts = ring;

        auto area = get_area(pts, dir);
        if (area < 0) {
            reverse(pts.begin(), pts.end());
        }

        auto sum = 0.;
        for (int i = 0, s = pts.size(); i < s; i ++) {
            auto &a = pts[i], &b = pts[(i + 1) % s];
            if (area < 0) swap_xy(a.f);
            auto &n = mesh.normals[a.f.x];
            sum += dot(cross(a.p - b.p, n), norm);
        }

        ret.push_back({ pts, area, sum > 0 });
    }
    return ret;
}

inline auto pt_2d(Point3D pt, int dir) {
    auto &p = pt.p;
    return dir == DIR_X ? Point { p.y, p.z } :
           dir == DIR_Y ? Point { p.z, p.x } :
                          Point { p.x, p.y };
}

auto make_polys(Loop &loop, int dir) {
    stl::Polygon poly;
    for (auto &pt : loop.pts) {
        bg::append(poly, pt_2d(pt, dir));
    }
    bg::append(poly, pt_2d(loop.pts[0], dir));
    return poly;
}

stl::Point operator+ (stl::Point &a, stl::Point &b) {
    return stl::Point { a.x() + b.x(), a.y() + b.y() };
}

stl::Point operator- (stl::Point &a, stl::Point &b) {
    return stl::Point { a.x() - b.x(), a.y() - b.y() };
}

stl::Point operator* (stl::Point &a, double f) {
    return stl::Point { a.x() * f, a.y() * f };
}

stl::Point operator/ (stl::Point &a, double f) {
    return stl::Point { a.x() / f, a.y() / f };
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

template <typename T> auto operator-(Box &shape, T &poly) {
    MultiPolygon input;
    bg::difference(shape, poly, input);
    return input;
}

template <typename T> auto operator*(MultiPolygon &shape, T &poly) {
    MultiPolygon input;
    bg::intersection(shape, poly, input);
    return input;
}

inline auto operator+(RTree &a, RTree &b) {
    RTree tree = a;
    for (auto &item : b) {
        tree.insert(item);
    }
    return tree;
}

inline auto operator+(Shape &a, Shape &b) {
    Shape shape;
    shape.polys = a.polys + b.polys;
    bg::envelope(shape.polys, shape.bound);
    shape.tree = a.tree + b.tree;
    shape.order = max(a.order, b.order);
    return shape;
}

inline auto operator-(Shape &a, Shape &b) {
    Shape shape;
    shape.polys = a.polys - b.polys;
    bg::envelope(shape.polys, shape.bound);
    shape.tree = a.tree + b.tree;
    shape.order = max(a.order, b.order);
    return shape;
}

inline auto operator*(Shape &a, Box &b) {
    Shape shape;
    if (a.polys.size() && bg::intersects(a.bound, b)) {
        shape.polys = a.polys * b;
        a.tree.query(bgi::intersects(b), bgi::inserter(shape.tree));
        bg::envelope(shape.polys, shape.bound);
    }
    shape.order = a.order;
    return shape;
}

// Note: structured binding not working in nvcc
map<int, Shape>& operator+=(map<int, Shape> &a, map<int, Shape> &b) {
    for (auto &pair : b) {
        auto i = pair.first;
        a[i] = a.count(i) ? a[i] + b[i] : b[i];
    }
    return a;
}

map<int, Shape>& operator-=(map<int, Shape> &a, map<int, Shape> &b) {
    auto removed = vector<int>();
    for (auto &pair : b) {
        auto i = pair.first;
        if (a.count(i)) {
            a[i] = a[i] - b[i];
            if (!a[i].polys.size()) {
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

inline auto x_or_y_inside(Point &pt, Point &min, Point &max, double ext = 0) {
    auto x = pt.x(), y = pt.y();
    return (x > min.x() + ext && x < max.x() - ext) || (y > min.y() + ext && y < max.y() - ext);
}

inline auto get_normal(Shape &shape, Point &pt, double tol) {
    vector<RTValue> norms;
    shape.tree.query(bgi::intersects(pt), back_inserter(norms));
    if (!norms.size()) {
        vector <RTValue> nearest;
        shape.tree.query(bgi::nearest(pt, 1), back_inserter(nearest));
        if (nearest.size()) {
            auto &pair = nearest[0];
            auto dist = bg::distance(pair.first, pt);
            if (dist < tol) {
                norms.push_back(pair);
            }
        }
    }

    double3 ret = { 0, 0, 0 };
    for (auto &pair : norms) {
        ret = ret + pair.second;
    }
    if (norms.size()) {
        ret = ret / norms.size();
    }
    return ret;
}

auto extract_patch(Shape &shape, Point &min, Point &max, int dir, double ext, double len) {
    auto area = bg::area(shape.polys);
    Shape ret;
    if (area < (max.x() - min.x()) * (max.y() - min.y()) * ext * 2) {
        return ret;
    }
    auto tol = fmin(max.x() - min.x(), max.y() - min.y()) * 1e-2;
    //ret.polys = Box(min, max) - shape.polys;
    for (auto &poly : shape.polys) {
        auto outer = poly.outer();
        for (int i = 0, n = outer.size(); i < n - 1; i ++) {
            auto &a = outer[i], &b = outer[i + 1], c = (a + b) / 2;
            if (x_or_y_inside(c, min, max)) {
                auto norm = get_normal(shape, c, tol);
                auto v = pt_2d({ norm }, dir);
                auto l = sqrt(v.x() * v.x() + v.y() * v.y());
                if (l) {
                    auto n = Point { v.x() / l, v.y() / l },
                        e = Point { n.y(), -n.x() };
                    stl::Polygon arrow;
                    bg::append(arrow, c);
                    bg::append(arrow, c - (n - e * 0.1) * len);
                    bg::append(arrow, c - (n + e * 0.1) * len);
                    bg::append(arrow, c);
                    bg::correct(arrow);
                    stl::Polygon edge;
                    bg::append(edge, a);
                    bg::append(edge, a + n * len * 0.1);
                    bg::append(edge, b + n * len * 0.1);
                    bg::append(edge, b);
                    bg::append(edge, a);
                    bg::correct(edge);
                    ret.polys = ret.polys + arrow + edge;
                }
            }
        }
    }
    return ret;
}

map<int, Shape> stl::extract_boundary(map<int, Shape> &frags, grid::Grid &grid, int dir, double ext, double len) {
    map<int, Shape> ret;
    for (auto &pair : frags) {
        auto idx = grid.GetIndex(pair.first);
        auto g = int3 { idx.x, idx.y, idx.z };
        auto min = grid.At(g) - ext / 2, max = grid.At(g + 1) + ext / 2;
        auto shape = extract_patch(pair.second, pt_2d({ min }, dir), pt_2d({ max }, dir), dir, ext, len);
        if (shape.polys.size()) {
            ret[pair.first] = shape;
        }
    }
    return ret;
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
    shape.order = mesh.order;
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

    sort(loops.begin(), loops.end(), [](Loop &a, Loop &b) {
        return abs(a.area) > abs(b.area);
    });
    for (auto &loop : loops) {
        if (loop.hole) {
            shape.polys = shape.polys - make_polys(loop, dir);
        } else {
            shape.polys = shape.polys + make_polys(loop, dir);
        }
    }

    for (auto &loop : loops) {
        for (int i = 0, n = loop.pts.size(); i < n; i ++) {
            auto &a = loop.pts[i], &b = loop.pts[(i + 1) % n];
            auto s = Segment(pt_2d(a, dir), pt_2d(b, dir));
            shape.tree.insert({ s, mesh.normals[a.f.x] });
        }
    }

    return shape;
}

void stl::Spliter::SliceX(grid::Grid &grid, Mesh &mesh, int i) {
    auto a = Slice(mesh, xs[i] - tol/2, DIR_X),
        b = Slice(mesh, xs[i] + tol/2, DIR_X),
        shape = a + b;
    for (int j = 0; j < ny - 1; j ++) {
        auto dy = ys[j + 1] - ys[j];
        auto stride = shape * Box({ ys[j] - ext, min.z }, { ys[j + 1] + ext, max.z });
        for (int k = 0; k < nz - 1; k ++) {
            auto patch = stride * Box({ min.y, zs[k] - ext }, { max.y, zs[k + 1] + ext });
            if (patch.polys.size()) {
                lock_guard<mutex> guard(locks.x);
                fragments.x[grid.GetIndex(i, j, k)] = patch;
            }
        }
    }
}

void stl::Spliter::SliceY(grid::Grid &grid, Mesh &mesh, int j) {
    auto a = Slice(mesh, ys[j] - tol/2, DIR_Y),
        b = Slice(mesh, ys[j] + tol/2, DIR_Y),
        shape = a + b;
    for (int k = 0; k < nz - 1; k ++) {
        auto stride = shape * Box({ zs[k] - ext, min.x }, { zs[k + 1] + ext, max.x });
        for (int i = 0; i < nx - 1; i ++) {
            auto patch = stride * Box({ min.z, xs[i] - ext }, { max.z, xs[i + 1] + ext });
            if (patch.polys.size()) {
                lock_guard<mutex> guard(locks.y);
                fragments.y[grid.GetIndex(i, j, k)] = patch;
            }
        }
    }
}

void stl::Spliter::SliceZ(grid::Grid &grid, Mesh &mesh, int k) {
    auto a = Slice(mesh, zs[k] - tol/2, DIR_Z),
        b = Slice(mesh, zs[k] + tol/2, DIR_Z),
        shape = a + b;
    for (int i = 0; i < nx - 1; i ++) {
        auto stride = shape * Box({ xs[i] - ext, min.y }, { xs[i + 1] + ext, max.y });
        for (int j = 0; j < ny - 1; j ++) {
            auto patch = stride * Box({ min.x, ys[j] - ext }, { max.x, ys[j + 1] + ext });
            if (patch.polys.size()) {
                lock_guard<mutex> guard(locks.z);
                fragments.z[grid.GetIndex(i, j, k)] = patch;
            }
        }
    }
}

stl::Spliter::Spliter(grid::Grid &grid, Mesh &mesh, double tol, double ext) {
    this->tol = tol;
    this->ext = ext;
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
