#include <fstream>
#include <iostream>
#include <map>
#include <set>

#include "stl.h"
#include "utils.h"

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

__host__ __device__ __forceinline__ Splited split_face(
    double3 *vertices, int3 *faces, int k, double val, int dir) {
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
    return u.x > 0 && u.x < 1 && u.y > 0 && u.y < 1 ?
                Splited { k, f.x, f.y, lerp(a, b, u.x), f.y, f.z, lerp(b, c, u.y) } :
           u.y > 0 && u.y < 1 && u.z > 0 && u.z < 1 ?
                Splited { k, f.y, f.z, lerp(b, c, u.y), f.z, f.x, lerp(c, a, u.z) } :
                Splited { k, f.z, f.x, lerp(c, a, u.z), f.x, f.y, lerp(a, b, u.x) };
}

__global__ void kernel_split(
    double3 *vertices, int3 *faces, int* bucket, int binLen, double val, int dir,
    Splited *out, int *idx) {
    for (int i = cuIdx(x); i < binLen; i += cuDim(x)) {
        out[atomicAdd(idx, 1)] = split_face(vertices, faces, bucket[i], val, dir);
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

auto get_rings(Mesh &mesh, vector<Splited> &splited) {
    auto maxVert = mesh.vertices.size() + 1;
    auto conns = map<int, Conn>();
    auto coord = map<int, double3>();
    for (int i = 0, n = splited.size(); i < n; i ++) {
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

auto get_loops(Mesh &mesh, vector<Splited> &splited, int dir) {
    auto rings = get_rings(mesh, splited);
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

auto round_by_y(Shape &a, double tol) {
    for (auto &poly : a.polys) {
        for (auto &pt : poly.outer()) {
            pt.set<1>(round_by(pt.get<1>(), tol));
        }
        for (auto &inner : poly.inners()) {
            for (auto &pt : inner) {
                pt.set<1>(round_by(pt.get<1>(), tol));
            }
        }
    }
}

Shape clip_with(Shape &a, Box &b, int dir) {
    Shape shape;
    if (a.polys.size() && bg::intersects(a.bound, b)) {
        auto min = dir == 0 ? b.min_corner().x() : b.min_corner().y(),
            max = dir == 0 ? b.max_corner().x() : b.max_corner().y();
        shape.polys = clip(a.polys, dir, min, max);
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

void stl::export_svg(string file, map<int, Shape> &shapes, function<bool(int)> test) {
    auto list = vector<stl::Shape *>();
    for (auto &pair : shapes) {
        auto &n = pair.first; auto &s = pair.second;
        if (test(n)) {
            list.push_back(&s);
        }
    }
    if (list.size()) {
        ofstream fn(file);
        stl::bg::svg_mapper<stl::Point> map(fn, 500, 500);
        for (auto ptr : list) {
            map.add(ptr->polys);
        }
        for (auto ptr : list) {
            map.map(ptr->polys, "fill:blue;stroke:black;stroke-width:0.1");
        }
    }
}

void Fragments::Dump(string prefix, grid::Grid &grid) {
    for (int i = 0; i < grid.nx; i ++) {
        export_svg(prefix + ".x-" + to_string(i) + "-" + to_string(grid.xs[i]) + ".svg",
            x, [&, i](int n) { return grid.GetIndex(n).x == i; });
    }
    for (int j = 0; j < grid.ny; j ++) {
        export_svg(prefix + ".y-" + to_string(j) + "-" + to_string(grid.ys[j]) + ".svg",
            y, [&, j](int n) { return grid.GetIndex(n).y == j; });
    }
    for (int k = 0; k < grid.nz; k ++) {
        export_svg(prefix + ".z-" + to_string(k) + "-" + to_string(grid.zs[k]) + ".svg",
            z, [&, k](int n) { return grid.GetIndex(n).z == k; });
    }
}

Fragments Fragments::GetBoundary(grid::Grid &grid, double tol, double len) {
    return Fragments {
        extract_boundary(x, grid, DIR_X, tol, len),
        extract_boundary(y, grid, DIR_Y, tol, len),
        extract_boundary(z, grid, DIR_Z, tol, len),
    };
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

inline auto x_or_y_inside(const Point &pt, const Point &min, const Point &max, double ext = 0) {
    auto x = pt.x(), y = pt.y();
    return (x > min.x() + ext && x < max.x() - ext) || (y > min.y() + ext && y < max.y() - ext);
}

inline auto get_normal(const Shape &shape, Point &pt, double tol) {
    vector<RTPair> norms;
    shape.tree.query(bgi::intersects(pt), back_inserter(norms));

    int count = 0;
    auto ret = double3 { 0, 0, 0 };
    for (auto &pair : norms) {
        auto dist = bg::distance(pair.second.segment, pt);
        if (dist < tol) {
            count ++;
            ret = ret + pair.second.normal;
        }
    }
    return count > 1 ? ret / count : ret;
}

auto extract_patch(const Shape &shape, const Point &min, const Point &max, int dir, double ext, double len) {
    auto area = bg::area(shape.polys);
    Shape ret;
    if (area < (max.x() - min.x()) * (max.y() - min.y()) * ext * 2) {
        return ret;
    }
    auto tol = fmin(max.x() - min.x(), max.y() - min.y()) * 1e-2;
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

vector<Splited> stl::Spliter::Split(Mesh &mesh, double pos, int dir) {
    auto vertices = mesh.vertices.data();
    auto faces = mesh.faces.data();
    auto splitedArr = vector<Splited>();
    if (faceNum < 100) {
        for (int k = 0, n = mesh.faces.size(); k < n; k ++) {
            auto &b = mesh.bounds[k];
            auto min = dir == DIR_X ? b.min.x : dir == DIR_Y ? b.min.y : b.min.z,
                max = dir == DIR_X ? b.max.x : dir == DIR_Y ? b.max.y : b.max.z;
            if (min < pos && pos < max) {
                splitedArr.push_back(split_face(vertices, faces, k, pos, dir));
            }
        }
    } else {
        auto bucket = malloc_device<int>(faceNum);
        auto idx = malloc_device<int>(1);
        auto splited = malloc_device<Splited>(faceNum);

        int bucketLen = 0;
        to_device(&bucketLen, 1, idx);
        kernel_group CU_DIM(256, 128) (vertices, faces, faceNum, bounds, pos, dir, bucket, idx);
        CU_ASSERT(cudaGetLastError());
        from_device(idx, 1, &bucketLen);

        if (bucketLen > 0) {
            int splitedLen = 0;
            to_device(&splitedLen, 1, idx);
            kernel_split CU_DIM(256, 128) (vertices, faces, bucket, bucketLen, pos, dir, splited, idx);
            CU_ASSERT(cudaGetLastError());
            from_device(idx, 1, &splitedLen);
            splitedArr.resize(splitedLen);
            from_device(splited, splitedLen, splitedArr.data());
        }

        cudaFree(bucket);
        cudaFree(idx);
        cudaFree(splited);
    }
    return splitedArr;
}

Shape stl::Spliter::Slice(Mesh &mesh, double pos, int dir) {
    auto splited = Split(mesh, pos, dir);
    auto loops = get_loops(mesh, splited, dir);

    Shape shape;
    shape.order = mesh.order;
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
            auto u = pt_2d(a, dir), v = pt_2d(b, dir);
            auto r = Box(
                { fmin(u.x(), v.x()) - tol, fmin(u.y(), v.y()) - tol },
                { fmax(u.x(), v.x()) + tol, fmax(u.y(), v.y()) + tol }
            );
            shape.tree.insert({ r, { Segment(u, v), mesh.normals[a.f.x] } });
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
        auto stride = clip_with(shape, Box({ ys[j] - tol/2, min.z }, { ys[j + 1] + tol/2, max.z }), 0); round_by_y(stride, tol);
        //auto stride = shape * Box({ ys[j] - ext, min.z }, { ys[j + 1] + ext, max.z });
        for (int k = 0; k < nz - 1; k ++) {
            auto patch = clip_with(stride, Box({ min.y, zs[k] - tol/2 }, { max.y, zs[k + 1] + tol/2 }), 1);
            //auto patch = stride * Box({ min.y, zs[k] - ext }, { max.y, zs[k + 1] + ext });
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
        auto stride = clip_with(shape, Box({ zs[k] - tol/2, min.x }, { zs[k + 1] + tol/2, max.x }), 0); round_by_y(stride, tol);
        //auto stride = shape * Box({ zs[k] - ext, min.x }, { zs[k + 1] + ext, max.x });
        for (int i = 0; i < nx - 1; i ++) {
            auto patch = clip_with(stride, Box({ min.z, xs[i] - tol/2 }, { max.z, xs[i + 1] + tol/2 }), 1);
            //auto patch = stride * Box({ min.z, xs[i] - ext }, { max.z, xs[i + 1] + ext });
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
        auto stride = clip_with(shape, Box({ xs[i] - tol/2, min.y }, { xs[i + 1] + tol/2, max.y }), 0); round_by_y(stride, tol);
        //auto stride = shape * Box({ xs[i] - ext, min.y }, { xs[i + 1] + ext, max.y });
        for (int j = 0; j < ny - 1; j ++) {
            auto patch = clip_with(stride, Box({ min.x, ys[j] - tol/2 }, { max.x, ys[j + 1] + tol/2 }), 1);
            //auto patch = stride * Box({ min.x, ys[j] - ext }, { max.x, ys[j + 1] + ext });
            if (patch.polys.size()) {
                lock_guard<mutex> guard(locks.z);
                fragments.z[grid.GetIndex(i, j, k)] = patch;
            }
        }
    }
}

stl::Spliter::Spliter(grid::Grid &grid, Mesh &mesh, ctpl::thread_pool &pool, double tol, double ext) {
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

    auto futures = vector<future<void>>();
    for (int i = 0; i < nx; i ++) {
        if (mesh.min.x <= grid.xs[i] && grid.xs[i] <= mesh.max.x) {
            futures.push_back(pool.push([&, i](int id) { SliceX(grid, mesh, i); }));
        }
    }
    for (int j = 0; j < ny; j ++) {
        if (mesh.min.y <= grid.ys[j] && grid.ys[j] <= mesh.max.y) {
            futures.push_back(pool.push([&, j](int id) { SliceY(grid, mesh, j); }));
        }
    }
    for (int k = 0; k < nz; k ++) {
        if (mesh.min.z <= grid.zs[k] && grid.zs[k] <= mesh.max.z) {
            futures.push_back(pool.push([&, k](int id) { SliceZ(grid, mesh, k); }));
        }
    }
    for (auto &future : futures) {
        future.wait();
    }
}

stl::Spliter::~Spliter() {
    cudaFree(faces);
    cudaFree(vertices);
    cudaFree(normals);
    cudaFree(bounds);
}

struct ClipJoint {
    int i, j, k;
    Point p;
    bool removed;
    ClipJoint *next;
};
struct Ending {
    ClipJoint *a, *b;
};
struct PointArray {
    Point *ptr;
    size_t num;
};
auto slice(stl::Polygon &poly, int dir, double pos, int side) {
    vector<PointArray> rings(1 + poly.inners().size());
    auto &outer = poly.outer();
    rings[0] = { outer.data(), outer.size() - 1 };
    auto &inners = poly.inners();
    for (int k = 0; k < inners.size(); k ++) {
        rings[k + 1] = { inners[k].data(), inners[k].size() - 1 };
    }

    vector<ClipJoint *> joints;
    vector<vector<Ending>> ends;
    vector<int> outside;
    for (auto &ring : rings) {
        ends.push_back(vector<Ending>(ring.num));
    }
    for (int k = 0; k < rings.size(); k ++) {
        auto &r = rings[k];
        auto &e = ends[k];
        auto is_sliced = false;
        for (int i = 0, n = r.num; i < n; i ++) {
            auto j = (i + 1) % n;
            auto &a = r.ptr[i], &b = r.ptr[j];
            auto test = dir == 0 ? (a.x() - pos) * (b.x() - pos)   : (a.y() - pos) * (b.y() - pos),
                 fac =  dir == 0 ? (a.x() - pos) / (a.x() - b.x()) : (a.y() - pos) / (a.y() - b.y());
            if (test < 0) {
                auto p = a * (1 - fac) + b * fac;
                auto joint = new ClipJoint { i, j, k, p, false, NULL };
                joints.push_back(joint);
                e[i].a == NULL ? (e[i].a = joint) : (e[i].b = joint);
                e[j].a == NULL ? (e[j].a = joint) : (e[j].b = joint);
                is_sliced = true;
            }
        }
        if (!is_sliced && r.num > 0) {
            auto p = r.ptr[0];
            auto position = dir == 0 ? p.x() : p.y();
            if ((position - pos) * side > 0) {
                outside.push_back(k);
            }
        }
    }
    if (joints.size() % 2 != 0) {
        throw runtime_error(string("slice at ") +
            (dir ? "x" : "y") + "=" + to_string(pos) + " failed " +
            "(" + to_string(joints.size()) + " joints)");
    }
    sort(joints.begin(), joints.end(), [dir](ClipJoint *a, ClipJoint *b) {
        return dir == 0 ? a->p.y() < b->p.y() : a->p.x() < b->p.x();
    });

    MultiPolygon ret;
    if (!joints.size()) {
        if (outer.size()) {
            auto p = outer[0];
            auto position = dir == 0 ? p.x() : p.y();
            if ((position - pos) * side > 0) {
                ret.push_back(poly);
            }
        }
        return ret;
    }
    for (int i = 0; i < joints.size() - 1; i ++) {
        joints[i]->next = joints[i + 1];
    }

    auto iter = find_if(joints.begin(), joints.end(), [](ClipJoint *a) { return !a->removed; });
    while (iter != joints.end()) {
        auto idx = distance(joints.begin(), iter);
        auto a = joints[idx], b = a->next;

        stl::Polygon poly;
        auto &ring = poly.outer();
        while (a && b) {
            ring.push_back(a->p);
            ring.push_back(b->p);
            a->removed = b->removed = true;

            auto &r = rings[b->k];
            auto &e = ends[b->k];
            auto &p = r.ptr[b->i];
            auto inversed = ((dir == 0 ? p.x() : p.y()) - pos) * side > 0;
            auto i = inversed ? b->i : b->j;
            auto d = inversed ? -1 : 1;
            while (true) {
                ring.push_back(r.ptr[i]);
                auto &ending = e[i];
                auto joint = ending.a == b ? ending.b : ending.a;
                if (joint && joint != b) {
                    if (joint->removed) {
                        a = b = NULL;
                    } else {
                        a = joint;
                        b = joint->next;
                    }
                    break;
                } else {
                    i = (i + d + r.num) % r.num;
                }
            }
        }
        bg::correct(poly);
        ret.push_back(poly);

        iter = find_if(joints.begin(), joints.end(), [](ClipJoint *a) { return !a->removed; });
    }
    for (auto joint : joints) {
        delete joint;
    }

    for (auto k : outside) {
        auto p = stl::Polygon();
        auto &out = p.outer();
        out.insert(out.begin(), inners[k - 1].begin(), inners[k - 1].end());
        bg::correct(p);
        ret = ret - p;
    }
    return ret;
}

MultiPolygon stl::clip(MultiPolygon &c, int dir, double min, double max) {
    MultiPolygon e, f;
    for (auto &p : c) {
        e = e + slice(p, dir, min, 1);
    }
    for (auto &p : e) {
        f = f + slice(p, dir, max, -1);
    }
    return f;
}
