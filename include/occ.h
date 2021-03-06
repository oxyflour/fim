#include <map>
#include <mutex>
#include <TopoDS_Shape.hxx>

#include "grid.h"

#ifndef OCC_H
#define OCC_H

namespace occ {
    using namespace std;

    struct bound_type {
        double xmin, ymin, zmin, xmax, ymax, zmax;
        bool intersects(bound_type &bound, double tol);
    };

    struct Step {
        static TopoDS_Shape load(string file);
        static void save(string file, TopoDS_Shape &shape);
    };

    struct Bool {
        static TopoDS_Shape common(const TopoDS_Shape &a, const TopoDS_Shape &b);
        static TopoDS_Shape fuse(const TopoDS_Shape &a, const TopoDS_Shape &b);
        static TopoDS_Shape cut(const TopoDS_Shape &a, const TopoDS_Shape &b);
    };

    struct Builder {
        static TopoDS_Shape sphere(float3 &position, float radius);
        static TopoDS_Shape box(float3 &min, float3 &max);
        static TopoDS_Shape line(float3 &from, float3 &to);
        static TopoDS_Shape plane(float3 &pos, float3 &dir);
        static TopoDS_Shape component(vector<TopoDS_Shape> &shapes);
    };

    struct Mesh {
        vector<double3> verts;
        vector<int3> faces;
        static Mesh triangulate(TopoDS_Shape &shape);
    };

    struct Shape {
        static bound_type bound(const TopoDS_Shape &shape);
        static vector<TopoDS_Shape> find(TopoDS_Shape &shape, TopAbs_ShapeEnum type);
    };

    struct Mesher {
        Mesher(grid::Grid &grid, TopoDS_Shape &shape, float unit = 1.f);
        void Save(string file);
        static void Merge(grid::Grid &grid, vector<Mesher> &mats, vector<int> &priority);

        vector<float> sx, sy, sz, lx, ly, lz;
        map<int, TopoDS_Shape> msx, msy, msz, mlx, mly, mlz;
    private:
        int nx, ny, nz;
        vector<double> xs, ys, zs;

        TopoDS_Shape shape;
        TopoDS_Shape faces;
        double xmin, ymin, zmin, xmax, ymax, zmax;

        mutex lock;
        void MeshX(int i0, int i1);
        void MeshY(int j0, int j1);
        void MeshZ(int k0, int k1);
    };
};

#endif
