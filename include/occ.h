#include "grid.h"

#include <map>
#include <TopoDS_Shape.hxx>

#ifndef OCC_H
#define OCC_H

namespace occ {
    typedef struct bound_type {
        double xmin, ymin, zmin, xmax, ymax, zmax;
        bool intersects(bound_type &bound, double tol);
    } bound_type;

    typedef struct Step {
        static TopoDS_Shape load(std::string &file);
        static void save(std::string &file, TopoDS_Shape &shape);
    } Step;

    typedef struct Bool {
        static TopoDS_Shape common(const TopoDS_Shape &a, const TopoDS_Shape &b);
    } Bool;

    typedef struct Builder {
        static TopoDS_Shape sphere(float3 &position, float radius);
        static TopoDS_Shape box(float3 &min, float3 &max);
        static TopoDS_Shape line(float3 &from, float3 &to);
        static TopoDS_Shape plane(float3 &pos, float3 &dir);
        static TopoDS_Shape component(std::vector<TopoDS_Shape> &shapes);
    } Builder;

    typedef struct Shape {
        static bound_type bound(const TopoDS_Shape &shape);
        static std::vector<TopoDS_Shape> find(TopoDS_Shape &shape, TopAbs_ShapeEnum type);
    } Shape;

    typedef struct Mesher {
        Mesher(Grid &grid, std::string &file, float unit = 1.f);
    private:
        int nx, ny, nz;
        std::vector<double> xs, ys, zs;
        std::vector<float> sx, sy, sz, lx, ly, lz;
        std::map<int, TopoDS_Shape> msx, msy, msz, mlx, mly, mlz;
        TopoDS_Shape shape;
        TopoDS_Shape faces;
        double xmin, ymin, zmin, xmax, ymax, zmax;
        void MeshX();
        void MeshY();
        void MeshZ();
    } Mesher;
};

#endif
