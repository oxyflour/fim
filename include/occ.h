#include "grid.h"

#include <TopoDS_Shape.hxx>

#ifndef OCC_H
#define OCC_H

namespace occ {
    typedef struct bound_type {
        double xmin, ymin, zmin, xmax, ymax, zmax;
    } bound_type;

    typedef struct Step {
        static TopoDS_Shape load(std::string &file);
        static void save(std::string &file, TopoDS_Shape &shape);
    } Step;

    typedef struct Mesher {
        Mesher(Grid &grid, float unit, std::string &file);
    private:
        int nx, ny, nz;
        std::vector<double> xs, ys, zs;
        float unit;
        TopoDS_Shape shape;
        void MeshX();
    } Mesher;
};

#endif
