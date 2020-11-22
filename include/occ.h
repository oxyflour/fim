#include "grid.h"

#include <TopoDS.hxx>

#ifndef OCC_H
#define OCC_H

namespace occ {
    typedef struct Mesher {
        Mesher(Grid &grid, std::string &file);
    private:
        Grid *grid;
        TopoDS_Shape shape;
        std::vector<TopoDS_Shape> MeshX();
    } Mesher;
};

#endif
