#include <map>
#include <vector>

#include <windows.h>

#include "vendor/cst.h"
#include "cuda.h"
#include "utils.h"
#include "grid.h"

#ifndef CST_H
#define CST_H

namespace cst {
    typedef struct port_type {
        float3 src, dst;
    } port_type;

    typedef struct excitation_type {
        std::vector<float> x, y;
    } excitation_type;

    typedef struct solid_type {
        std::string name, material;
    } solid_type;

    typedef struct units_type {
        float geometry, time, frequency;
    } units_type;

    typedef struct Project {
        Project(std::string path, std::string version, bool keepCache = false);
        ~Project();

        grid::Grid GetHexGrid();
        float *GetMatrix(int mat);
        double *Get1DResult(std::string tree, int num, int type);

        // will be filled by MakeCacheAndLoadSettings
        float dt = -1;
        std::vector<port_type> ports;
        std::vector<solid_type> solids;
        excitation_type excitation;
        units_type units;

private:
        std::string path;
        std::string version;
        utils::DLL *dll;
        CSTProjHandle handle;

        bool keepCache = false;
        std::string cachePath;
        std::string MakeCacheAndLoadSettings(std::string cst);
    } Project;
}

#endif
