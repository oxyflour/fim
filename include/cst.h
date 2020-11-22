#include <map>
#include <vector>

#include <windows.h>

#include "utils/cst_interface.h"
#include "utils/assert.h"

#include "grid.h"
#include "utils.h"

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
        Project(std::string &path, std::string &version);
        ~Project();

        Grid GetHexGrid();
        std::vector<float> GetMatrix(int mat);
        std::vector<double> Get1DResult(std::string tree, int num, int type);

        // will be filled by ForkAndExportSettings
        std::vector<port_type> ports;
        std::vector<solid_type> solids;
        excitation_type excitation;
        units_type units;

private:
        std::string path;
        std::string version;
        utils::DLL *dll;
        CSTProjHandle handle;

        void ForkAndExportSettings(std::string cst);
    } Project;
}

#endif
