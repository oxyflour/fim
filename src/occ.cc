#include "occ.h"
#include "utils/assert.h"

#include <vector>
#include <algorithm>

#include <gp_Pln.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>

#include <TopoDS_Face.hxx>
#include <TopExp_Explorer.hxx>

#include <BRepAlgoAPI_Common.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepBndLib.hxx>

#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>

using namespace std;

TopoDS_Shape occ::Step::load(string &file) {
    STEPControl_Reader reader;
    auto stat = reader.ReadFile(file.c_str());
    reader.TransferRoot();
    ASSERT(stat == IFSelect_RetDone, "load from file " + file + " failed");
    return reader.Shape();
}

void occ::Step::save(string &file, TopoDS_Shape &shape) {
    STEPControl_Writer writer;
    auto stat = writer.Transfer(shape, STEPControl_StepModelType::STEPControl_AsIs);
    ASSERT(stat == IFSelect_RetDone, "save to file " + file + " failed");
    writer.Write(file.c_str());
}

static auto bound(TopoDS_Shape &shape) {
    Bnd_Box box;
    BRepBndLib::Add(shape, box);
    occ::bound_type ret;
    box.Get(ret.xmin, ret.ymin, ret.zmin, ret.xmax, ret.ymax, ret.zmax);
    return ret;
}

static auto find(TopoDS_Shape &shape, TopAbs_ShapeEnum type) {
    TopExp_Explorer exp;
    vector<TopoDS_Shape> ret;
    for (exp.Init(shape, type); exp.More(); exp.Next()) {
        ret.push_back(exp.Current());
    }
    return ret;
}

static auto component(vector<TopoDS_Shape> &shapes) {
    TopoDS_Compound ret;
    BRep_Builder builder;
    builder.MakeCompound(ret);
    for (auto shape : shapes) {
        builder.Add(ret, shape);
    }
    return ret;
}

static auto common(const TopoDS_Shape &a, const TopoDS_Shape &b) {
    BRepAlgoAPI_Common api;
    TopTools_ListOfShape args, tools;
    args.Append(a);
    api.SetArguments(args);
    tools.Append(b);
    api.SetTools(tools);
    api.Build();
    return api.Shape();
}

static auto getLinearProps(const TopoDS_Shape &shape) {
    GProp_GProps props;
    BRepGProp::LinearProperties(shape, props);
    return props;
}

static auto getSurfaceProps(const TopoDS_Shape &shape) {
    GProp_GProps props;
    BRepGProp::SurfaceProperties(shape, props);
    return props;
}

occ::Mesher::Mesher(Grid &grid, float unit, string &file) {
    xs = grid.xs; ys = grid.ys; zs = grid.zs;
    nx = xs.size(); ny = ys.size(); nz = zs.size();
    for (int i = 0; i < nx; i ++) xs[i] /= unit;
    for (int i = 0; i < ny; i ++) ys[i] /= unit;
    for (int i = 0; i < nz; i ++) zs[i] /= unit;
    shape = Step::load(file);
    MeshX();
}

typedef struct mesh_part {
    TopoDS_Shape shape;
    int i, j, k;
} mesh_part;

void occ::Mesher::MeshX() {
    auto areaArr = new float[nx * ny * nz];
    auto lengthArr = new float[nx * ny * nz];

    auto b = bound(shape);
    vector<TopoDS_Shape> shapes;

    auto boxY = vector<TopoDS_Shape>(ny);
    for (int j = 0; j < ny - 1; j ++) {
        boxY[j] = BRepPrimAPI_MakeBox(gp_Pnt(b.xmin, ys[j], b.zmin), gp_Pnt(b.xmax, ys[j + 1], b.zmax)).Shape();
    }

    for (int i = 0; i < nx - 1; i ++) {
        auto x0 = xs[i], x1 = xs[i + 1];
        if (!(b.xmin < x1 && x0 < b.xmax) || i != nx / 2) {
            continue;
        }

        TopoDS_Face plane = BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(x0, 0, 0), gp_Dir(1, 0, 0))).Face();
        auto c1 = common(plane, shape);
        auto b1 = bound(c1);
        auto w1 = find(c1, TopAbs_ShapeEnum::TopAbs_WIRE);
        if (w1.size() == 0) {
            continue;
        }

        for (int j = 0; j < ny - 1; j ++) {
            auto y0 = ys[j], y1 = ys[j + 1];
            if (!(b1.ymin < y1 && y0 < b1.ymax)) {
                continue;
            }

            auto c2 = common(c1, boxY[j]);
            auto b2 = bound(c2);
            auto w2 = vector<TopoDS_Shape>();
            for (auto w : w1) {
                w2.push_back(common(w, boxY[j]));
            }
            for (int k = 0; k < nz - 1; k ++) {
                auto z0 = zs[k], z1 = zs[k + 1];
                if (!(b2.zmin < z1 && z0 < b2.zmax)) {
                    continue;
                }

                auto idx = i + j * nx + k * nx * ny;
                auto box = BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape();
                auto area = common(c2, box);
                areaArr[idx] = getSurfaceProps(area).Mass();

                auto line = BRepBuilderAPI_MakeEdge(gp_Pnt(x0, y0, z0), gp_Pnt(x0, y1, z0)).Shape();
                auto length = common(c2, line);
                lengthArr[idx] = getLinearProps(length).Mass();

                // debug
                //shapes.push_back(area);
                //shapes.push_back(length);
                for (auto wire : w2) {
                    shapes.push_back(common(wire, box));
                }
                for (auto wire : w2) {
                    shapes.push_back(common(wire, line));
                }
            }
        }
    }
    // FIXME:
    Step::save(string("build/test.stp"), component(shapes));
}
