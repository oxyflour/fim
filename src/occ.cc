#include "occ.h"
#include "utils/assert.h"

#include <vector>
#include <algorithm>

#include <gp_Pln.hxx>
#include <Bnd_Box.hxx>

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

static TopoDS_Shape load(string &file) {
    STEPControl_Reader reader;
    auto stat = reader.ReadFile(file.c_str());
    reader.TransferRoot();
    ASSERT(stat == IFSelect_RetDone, "load from file " + file + " failed");
    return reader.Shape();
}

static void save(string &file, TopoDS_Shape &shape) {
    STEPControl_Writer writer;
    auto stat = writer.Transfer(shape, STEPControl_StepModelType::STEPControl_AsIs);
    ASSERT(stat == IFSelect_RetDone, "save to file " + file + " failed");
    writer.Write(file.c_str());
}

typedef struct bound_type {
    double xmin, ymin, zmin, xmax, ymax, zmax;
} bound_type;

static Bnd_Box bound(TopoDS_Shape &shape) {
    Bnd_Box box;
    BRepBndLib::Add(shape, box);
    bound_type ret;
    box.Get(ret.xmin, ret.ymin, ret.zmin, ret.xmax, ret.ymax, ret.zmax);
    return ret;
}

static TopoDS_Compound component(vector<TopoDS_Shape> &shapes) {
    TopoDS_Compound ret;
    BRep_Builder builder;
    builder.MakeCompound(ret);
    for (auto shape : shapes) {
        builder.Add(ret, shape);
    }
    return ret;
}

static TopoDS_Shape common(const TopoDS_Shape &a, const TopoDS_Shape &b) {
    BRepAlgoAPI_Common api;
    TopTools_ListOfShape args, tools;
    args.Append(a);
    api.SetArguments(args);
    tools.Append(b);
    api.SetTools(tools);
    api.Build();
    return api.Shape();
}

static vector<TopoDS_Shape> find(const TopoDS_Shape shape, TopAbs_ShapeEnum type) {
    TopExp_Explorer exp;
    vector<TopoDS_Shape> ret;
    for (exp.Init(shape, type); exp.More(); exp.Next()) {
        ret.push_back(exp.Current());
    }
    return ret;
}

occ::Mesher::Mesher(Grid &grid, string &file) {
    this->grid = &grid;
    shape = load(file);
    save(file + ".split.stp", component(MeshX()));
}

vector<TopoDS_Shape> occ::Mesher::MeshX() {
    auto xs = grid->xs, ys = grid->ys, zs = grid->zs;
    for (int i = 0, n = xs.size(); i < n; i ++) xs[i] /= grid->unit;
    for (int i = 0, n = ys.size(); i < n; i ++) ys[i] /= grid->unit;
    for (int i = 0, n = zs.size(); i < n; i ++) zs[i] /= grid->unit;

    auto b1 = bound(shape);
    vector<TopoDS_Shape> shapes;
    for (int i = 0, nx = xs.size(); i < nx - 1; i ++) {
        auto x0 = xs[i], x1 = xs[i + 1];
        if (!(b1.xmin < x1 && x0 < b1.xmax) || i != nx / 2) {
            continue;
        }

        TopoDS_Face plane = BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(x0, 0, 0), gp_Dir(1, 0, 0))).Face();
        auto c1 = common(plane, shape);
        auto b2 = bound(c1);
        auto wires = find(c1, TopAbs_ShapeEnum::TopAbs_WIRE);
        if (wires.size() == 0) {
            continue;
        }

        for (int j = 0, ny = ys.size(); j < ny - 1; j ++) {
            auto y0 = ys[j], y1 = ys[j + 1];
            if (!(b2.ymin < y1 && y0 < b2.ymax)) {
                continue;
            }

            auto box = BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, b2.zmin), gp_Pnt(x1, y1, b2.zmax)).Shape();
            auto c2 = common(c1, box);
            auto b3 = bound(c2);

            for (int k = 0, nz = zs.size(); k < nz - 1; k ++) {
                auto z0 = zs[k], z1 = zs[k + 1];
                if (!(b3.zmin < z1 && z0 < b3.zmax)) {
                    continue;
                }

                auto box = BRepPrimAPI_MakeBox(gp_Pnt(x0, y0, z0), gp_Pnt(x1, y1, z1)).Shape();
                shapes.push_back(common(c2, box));
                for (auto wire : wires) {
                    shapes.push_back(common(wire, box));
                }

                auto line = BRepBuilderAPI_MakeEdge(gp_Pnt(x0, y0, z0), gp_Pnt(x0, y1, z0)).Shape();
                shapes.push_back(common(c2, line));
                for (auto wire : wires) {
                    shapes.push_back(common(wire, line));
                }
            }
        }
    }
    return shapes;
}
