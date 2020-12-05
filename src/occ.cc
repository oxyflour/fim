#include "occ.h"
#include "utils/assert.h"
#include "utils.h"

#include <vector>
#include <algorithm>

#include <gp_Pln.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>

#include <BRepAlgoAPI_Common.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepBndLib.hxx>

#include <BOPTools_AlgoTools3D.hxx>

#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>
#include <Message_Messenger.hxx>
#include <Message.hxx>
#include <Message_PrinterOStream.hxx>
#include <XSControl_WorkSession.hxx>

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>

using namespace std;
using namespace occ;

// https://stackoverflow.com/questions/64912439/preventing-opencascade-to-write-to-console
//Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));

TopoDS_Shape Step::load(string &file) {
    STEPControl_Reader reader;
    auto stat = reader.ReadFile(file.c_str());
    reader.TransferRoot();
    ASSERT(stat == IFSelect_RetDone, "load from file " + file + " failed");
    return reader.Shape();
}

void Step::save(string &file, TopoDS_Shape &shape) {
    STEPControl_Writer writer;
    auto stat = writer.Transfer(shape, STEPControl_StepModelType::STEPControl_AsIs);
    ASSERT(stat == IFSelect_RetDone, "save to file " + file + " failed");
    writer.Write(file.c_str());
}

bound_type Shape::bound(const TopoDS_Shape &shape) {
    Bnd_Box box;
    BRepBndLib::Add(shape, box);
    bound_type ret;
    box.Get(ret.xmin, ret.ymin, ret.zmin, ret.xmax, ret.ymax, ret.zmax);
    return ret;
}

vector<TopoDS_Shape> Shape::find(TopoDS_Shape &shape, TopAbs_ShapeEnum type) {
    TopExp_Explorer exp;
    vector<TopoDS_Shape> ret;
    for (exp.Init(shape, type); exp.More(); exp.Next()) {
        ret.push_back(exp.Current());
    }
    return ret;
}

TopoDS_Shape Builder::component(vector<TopoDS_Shape> &shapes) {
    TopoDS_Compound ret;
    BRep_Builder builder;
    builder.MakeCompound(ret);
    for (auto shape : shapes) {
        builder.Add(ret, shape);
    }
    return ret;
}

TopoDS_Shape Bool::common(const TopoDS_Shape &a, const TopoDS_Shape &b) {
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

Mesher::Mesher(Grid &grid, string &file, float unit) {
    xs = grid.xs; ys = grid.ys; zs = grid.zs;
    nx = xs.size(); ny = ys.size(); nz = zs.size();
    for (int i = 0; i < nx; i ++) xs[i] /= unit;
    for (int i = 0; i < ny; i ++) ys[i] /= unit;
    for (int i = 0; i < nz; i ++) zs[i] /= unit;
    shape = Step::load(file);
    MeshX();
}

template <typename T> auto intersects(T a0, T a1, T b0, T b1, T tol = 1e-3) {
    auto v = (a1 - a0) * tol;
    return a0 <= b1 + v && b0 <= a1 + v;
}

bool bound_type::intersects(bound_type &b, double tol = 1e-3) {
    return ::intersects(xmin, xmax, b.xmin, b.xmax, tol) &&
        ::intersects(xmin, xmax, b.xmin, b.xmax, tol) &&
        ::intersects(xmin, xmax, b.xmin, b.xmax, tol);
}

TopoDS_Shape Builder::sphere(float3 &position, float radius) {
    return BRepPrimAPI_MakeSphere(gp_Pnt(position.x, position.y, position.z), radius).Shape();
}

TopoDS_Shape Builder::plane(float3 &pos, float3 &dir) {
    return BRepBuilderAPI_MakeFace(gp_Pln(gp_Pnt(pos.x, pos.y, pos.z), gp_Dir(dir.x, dir.y, dir.z))).Face();
}

TopoDS_Shape Builder::box(float3 &min, float3 &max) {
    return BRepPrimAPI_MakeBox(gp_Pnt(min.x, min.y, min.z), gp_Pnt(max.x, max.y, max.z)).Shape();
}

TopoDS_Shape Builder::line(float3 &from, float3 &to) {
    return BRepBuilderAPI_MakeEdge(gp_Pnt(from.x, from.y, from.z), gp_Pnt(to.x, to.y, to.z)).Shape();
}

typedef struct shape_bounds {
    const TopoDS_Shape shape;
    bound_type bound;
} shape_bounds;

void calcWireNormal(shape_bounds &box, vector<shape_bounds> &faces, int dir) {
    for (auto &fv : faces) {
        if (fv.bound.intersects(box.bound)) {
            for (auto &fc : Shape::find(Bool::common(fv.shape, box.shape), TopAbs_FACE)) {
                auto face = TopoDS::Face(fc);
                for (auto &ed : Shape::find(fc, TopAbs_EDGE)) {
                    auto bound = Shape::bound(ed);
                    if ((dir == 0 &&
                            abs(bound.xmin - box.bound.xmin) < (box.bound.xmax - box.bound.xmin) * 1e-3 &&
                            abs(bound.xmax - box.bound.xmin) < (box.bound.xmax - box.bound.xmin) * 1e-3)) {
                        auto edge = TopoDS::Edge(ed);
                        gp_Dir dir;
                        BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(edge, face, 0.5, dir);
                        printf("%f, %f, %f / %f, %f, %f\n",
                            box.bound.xmin, box.bound.ymin, box.bound.zmin, dir.X(), dir.Y(), dir.Z());
                    }
                }
            }
        }
    }
}

void Mesher::MeshX() {
    auto areaArr = vector<float>(nx * ny * nz);
    auto lengthArr = vector<float>(nx * ny * nz);

    auto surfMap = map<int, TopoDS_Shape>();
    auto edgeMap = map<int, TopoDS_Shape>();

    auto faceBounds = utils::map(Shape::find(shape, TopAbs_FACE), [](const TopoDS_Shape &shape){
        return shape_bounds{ shape, Shape::bound(shape) };
    });

    auto b = Shape::bound(shape);
    auto shapes = vector<TopoDS_Shape>();
    for (int i = 0; i < nx - 1; i ++) {
        auto x0 = xs[i], x1 = xs[i + 1], dx = x1 - x0;
        if (!intersects(x0, x1, b.xmin, b.xmax)) {
            continue;
        }

        auto plane = Builder::plane(float3::create(x0, 0, 0), float3::create(1, 0, 0));
        auto c1 = Bool::common(plane, shape);
        if (BOPTools_AlgoTools3D::IsEmptyShape(c1)) {
            continue;
        }

        auto b1 = Shape::bound(c1);
        for (int j = 0; j < ny - 1; j ++) {
            auto y0 = ys[j], y1 = ys[j + 1], dy = y1 - y0;
            if (!intersects(y0, y1, b1.ymin, b1.ymax)) {
                continue;
            }

            auto box = Builder::box(float3::create(x0, y0, b.zmin), float3::create(x1, y1, b.zmax));
            auto c2 = Bool::common(c1, box);
            if (BOPTools_AlgoTools3D::IsEmptyShape(c2)) {
                continue;
            }

            auto b2 = Shape::bound(c2);
            for (int k = 0; k < nz - 1; k ++) {
                auto z0 = zs[k], z1 = zs[k + 1], dz = z1 - z0;
                if (!intersects(z0, z1, b2.zmin, b2.zmax)) {
                    continue;
                }

                auto calcWire = false;

                auto idx = i + j * nx + k * nx * ny;
                auto box = Builder::box(float3::create(x0, y0, z0), float3::create(x1, y1, z1));
                auto surface = Bool::common(c2, box);
                auto area = areaArr[idx] = getSurfaceProps(surface).Mass();
                if (area > 0.001 * dx * dy && area < 0.999 * dx * dy) {
                    surfMap[idx] = surface;
                    calcWire = true;
                }

                auto line = Builder::line(float3::create(x0, y0, z0), float3::create(x0, y1, z0));
                auto edge = Bool::common(c2, line);
                auto length = lengthArr[idx] = getLinearProps(edge).Mass();
                if (length > 0.01 * dy && length < 0.99 * dy) {
                    edgeMap[idx] = edge;
                    calcWire = true;
                }

                if (calcWire) {
                    auto bound = bound_type { x0, y0, z0, x1, y1, z1 };
                    auto boxBound = shape_bounds { box, bound };
                    calcWireNormal(boxBound, faceBounds, 0);
                }
            }
        }
    }
    /*
    for (auto pair : edgeMap) {
        shapes.push_back(pair.second);
    }
    for (auto pair : surfMap) {
        shapes.push_back(pair.second);
    }
     */
    Step::save(string("E:\\out.stp"), Builder::component(shapes));
}
