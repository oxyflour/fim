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
    auto bound = Shape::bound(shape);
    xmin = bound.xmin; ymin = bound.ymin; zmin = bound.zmin;
    xmax = bound.xmax; ymax = bound.ymax; zmax = bound.zmax;

    faceBounds = utils::map(Shape::find(shape, TopAbs_FACE), [](const TopoDS_Shape &shape){
        return shape_bounds{ shape, Shape::bound(shape) };
    });

    auto nxyz = nx * ny * nz;
    sx.resize(nxyz); std::fill(sx.begin(), sx.end(), 0);
    sy.resize(nxyz); std::fill(sy.begin(), sy.end(), 0);
    sz.resize(nxyz); std::fill(sz.begin(), sz.end(), 0);
    lx.resize(nxyz); std::fill(lx.begin(), lx.end(), 0);
    ly.resize(nxyz); std::fill(ly.begin(), ly.end(), 0);
    lz.resize(nxyz); std::fill(lz.begin(), lz.end(), 0);

    //MeshX();
    MeshY();
    //MeshZ();

    vector<TopoDS_Shape> shapes;
    for (auto pair : msx) shapes.push_back(pair.second);
    for (auto pair : msy) shapes.push_back(pair.second);
    for (auto pair : msz) shapes.push_back(pair.second);
    for (auto pair : mlx) shapes.push_back(pair.second);
    for (auto pair : mly) shapes.push_back(pair.second);
    for (auto pair : mlz) shapes.push_back(pair.second);
    Step::save(string("E:\\out.stp"), Builder::component(shapes));
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

void calcWireNormal(bound_type &bound, vector<shape_bounds> &faces, int dir) {
    auto xmin = bound.xmin, ymin = bound.ymin, zmin = bound.zmin,
        xmax = bound.xmax, ymax = bound.ymax, zmax = bound.zmax;
    auto box = Builder::box(float3(xmin, ymin, zmin), float3(xmax, ymax, zmax));
    for (auto &fv : faces) {
        if (fv.bound.intersects(bound)) {
            for (auto &fc : Shape::find(Bool::common(fv.shape, box), TopAbs_FACE)) {
                auto face = TopoDS::Face(fc);
                for (auto &ed : Shape::find(fc, TopAbs_EDGE)) {
                    auto bound = Shape::bound(ed);
                    if ((dir == 0 &&
                            abs(bound.xmin - xmin) < (xmax - xmin) * 1e-3 &&
                            abs(bound.xmax - xmin) < (xmax - xmin) * 1e-3) ||
                        (dir == 1 &&
                            abs(bound.ymin - ymin) < (ymax - ymin) * 1e-3 &&
                            abs(bound.ymax - ymin) < (ymax - ymin) * 1e-3) ||
                        (dir == 2 &&
                            abs(bound.zmin - zmin) < (zmax - zmin) * 1e-3 &&
                            abs(bound.zmax - zmin) < (zmax - zmin) * 1e-3)) {
                        auto edge = TopoDS::Edge(ed);
                        gp_Dir dir;
                        BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(edge, face, 0.5, dir);
                        // TODO: output
                        //printf("%f, %f, %f / %f, %f, %f\n", xmin, ymin, zmin, dir.X(), dir.Y(), dir.Z());
                    }
                }
            }
        }
    }
}

void Mesher::MeshX() {
    auto &areaArr = sx, &lengthArr = ly;
    auto &surfMap = msx, &edgeMap = mly;

    for (int i = 0; i < nx - 1; i ++) {
        auto x0 = xs[i], x1 = xs[i + 1], dx = x1 - x0;
        if (!intersects(x0, x1, xmin, xmax)) continue;

        auto cyz = Bool::common(shape, Builder::plane(float3(x0, 0, 0), float3(1, 0, 0)));
        if (BOPTools_AlgoTools3D::IsEmptyShape(cyz)) continue;

        auto byz = Shape::bound(cyz);
        for (int j = 0; j < ny - 1; j ++) {
            auto y0 = ys[j], y1 = ys[j + 1], dy = y1 - y0;
            if (!intersects(y0, y1, byz.ymin, byz.ymax)) continue;

            auto cy = Bool::common(cyz, Builder::box(float3(x0, y0, zmin), float3(x1, y1, zmax)));
            if (BOPTools_AlgoTools3D::IsEmptyShape(cy)) continue;

            auto by = Shape::bound(cy);
            for (int k = 0; k < nz - 1; k ++) {
                auto z0 = zs[k], z1 = zs[k + 1], dz = z1 - z0;
                if (!intersects(z0, z1, by.zmin, by.zmax)) continue;

                auto idx = i + j * nx + k * nx * ny;
                auto calcNorm = false;

                auto surf = Bool::common(cy, Builder::box(float3(x0, y0, z0), float3(x1, y1, z1)));
                auto area = areaArr[idx] = getSurfaceProps(surf).Mass() / dy / dz;
                if (area > 0.001 && area < 0.999) {
                    surfMap[idx] = surf;
                    calcNorm = true;
                }

                auto edge = Bool::common(cy, Builder::line(float3(x0, y0, z0), float3(x0, y1, z0)));
                auto length = lengthArr[idx] = getLinearProps(edge).Mass() / dy;
                if (length > 0.01 && length < 0.99) {
                    edgeMap[idx] = edge;
                    calcNorm = true;
                }

                if (calcNorm) {
                    calcWireNormal(bound_type { x0, y0, z0, x1, y1, z1 }, faceBounds, 0);
                }
            }
        }
    }
}

void Mesher::MeshY() {
    auto &areaArr = sy, &lengthArr = lz;
    auto &surfMap = msy, &edgeMap = mlz;

    for (int j = 0; j < ny - 1; j ++) {
        auto y0 = ys[j], y1 = ys[j + 1], dy = y1 - y0;
        if (!intersects(y0, y1, ymin, ymax)) continue;

        auto czx = Bool::common(shape, Builder::plane(float3(0, y0, 0), float3(0, 1, 0)));
        if (BOPTools_AlgoTools3D::IsEmptyShape(czx)) continue;

        auto bzx = Shape::bound(czx);
        for (int k = 0; k < nz - 1; k ++) {
            auto z0 = zs[k], z1 = zs[k + 1], dz = z1 - z0;
            if (!intersects(z0, z1, bzx.zmin, bzx.zmax)) continue;

            auto cz = Bool::common(czx, Builder::box(float3(xmin, y0, z0), float3(xmax, y1, z1)));
            if (BOPTools_AlgoTools3D::IsEmptyShape(cz)) continue;

            auto bz = Shape::bound(cz);
            for (int i = 0; i < nx - 1; i ++) {
                auto x0 = xs[i], x1 = xs[i + 1], dx = x1 - x0;
                if (!intersects(x0, x1, bz.xmin, bz.xmax)) continue;

                auto idx = i + j * nx + k * nx * ny;
                auto calcNorm = false;

                auto surf = Bool::common(cz, Builder::box(float3(x0, y0, z0), float3(x1, y1, z1)));
                auto area = areaArr[idx] = getSurfaceProps(surf).Mass() / dz / dx;
                if (area > 0.001 && area < 0.999) {
                    surfMap[idx] = surf;
                    calcNorm = true;
                }

                auto edge = Bool::common(cz, Builder::line(float3(x0, y0, z0), float3(x0, y0, z1)));
                auto length = lengthArr[idx] = getLinearProps(edge).Mass() / dz;
                if (length > 0.01 && length < 0.99) {
                    edgeMap[idx] = edge;
                    calcNorm = true;
                }

                if (calcNorm) {
                    calcWireNormal(bound_type { x0, y0, z0, x1, y1, z1 }, faceBounds, 1);
                }
            }
        }
    }
}

void Mesher::MeshZ() {
    auto &areaArr = sz, &lengthArr = lx;
    auto &surfMap = msz, &edgeMap = mlx;

    for (int k = 0; k < nz - 1; k ++) {
        auto z0 = zs[k], z1 = zs[k + 1], dz = z1 - z0;
        if (!intersects(z0, z1, zmin, zmax)) continue;

        auto cxy = Bool::common(shape, Builder::plane(float3(0, 0, z0), float3(0, 0, 1)));
        if (BOPTools_AlgoTools3D::IsEmptyShape(cxy)) continue;

        auto bxy = Shape::bound(cxy);
        for (int i = 0; i < nx - 1; i ++) {
            auto x0 = xs[i], x1 = xs[i + 1], dx = x1 - x0;
            if (!intersects(x0, x1, bxy.xmin, bxy.xmax)) continue;

            auto cx = Bool::common(cxy, Builder::box(float3(x0, ymin, z0), float3(x1, ymax, z1)));
            if (BOPTools_AlgoTools3D::IsEmptyShape(cx)) continue;

            auto bx = Shape::bound(cx);
            for (int j = 0; j < ny - 1; j ++) {
                auto y0 = ys[j], y1 = ys[j + 1], dy = y1 - y0;
                if (!intersects(y0, y1, bx.ymin, bx.ymax)) continue;

                auto idx = i + j * nx + k * nx * ny;
                auto calcNorm = false;

                auto surf = Bool::common(cx, Builder::box(float3(x0, y0, z0), float3(x1, y1, z1)));
                auto area = areaArr[idx] = getSurfaceProps(surf).Mass() / dx / dy;
                if (area > 0.001 && area < 0.999) {
                    surfMap[idx] = surf;
                    calcNorm = true;
                }

                auto edge = Bool::common(cx, Builder::line(float3(x0, y0, z0), float3(x1, y0, z0)));
                auto length = lengthArr[idx] = getLinearProps(edge).Mass() / dx;
                if (length > 0.01 && length < 0.99) {
                    edgeMap[idx] = edge;
                    calcNorm = true;
                }

                if (calcNorm) {
                    calcWireNormal(bound_type { x0, y0, z0, x1, y1, z1 }, faceBounds, 2);
                }
            }
        }
    }
}
