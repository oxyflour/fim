#include <ctpl_stl.h>

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

#include "cuda.h"
#include "occ.h"
#include "utils.h"
#include "utils/check.h"

using namespace std;
using namespace occ;

// https://stackoverflow.com/questions/64912439/preventing-opencascade-to-write-to-console
//Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));

TopoDS_Shape Step::load(string file) {
    STEPControl_Reader reader;
    auto stat = reader.ReadFile(file.c_str());
    reader.TransferRoot();
    CHECK(stat == IFSelect_RetDone, "load from file " + file + " failed");
    return reader.Shape();
}

void Step::save(string file, TopoDS_Shape &shape) {
    STEPControl_Writer writer;
    auto stat = writer.Transfer(shape, STEPControl_StepModelType::STEPControl_AsIs);
    CHECK(stat == IFSelect_RetDone, "save to file " + file + " failed");
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

Mesher::Mesher(grid::Grid &grid, TopoDS_Shape &shape, float unit) {
    xs = grid.xs; ys = grid.ys; zs = grid.zs;

    nx = xs.size(); ny = ys.size(); nz = zs.size();
    for (int i = 0; i < nx; i ++) xs[i] /= unit;
    for (int i = 0; i < ny; i ++) ys[i] /= unit;
    for (int i = 0; i < nz; i ++) zs[i] /= unit;

    this->shape = TopoDS_Shape(shape);
    faces = Builder::component(Shape::find(shape, TopAbs_FACE));
    auto bound = Shape::bound(shape);
    xmin = bound.xmin; ymin = bound.ymin; zmin = bound.zmin;
    xmax = bound.xmax; ymax = bound.ymax; zmax = bound.zmax;

    auto nxyz = nx * ny * nz;
    sx.resize(nxyz); std::fill(sx.begin(), sx.end(), 0.f);
    sy.resize(nxyz); std::fill(sy.begin(), sy.end(), 0.f);
    sz.resize(nxyz); std::fill(sz.begin(), sz.end(), 0.f);
    lx.resize(nxyz); std::fill(lx.begin(), lx.end(), 0.f);
    ly.resize(nxyz); std::fill(ly.begin(), ly.end(), 0.f);
    lz.resize(nxyz); std::fill(lz.begin(), lz.end(), 0.f);

    auto pool = ctpl::thread_pool(thread::hardware_concurrency());
    for (int i = 0; i < nx - 1; i ++) {
        pool.push([&, i](int id) { MeshX(i, i + 1); });
    }
    pool.stop(true);
}

void Mesher::Save(string file) {
    vector<TopoDS_Shape> shapes;
    for (auto [id, shape] : msx) shapes.push_back(shape);
    for (auto [id, shape] : msy) shapes.push_back(shape);
    for (auto [id, shape] : msz) shapes.push_back(shape);
    for (auto [id, shape] : mlx) shapes.push_back(shape);
    for (auto [id, shape] : mly) shapes.push_back(shape);
    for (auto [id, shape] : mlz) shapes.push_back(shape);
    Step::save(file, Builder::component(shapes));
}

const int METAL = 1000;

void MergeCell(grid::Grid &grid, int i, int j, int k,
        vector<int> &idx, vector<Mesher> &meshers, vector<int> &priority) {
/*
    auto n = grid.GetIndex(i, j, k);
    {
        auto i = idx[idx.size() - 1];
        auto &mesher = meshers[i];
        if (priority[i] > METAL) {
            if (mesher.sx[n] > 0.999) {
                // filled with metal
            }
            if (mesher.sy[n] > 0.999) {
                // filled with metal
            }
            if (mesher.sz[n] > 0.999) {
                // filled with metal
            }
        }
    }
    int lastPriority = 0;
    auto material = std::map<int, Cell>();
    for (auto id : idx) {
        auto &mesher = meshers[id];
        // merge x
        if (mesher.msx.count(n)) {
            if (priority[id] > lastPriority) {
            } else {
            }
        }
        // merge y
        // merge z
        lastPriority = priority[id];
    }
 */
}

void Mesher::Merge(grid::Grid &grid, vector<Mesher> &meshers, vector<int> &priority) {
    auto idx = utils::range((int) meshers.size());
    std::sort(idx.begin(), idx.end(), [&](int i, int j) { return priority[i] - priority[j]; });

    auto pool = ctpl::thread_pool(thread::hardware_concurrency());
    for (int i = 0; i < grid.nx - 1; i ++) {
        for (int j = 0; j < grid.ny - 1; j ++) {
            for (int k = 0; k < grid.nz - 1; k ++) {
                pool.push([&](int id) { MergeCell(grid, i, j, k, idx, meshers, priority); });
            }
        }
    }
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

void calcWireNormal(bound_type &bound, vector<TopoDS_Shape> &faces, int dir) {
    auto xmin = bound.xmin, ymin = bound.ymin, zmin = bound.zmin,
        xmax = bound.xmax, ymax = bound.ymax, zmax = bound.zmax;
    for (auto &fc : faces) {
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
                gp_Pnt2d pt2;
                gp_Pnt pos;
                BOPTools_AlgoTools3D::PointNearEdge(edge, face, 0.5, 0, pt2, pos);
                gp_Dir dir;
                BOPTools_AlgoTools3D::GetNormalToFaceOnEdge(edge, face, 0.5, dir);
                // TODO: output
                printf("%f, %f, %f / %f, %f, %f\n", pos.X(), pos.Y(), pos.Z(), dir.X(), dir.Y(), dir.Z());
            }
        }
    }
}

void Mesher::MeshX(int i0, int i1) {
    auto &areaArr = sx, &lengthArr = ly;
    auto &surfMap = msx, &edgeMap = mly;

    for (int i = max(i0, 0); i < min(i1, nx - 1); i ++) {
        auto x0 = xs[i], x1 = xs[i + 1], dx = x1 - x0;
        if (!intersects(x0, x1, xmin, xmax)) continue;

        auto box = Builder::box(make_float3(x0, ymin, zmin), make_float3(x1, ymax, zmax));
        auto cyz = Bool::common(shape, Builder::plane(make_float3(x0, 0, 0), make_float3(1, 0, 0)));
        auto myz = Bool::common(faces, box);
        if (BOPTools_AlgoTools3D::IsEmptyShape(cyz)) continue;

        auto byz = Shape::bound(cyz);
        for (int j = 0; j < ny - 1; j ++) {
            auto y0 = ys[j], y1 = ys[j + 1], dy = y1 - y0;
            if (!intersects(y0, y1, byz.ymin, byz.ymax)) continue;

            auto box = Builder::box(make_float3(x0, y0, zmin), make_float3(x1, y1, zmax));
            auto cy = Bool::common(cyz, box);
            auto my = Bool::common(myz, box);
            if (BOPTools_AlgoTools3D::IsEmptyShape(cy)) continue;

            auto by = Shape::bound(cy);
            for (int k = 0; k < nz - 1; k ++) {
                auto z0 = zs[k], z1 = zs[k + 1], dz = z1 - z0;
                if (!intersects(z0, z1, by.zmin, by.zmax)) continue;

                auto idx = i + j * nx + k * nx * ny;

                auto box = Builder::box(make_float3(x0, y0, z0), make_float3(x1, y1, z1));
                auto surf = Bool::common(cy, box);
                auto area = areaArr[idx] = getSurfaceProps(surf).Mass() / dy / dz;
                if (area > 0.001 && area < 0.999) {
                    lock_guard guard(lock);
                    surfMap[idx] = surf;
                }

                auto line = Builder::line(make_float3(x0, y0, z0), make_float3(x0, y1, z0));
                auto edge = Bool::common(cy, line);
                auto length = lengthArr[idx] = getLinearProps(edge).Mass() / dy;
                if (length > 0.01 && length < 0.99) {
                    lock_guard guard(lock);
                    edgeMap[idx] = edge;
                }

                auto faces = Shape::find(Bool::common(my, box), TopAbs_FACE);
                if (faces.size()) {
                    calcWireNormal(bound_type { x0, y0, z0, x1, y1, z1 }, faces, 0);
                }
            }
        }
    }
}

void Mesher::MeshY(int j0, int j1) {
    auto &areaArr = sy, &lengthArr = lz;
    auto &surfMap = msy, &edgeMap = mlz;

    for (int j = max(0, j0); j < min(ny - 1, j1); j ++) {
        auto y0 = ys[j], y1 = ys[j + 1], dy = y1 - y0;
        if (!intersects(y0, y1, ymin, ymax)) continue;

        auto box = Builder::box(make_float3(xmin, y0, zmin), make_float3(xmax, y1, zmax));
        auto czx = Bool::common(shape, Builder::plane(make_float3(0, y0, 0), make_float3(0, 1, 0)));
        auto mzx = Bool::common(faces, box);
        if (BOPTools_AlgoTools3D::IsEmptyShape(czx)) continue;

        auto bzx = Shape::bound(czx);
        for (int k = 0; k < nz - 1; k ++) {
            auto z0 = zs[k], z1 = zs[k + 1], dz = z1 - z0;
            if (!intersects(z0, z1, bzx.zmin, bzx.zmax)) continue;

            auto box = Builder::box(make_float3(xmin, y0, z0), make_float3(xmax, y1, z1));
            auto cz = Bool::common(czx, box);
            auto mz = Bool::common(mzx, box);
            if (BOPTools_AlgoTools3D::IsEmptyShape(cz)) continue;

            auto bz = Shape::bound(cz);
            for (int i = 0; i < nx - 1; i ++) {
                auto x0 = xs[i], x1 = xs[i + 1], dx = x1 - x0;
                if (!intersects(x0, x1, bz.xmin, bz.xmax)) continue;

                auto idx = i + j * nx + k * nx * ny;

                auto box = Builder::box(make_float3(x0, y0, z0), make_float3(x1, y1, z1));
                auto surf = Bool::common(cz, box);
                auto area = areaArr[idx] = getSurfaceProps(surf).Mass() / dz / dx;
                if (area > 0.001 && area < 0.999) {
                    surfMap[idx] = surf;
                }

                auto line = Builder::line(make_float3(x0, y0, z0), make_float3(x0, y0, z1));
                auto edge = Bool::common(cz, line);
                auto length = lengthArr[idx] = getLinearProps(edge).Mass() / dz;
                if (length > 0.01 && length < 0.99) {
                    edgeMap[idx] = edge;
                }

                auto faces = Shape::find(Bool::common(mz, box), TopAbs_FACE);
                if (faces.size()) {
                    calcWireNormal(bound_type { x0, y0, z0, x1, y1, z1 }, faces, 1);
                }
            }
        }
    }
}

void Mesher::MeshZ(int j0, int j1) {
    auto &areaArr = sz, &lengthArr = lx;
    auto &surfMap = msz, &edgeMap = mlx;

    for (int k = max(0, j0); k < min(nz - 1, j1); k ++) {
        auto z0 = zs[k], z1 = zs[k + 1], dz = z1 - z0;
        if (!intersects(z0, z1, zmin, zmax)) continue;

        auto box = Builder::box(make_float3(xmin, ymin, z0), make_float3(xmax, ymax, z1));
        auto cxy = Bool::common(shape, Builder::plane(make_float3(0, 0, z0), make_float3(0, 0, 1)));
        auto mxy = Bool::common(faces, box);
        if (BOPTools_AlgoTools3D::IsEmptyShape(cxy)) continue;

        auto bxy = Shape::bound(cxy);
        for (int i = 0; i < nx - 1; i ++) {
            auto x0 = xs[i], x1 = xs[i + 1], dx = x1 - x0;
            if (!intersects(x0, x1, bxy.xmin, bxy.xmax)) continue;

            auto box = Builder::box(make_float3(x0, ymin, z0), make_float3(x1, ymax, z1));
            auto cx = Bool::common(cxy, box);
            auto mx = Bool::common(mxy, box);
            if (BOPTools_AlgoTools3D::IsEmptyShape(cx)) continue;

            auto bx = Shape::bound(cx);
            for (int j = 0; j < ny - 1; j ++) {
                auto y0 = ys[j], y1 = ys[j + 1], dy = y1 - y0;
                if (!intersects(y0, y1, bx.ymin, bx.ymax)) continue;

                auto idx = i + j * nx + k * nx * ny;

                auto box = Builder::box(make_float3(x0, y0, z0), make_float3(x1, y1, z1));
                auto surf = Bool::common(cx, box);
                auto area = areaArr[idx] = getSurfaceProps(surf).Mass() / dx / dy;
                if (area > 0.001 && area < 0.999) {
                    surfMap[idx] = surf;
                }

                auto line = Builder::line(make_float3(x0, y0, z0), make_float3(x1, y0, z0));
                auto edge = Bool::common(cx, line);
                auto length = lengthArr[idx] = getLinearProps(edge).Mass() / dx;
                if (length > 0.01 && length < 0.99) {
                    edgeMap[idx] = edge;
                }

                auto faces = Shape::find(Bool::common(mx, box), TopAbs_FACE);
                if (faces.size()) {
                    calcWireNormal(bound_type { x0, y0, z0, x1, y1, z1 }, faces, 2);
                }
            }
        }
    }
}
