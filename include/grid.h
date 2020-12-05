#include <vector>

#ifndef GRID_H
#define GRID_H

typedef struct float3 {
    float x, y, z;
    float3() { }
    template <typename A, typename B, typename C> float3(A x, B y, C z) : x(x), y(y), z(z) { }
    float3 operator- (const float3& first);
    float length();
} float3;

typedef struct int3 {
    int x, y, z;
    int3 operator- (const int3& first);
    int3 operator+ (const int3& first);
    int sum();
    int3 abs();
} int3;

typedef struct Grid {
    std::vector<double> xs, ys, zs;
    float3 At(int3 idx);
    int GetFlatIndex(int3 idx, int dir);
    int3 FindIndex(float3 pos, float epsi);
    std::vector<int3> ParsePort(float3 src, float3 dst, float epsi);
} Grid;

#endif
