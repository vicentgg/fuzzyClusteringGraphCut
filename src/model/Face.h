
#ifndef MESHSEGMENTATION_FACE_H
#define MESHSEGMENTATION_FACE_H

#include <vector>
#include "Point.h"
#include "Indices.h"
#include "../util/Util.h"
/** Structure to store neighbor`s faces*/
struct DualEdge
{
    int face;
    float weight;
    float angle;
    float angle_dist, geo_dist;

    DualEdge(int f, float a_dist, float g_dist, float ang)
            : face(f),
              angle_dist(a_dist),
              geo_dist(g_dist),
              angle(ang),
              weight(0.0f)
    {}

    friend std::ostream &operator<<(std::ostream &os, const DualEdge &de)
    {
        return os << de.face << ", " << de.weight << ", " << de.angle_dist << ", "
                  << de.geo_dist;
    }
};

class Face
{
public:

    Indices indices;
    std::vector<Face*> neighbors;
    std::vector<DualEdge> dedges;  // Neighbor faces for Graph Cut
    std::vector<double> Angs; // 各顶点的内角

    Point normal;
    Point center;

    int label = 0;
    bool marked = false;

    Face() = default;

    Face(int v1, int v2, int v3);

    Face(Indices index, const Point& px, const Point& py, const Point& pz);

    Face(const Face &an);

    friend bool operator==(const Face &a, const Face &b);

    friend bool operator<(const Face &a, const Face &b);
};

Face::Face(int v1, int v2, int v3)
{
    indices.x = v1;
    indices.y = v2;
    indices.z = v3;
}

bool operator==(const Face &a, const Face &b)
{
    if (a.indices.x == b.indices.x && a.indices.y == b.indices.y && a.indices.z == b.indices.z)
        return true;
    return false;
}

bool operator<(const Face &a, const Face &b)
{
    if (a.indices.x < b.indices.x && a.indices.y < b.indices.y && a.indices.z < b.indices.z)
        return true;
    return false;
}

Face::Face(const Face &an)
{
    indices.x = an.indices.x;
    indices.y = an.indices.y;
    indices.z = an.indices.z;

    for (auto& f : an.neighbors)
        neighbors.push_back(f);

    normal = Point(an.normal);
    center = Point(an.center);
    label = an.label;
}

Face::Face(Indices index, const Point& px, const Point& py, const Point& pz)
{
    indices = index;
    center = (px + py + pz) / 3;
    // std::cout << "center: " << center << std::endl;
    normal = util::cross(py - px, pz - px);
    util::normalizeV(normal);
}

#endif //MESHSEGMENTATION_FACE_H
