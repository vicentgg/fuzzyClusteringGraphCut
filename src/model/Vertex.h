#ifndef MESHSEGMENTATION_VERTEX_H
#define MESHSEGMENTATION_VERTEX_H


#include "Point.h"

class Vertex
{
public:
    Point p;

    Vertex() = default;

    Vertex(double xx, double y, double z);

    Vertex(const Vertex& an);
};

Vertex::Vertex(double xx, double y, double z)
{
    p.x[0] = xx;
    p.x[1] = y;
    p.x[2] = z;
}

Vertex::Vertex(const Vertex &an)
{
    p.x[0] = an.p.x[0];
    p.x[1] = an.p.x[1];
    p.x[2] = an.p.x[2];
}

#endif //MESHSEGMENTATION_VERTEX_H
