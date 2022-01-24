
#ifndef MESHSEGMENTATION_POINT_H
#define MESHSEGMENTATION_POINT_H

#include <iostream>

class Point
{
public:
    double x[3]{};

    Point(double xx, double yy, double zz);

    Point(const Point &p);

    Point operator+(const Point &p) const;

    Point operator-(const Point &p) const;

    Point operator*(double d);

    Point operator/(double d);

    Point() = default;

    friend std::ostream& operator<<(std::ostream& stream, const Point& p);
};

Point::Point(double xx, double yy, double zz)
{
    x[0] = xx; x[1] = yy; x[2] = zz;
}

Point Point::operator+(const Point &p) const
{
    Point res;
    for (int i = 0; i < 3; ++i)
        res.x[i] = x[i] + p.x[i];
    return res;
}

Point Point::operator-(const Point &p) const
{
    Point res;
    for (int i = 0; i < 3; ++i)
        res.x[i] = x[i] - p.x[i];
    return res;
}

Point Point::operator*(double d)
{
    Point res;
    for (int i = 0; i < 3; ++i)
        res.x[i] = x[i] * d;
    return res;
}

Point Point::operator/(double d)
{
    Point res;
    for (int i = 0; i < 3; ++i)
        res.x[i] = x[i] / d;
    return res;
}

Point::Point(const Point &p)
{
    x[0] = p.x[0]; x[1] = p.x[1]; x[2] = p.x[2];
}

std::ostream &operator<<(std::ostream &stream, const Point &p)
{
    stream << "(" << p.x[0] << "," << p.x[1] << "," << p.x[2] << ")";

    return stream;
}


#endif //MESHSEGMENTATION_POINT_H
