//
// Created by aburahman10 on 28/02/23.
//

#ifndef ENGINE_LINE2D_H
#define ENGINE_LINE2D_H

#include "Point2D.h"
#include "Colour.h"
#include "easy_image.h"
#include "vector"

using namespace std;
using namespace img;

class Line2D
{
public:
    Point2D p1;
    Point2D p2;
    Colour color;

    double z1;
    double z2;

    Line2D(const Point2D &p1, const Point2D &p2, const Colour &color);
};

using Lines2D = vector<Line2D>;

#endif //ENGINE_LINE2D_H
