//
// Created by aburahman10 on 28/02/23.
//

#include "Line2D.h"

Line2D::Line2D(const Point2D &p1, const Point2D &p2, const Colour &color) : p1(p1), p2(p2), color(color) {}

Line2D::Line2D(const Point2D &p1, const Point2D &p2, const Colour &color, double z1, double z2) : p1(p1), p2(p2),
                                                                                                  color(color), z1(z1),
                                                                                                  z2(z2) {}
