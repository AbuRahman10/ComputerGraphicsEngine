//
// Created by aburahman10 on 13/03/23.
//

#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H

#include "iostream"
#include "vector"
#include "Point2D.h"
#include "Colour.h"
#include "vector3d.h"
#include "Face.h"

using namespace std;

class Figure
{
public:
    vector<Vector3D> points;
    vector<Face> faces;
    Colour color;
};

typedef vector<Figure> Figures3D;


#endif //ENGINE_FIGURE_H
