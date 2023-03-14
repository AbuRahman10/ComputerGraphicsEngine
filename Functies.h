//
// Created by aburahman10 on 28/02/23.
//

#ifndef ENGINE_FUNCTIES_H
#define ENGINE_FUNCTIES_H

#include "iostream"
#include "easy_image.h"
#include "Line2D.h"
#include "Colour.h"
#include "Point2D.h"
#include "l_parser.h"
#include "ini_configuration.h"
#include "vector3d.h"
#include "Figure.h"

using namespace std;
using namespace img;
using namespace LParser;
using namespace ini;

class Functies
{
public:

    //DRAW 2D LINES
    static EasyImage draw2DLines(const Lines2D &lines, const int size, vector<double> lineColor, vector<double> backgroundColor);

    //L-SYSTEM
    static Lines2D drawLSystem(const LSystem2D &l_system);

    static void leesString(double starting_angle, double angle, Lines2D &lines, double &x, double &y, string string1);

    static void tekenReplace(string &initiator, double starting_angle, double angle, Lines2D &lines, double x, double y, vector<pair<char,string>> replacements, unsigned int iterations);

    //DRAW 3D LINE
    static Matrix scaleFigure(const double scale);

    static Matrix rotateX(const double angle);

    static Matrix rotateY(const double angle);

    static Matrix rotateZ(const double angle);

    static Matrix translate(const Vector3D &vector);

    static void applyTransformation(Figure &figure, const Matrix &matrix);

    static Matrix eyePointTrans(const Vector3D &eyepoint);

    static void toPolar(const Vector3D &point, double &theta, double &phi, double &r);

    static void applyTransformation(Figures3D &figure, const Matrix &matrix);

    static Lines2D doProjection(const Figures3D &figures3D);

    static Point2D doProjection(const Vector3D &point, const double d);

    static void pasFigure(Figures3D &figures3D, const Configuration &configuration, Colour colour);

    static Lines2D omzetDimensie_3D_2D(Figures3D &figures3D, const Configuration &configuration);
};

#endif //ENGINE_FUNCTIES_H
