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
#include "ZBuffer.h"

using namespace std;
using namespace img;
using namespace LParser;
using namespace ini;

class Functies
{
public:

    //DRAW 2D LINES
    static EasyImage draw2DLines(const Lines2D &lines, const int size, vector<double> backgroundColor, string type);

    static EasyImage draw2DLines(const Lines2D &lines, Figures3D &figures3D, const int size, vector<double> backgroundColor);

    //L-SYSTEM
    static Lines2D drawLSystem(const LSystem2D &l_system, vector<double> color);

    static void leesString(double starting_angle, double angle, Lines2D &lines, double &x, double &y, string string1, vector<double> color);

    static void tekenReplace(string &initiator, vector<pair<char,string>> replacements, unsigned int iterations);

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

    static Lines2D pasFigure(Figures3D &figures3D, const Configuration &configuration);

    // 3D LICHAMEN
    static Figure createCube();

    static Figure createTetrahedron();

    static Figure createOctahedron();

    static Figure createIcosahedron();

    static Figure createDodecahedron();

    static Figure createSphere(const double radius, const int n);

    static Figure createCone(const int n, const double h);

    static Figure createCylinder(const int n, const double h);

    static Figure createTorus(const double r, const double R, const int n, const int m);

    ///°°°°°°/// 3D L-SYSTEM
    static Figure drawLSystem3D(const LSystem3D &lSystem3D);

    static void leesString(const double angle, string string1, Figure& figure, const LSystem3D &lSystem3D);
};

#endif //ENGINE_FUNCTIES_H
