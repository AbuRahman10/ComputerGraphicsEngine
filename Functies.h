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

using namespace std;
using namespace img;
using namespace LParser;
using namespace ini;

class Functies
{
public:

    static EasyImage draw2DLines(const Lines2D &lines, const int size, vector<double> lineColor, vector<double> backgroundColor);

    static Lines2D drawLSystem(const LSystem2D &l_system);

    static void leesString(double starting_angle, double angle, Lines2D &lines, double &x, double &y, string string1);

    static void tekenReplace(string &initiator, double starting_angle, double angle, Lines2D &lines, double x, double y, vector<pair<char,string>> replacements, unsigned int iterations);
};

#endif //ENGINE_FUNCTIES_H
