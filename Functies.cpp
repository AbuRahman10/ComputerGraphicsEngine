//
// Created by aburahman10 on 28/02/23.
//

#include "Functies.h"
#include <fstream>
#include "l_parser.h"

#include "cmath"
#include "easy_image.h"

#include <iostream>
#include <string>

#include "Line2D.h"
#include "Point2D.h"
#include "Colour.h"

using namespace std;
using namespace img;
using namespace LParser;

EasyImage Functies::draw2DLines(const Lines2D &lines, const int size)
{
    // x.min; x.max; y.min; y.max;

    double x_min = lines[0].p1.x;
    double y_min = lines[0].p1.y;

    double x_max = lines[0].p1.x;
    double y_max = lines[0].p1.y;

    for (int i = 0; i < lines.size(); i++)
    {
        double x_P1 = lines[i].p1.x;
        double x_P2 = lines[i].p2.x;

        double y_P1 = lines[i].p1.y;
        double y_P2 = lines[i].p2.y;

        if (y_P1 < y_min)
        {
            y_min = y_P1;
        }
        if (y_P2 < y_P1)
        {
            y_min = y_P2;
        }

        if (x_P1 < x_min)
        {
            x_min = x_P1;
        }
        if (x_P2 < x_P1)
        {
            x_min = x_P2;
        }

        if (y_P1 > y_max)
        {
            y_max = y_P1;
        }
        if (y_P2 > y_P1)
        {
            y_max = y_P2;
        }

        if (x_P1 > x_max)
        {
            x_max = x_P1;
        }
        if (y_P2 > y_P1)
        {
            y_max = y_P2;
        }
    }

    // GROOTTE VAN DE IMAGE
    double x_range = x_max - x_min;
    double y_range = y_max - y_min;

    double max1 = max(x_range,y_range);

    double deling_x_range = (x_range) / max1;
    double deling_y_range = (y_range) / max1;
    double Image_x = size * deling_x_range;
    double Image_y = size * deling_y_range;

    // LIJNTEKENING SCHALEN
    double d = 0.95 * (Image_x/ x_range); // schaalfactor d;

    Lines2D new_lines = lines;
    for (int i = 0; i < lines.size(); i++)
    {
        new_lines[i].p1.x *= d;
        new_lines[i].p2.x *= d;

        new_lines[i].p1.y *= d;
        new_lines[i].p2.y *= d;
    }

    // LIJNTEKENING VERSCHUIVEN

    double DC_x = d * ((x_min + x_max) / 2);
    double DC_y = d * ((y_min + y_max) / 2);

    double d_x = (Image_x / 2) - DC_x;
    double d_y = (Image_y / 2) - DC_y;

    for (int i = 0; i < lines.size(); i++)
    {
        new_lines[i].p1.x += d_x;
        new_lines[i].p2.x += d_x;

        new_lines[i].p1.y += d_y;
        new_lines[i].p2.y += d_y;
    }

    EasyImage image(lround(Image_x),lround(Image_y));


    Color color(0, 0, 0);

    for (unsigned int i = 0; i < Image_y; i++)
    {
        for (unsigned int j = 0; j < Image_x; j++)
        {
            image(i, j).red = 255;
            image(i, j).green = 255;
            image(i, j).blue = 255;
        }
    }

    for (int i = 0; i < new_lines.size(); i++)
    {
        double xP1 = lround(new_lines[i].p1.x);
        double yP1 = lround(new_lines[i].p1.y);

        double xP2 = lround(new_lines[i].p2.x);
        double yP2 = lround(new_lines[i].p2.y);
        image.draw_line(xP1, yP1, xP2, yP2,color);
    }
    return image;
}

Lines2D Functies::drawLSystem(const LSystem2D &l_system)
{
    Lines2D lines; // De uiteindelijke resultaat

    double pi = 3.14159265358979323846;

    // DE COMPONENTEN
    set<char> alfabet = l_system.get_alphabet();
    string replacement = l_system.get_replacement(*alfabet.begin());
    string initiator = l_system.get_initiator();
    unsigned int iterations = l_system.get_nr_iterations();
    double starting_angle = l_system.get_starting_angle();
    double angle = l_system.get_angle();

    starting_angle *= (pi/180);
    angle *= (pi/180);

    double x = 0;
    double y = 0;

    for (char i : initiator)
    {
        if (i == '-')
        {
            starting_angle -= angle;
        }
        else if (i == '+')
        {
            starting_angle += angle;
        }
        else if (i == *initiator.begin())
        {
            Point2D point1(x,y);
            x += cos(starting_angle);
            y += sin(starting_angle);
            Point2D point2(x,y);
            Colour color(0,0,0);
            Line2D line(point1,point2,color);
            lines.push_back(line);
        }
    }
    return lines;
}