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
#include "stack"

using namespace std;
using namespace img;
using namespace LParser;

EasyImage Functies::draw2DLines(const Lines2D &lines, const int size, vector<double> lineColor, vector<double> backgroundColor)
{
    // x.min; x.max; y.min; y.max;

    double Xmin = lines[0].p1.x;
    double Ymin = lines[0].p1.y;

    double Xmax = lines[0].p1.x;
    double Ymax = lines[0].p1.y;

    for(const Line2D i : lines)
    {
        if(i.p1.x < Xmin)
        {
            Xmin = i.p1.x;
        }
        if(i.p2.x < Xmin)
        {
            Xmin = i.p2.x;
        }
        if(i.p1.y < Ymin)
        {
            Ymin = i.p1.y;
        }
        if(i.p2.y < Ymin)
        {
            Ymin = i.p2.y;
        }
        if(i.p1.x > Xmax)
        {
            Xmax = i.p1.x;
        }
        if(i.p2.x > Xmax)
        {
            Xmax = i.p2.x;
        }
        if(i.p1.y > Ymax)
        {
            Ymax = i.p1.y;
        }
        if(i.p2.y > Ymax)
        {
            Ymax = i.p2.y;
        }
    }

    // GROOTTE VAN DE IMAGE
    double x_range = Xmax - Xmin;
    double y_range = Ymax - Ymin;

    double max1 = max(x_range,y_range);

    double deling_x_range = (x_range) / max1;
    double deling_y_range = (y_range) / max1;
    double Image_x = size * deling_x_range;
    double Image_y = size * deling_y_range;

    // LIJNTEKENING SCHALEN
    double d = 0.95 * (Image_x/ x_range); // schaalfactor d;

    Lines2D new_lines = lines;

    int size_newLine = new_lines.size();

    for (int i = 0; i < lines.size(); i++)
    {
        new_lines[i].p1.x *= d;
        new_lines[i].p2.x *= d;

        new_lines[i].p1.y *= d;
        new_lines[i].p2.y *= d;
    }

    // LIJNTEKENING VERSCHUIVEN

    double DC_x = d * ((Xmin + Xmax) / 2);
    double DC_y = d * ((Ymin + Ymax) / 2);

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

    for (unsigned int i = 0; i < lround(Image_y); i++)
    {
        for (unsigned int j = 0; j < ::lround(Image_x); j++)
        {
            image(j, i).red = lround(backgroundColor[0] * 255);
            image(j, i).green = lround(backgroundColor[1] * 255);
            image(j, i).blue = lround(backgroundColor[2] * 255);
        }
    }

    Color color(lround(lineColor[0]*255), lround(lineColor[1]*255), lround(lineColor[2]*255));
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
    string initiator = l_system.get_initiator();
    unsigned int iterations = l_system.get_nr_iterations();
    double starting_angle = l_system.get_starting_angle();
    double angle = l_system.get_angle();

    //REPLACEMENTS
    vector<pair<char,string>> replacements;

    for (char i : alfabet)
    {
        string str = l_system.get_replacement(i);
        pair<char,string> letter(i,str);
        replacements.push_back(letter);
    }

    starting_angle *= (pi/180);
    angle *= (pi/180);

    double x = 0;
    double y = 0;

    tekenReplace(initiator,starting_angle,angle,lines,x,y,replacements, iterations);
    leesString(starting_angle,angle,lines,x,y,initiator);

    return lines;
}

void Functies::tekenReplace(string &initiator, double starting_angle, double angle, Lines2D &lines, double x, double y, vector<pair<char,string>> replacements, unsigned int iterations)
{

    for (int it = 0; it < iterations; it++)
    {
        string nieuw_initiator;
        for (char i : initiator)
        {
            if (i == '+')
            {
                nieuw_initiator += i;
            }
            else if (i == '-')
            {
                nieuw_initiator += i;
            }
            else if (i == '(')
            {
                nieuw_initiator += i;
            }
            else if (i == ')')
            {
                nieuw_initiator += i;
            }
            else
            {
                for (pair<char,string> alfabet : replacements)
                {
                    if (alfabet.first == i)
                    {
                        nieuw_initiator += alfabet.second;
                    }
                }
            }
        }
        initiator = nieuw_initiator;
    }
}

void Functies::leesString(double starting_angle, double angle, Lines2D &lines, double &x, double &y, string string1)
{
    stack<pair<pair<double,double>,double>> myStack;

    for (char k : string1)
    {
        if (k == '-')
        {
            starting_angle -= angle;
        }
        else if (k == '+')
        {
            starting_angle += angle;
        }
        else if (k == '(')
        {
            pair<double,double> position(x,y);
            pair<pair<double,double>,double> pos_EN_hoek(position,starting_angle);

            myStack.push(pos_EN_hoek);
        }
        else if (k == ')')
        {
            pair<pair<double,double>,double> pos_EN_hoek = myStack.top();

            x = pos_EN_hoek.first.first;
            y = pos_EN_hoek.first.second;

            starting_angle = pos_EN_hoek.second;

            myStack.pop();
        }
        else
        {
            Point2D point1(x,y);
            double cosa = cos(starting_angle);
            double sina = sin(starting_angle);
            x += cosa;
            y += sina;
            Point2D point2(x,y);
            Colour color(0,0,0);
            Line2D line(point1,point2,color);
            lines.push_back(line);
        }
    }
}
