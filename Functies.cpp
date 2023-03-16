//
// Created by aburahman10 on 28/02/23.
//

#include "Functies.h"
#include <fstream>
#include "l_parser.h"

#include "cmath"
#include "math.h"
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


    for (int i = 0; i < new_lines.size(); i++)
    {
        double xP1 = lround(new_lines[i].p1.x);
        double yP1 = lround(new_lines[i].p1.y);

        double xP2 = lround(new_lines[i].p2.x);
        double yP2 = lround(new_lines[i].p2.y);

        Color color(lround(new_lines[i].color.red*255), lround(new_lines[i].color.green*255), lround(new_lines[i].color.blue*255));

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

Matrix Functies::scaleFigure(const double scale)
{
    Matrix matrix;

    matrix(1,1) = scale;
    matrix(2,2) = scale;
    matrix(3,3) = scale;

    return matrix;
}

Matrix Functies::rotateX(const double angle)
{
    Matrix matrix;

    matrix(2,2) = cos(angle);
    matrix(3,3) = cos(angle);
    matrix(3,2) = -sin(angle);
    matrix(2,3) = sin(angle);

    return matrix;
}

Matrix Functies::rotateY(const double angle)
{
    Matrix matrix;

    matrix(1,1) = cos(angle);
    matrix(3,3) = cos(angle);
    matrix(3,1) = sin(angle);
    matrix(1,3) = -sin(angle);

    return matrix;
}

Matrix Functies::rotateZ(const double angle)
{
    Matrix matrix;

    matrix(1,1) = cos(angle);
    matrix(2,2) = cos(angle);
    matrix(2,1) = -sin(angle);
    matrix(1,2) = sin(angle);

    return matrix;
}

Matrix Functies::translate(const Vector3D &vector3D)
{
    Matrix matrix;

    matrix(4,1) = vector3D.x;
    matrix(4,2) = vector3D.y;
    matrix(4,3) = vector3D.z;

    return matrix;
}

void Functies::applyTransformation(Figure &figure, const Matrix &matrix)
{
    for (Vector3D &i : figure.points)
    {
        i *= matrix;
    }
}

Matrix Functies::eyePointTrans(const Vector3D &eyepoint)
{
    double r = 0;
    double theta = 0;
    double phi = 0;

    toPolar(eyepoint,theta,phi,r);

    Matrix matrix;

    matrix(1,1) = -sin(theta);
    matrix(1,2) = -cos(theta) * cos(phi);
    matrix(1,3) = cos(theta) * sin(phi);

    matrix(2,1) = cos(theta);
    matrix(2,2) = -sin(theta) * cos(phi);
    matrix(2,3) = sin(theta) * sin(phi);

    matrix(3,2) = sin(phi);
    matrix(3,3) = cos(phi);

    matrix(4,3) = -r;
    matrix(4,4) = 1;

    return matrix;
}

void Functies::toPolar(const Vector3D &point, double &theta, double &phi, double &r)
{
    r = sqrt(pow(point.x,2)+pow(point.y,2)+pow(point.z,2));
    theta = atan2(point.y,point.x);
    phi = acos(point.z/r);
}

void Functies::applyTransformation(Figures3D &figure, const Matrix &matrix)
{
    for (Figure &i : figure)
    {
        applyTransformation(i,matrix);
    }
}

Lines2D Functies::doProjection(const Figures3D &figures3D)
{
    Lines2D lines2D;

    for (Figure figure : figures3D)
    {
        for (Face face : figure.faces)
        {
            Point2D p1 = doProjection(figure.points[face.point_indexes[0]],1);
            Point2D p2 = doProjection(figure.points[face.point_indexes[1]],1);
            Line2D line2D(p1,p2,figure.color);
            lines2D.push_back(line2D);
        }
    }

    return lines2D;
}

Point2D Functies::doProjection(const Vector3D &point, const double d)
{
    double x = (d * point.x) / -point.z;
    double y = (d * point.y) / -point.z;

    Point2D point2D(x,y);
    return point2D;
}

Lines2D Functies::pasFigure(Figures3D &figures3D, const Configuration &configuration, Colour colour)
{
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();

    string figure = "Figure";
    for (int cnt = 0; cnt < nrFigures; cnt++)
    {
        Figure figure1;
        figure += to_string(cnt);

        int nrPoints = configuration[figure]["nrPoints"].as_int_or_die();
        int nrLines = configuration[figure]["nrLines"].as_int_or_die();
        vector<double> colors = configuration[figure]["color"].as_double_tuple_or_die();

        vector<Vector3D> points;
        string point = "point";
        for (int i = 0; i < nrPoints; i++)
        {
            point += to_string(i);
            vector<double> pts = configuration[figure][point].as_double_tuple_or_die();
            points.push_back(Vector3D::point(pts[0],pts[1],pts[2]));
            point = "point";
        }

        vector<Face> pts_collection;
        string line = "line";
        for (int i = 0; i < nrLines; i++)
        {
            line += to_string(i);
            vector<int> point_indexes = configuration[figure][line].as_int_tuple_or_die();
            Face face;
            face.point_indexes = point_indexes;
            pts_collection.push_back(face);
            line = "line";
        }

        figure1.faces = pts_collection;
        figure1.points = points;
        figure1.color.red = colors[0];
        figure1.color.green = colors[1];
        figure1.color.blue = colors[3];

        // OMZETTEN VAN 3D NAAR 2D

        vector<double> center = configuration[figure]["center"].as_double_tuple_or_die();
        double scale = configuration[figure]["scale"].as_double_or_die();
        double rotX = configuration[figure]["rotateX"].as_double_or_die();
        double rotY = configuration[figure]["rotateY"].as_double_or_die();
        double rotZ = configuration[figure]["rotateZ"].as_double_or_die();

        double pi = 3.14159265358979323846;

        // 1) MATRIX OPSTELLEN

        Matrix S,M,T;

        S = scaleFigure(scale);
        M = rotateX(rotX*(pi/180)) * rotateY(rotY*(pi/180)) * rotateZ(rotZ*(pi/180));
        T = translate(Vector3D::point(center[0],center[1],center[2]));

        Matrix omzetMatrix = S * M * T;
        applyTransformation(figure1,omzetMatrix);

        figure = "Figure";
        figures3D.push_back(figure1);
    }

    // 2) EYE_COORDINATEN

    vector<double> eyePoint = configuration["General"]["eye"].as_double_tuple_or_die();

    Matrix V = eyePointTrans(Vector3D::point(eyePoint[0],eyePoint[1],eyePoint[2]));
    applyTransformation(figures3D,V);

    // 3) 2D Lijntekening

    return doProjection(figures3D);
}