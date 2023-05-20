//
// Created by aburahman10 on 28/02/23.
//

#include "Functies.h"
#include <fstream>
#include "l_parser.h"
#include "cmath"
#include "easy_image.h"
#include <string>
#include <algorithm>
#include "Line2D.h"
#include "Point2D.h"
#include "Colour.h"
#include "stack"
#include "ZBuffer.h"
#include "limits"

using namespace std;
using namespace img;
using namespace LParser;

EasyImage Functies::draw2DLines(const Lines2D &lines, const int size, vector<double> backgroundColor, string type)
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

    ZBuffer zBuffer(lround(Image_x),lround(Image_y));

    for (int i = 0; i < new_lines.size(); i++)
    {
        double xP1 = lround(new_lines[i].p1.x);
        double yP1 = lround(new_lines[i].p1.y);

        double z1 = new_lines[i].z1;
        double z2 = new_lines[i].z2;

        double xP2 = lround(new_lines[i].p2.x);
        double yP2 = lround(new_lines[i].p2.y);

        Color color(lround(new_lines[i].color.red*255), lround(new_lines[i].color.green*255), lround(new_lines[i].color.blue*255));

        if (type == "ZBufferedWireframe")
        {
            zBuffer.draw_zbuf_line(image,xP1,yP1,z1,xP2,yP2,z2,color);
        }
        else
        {
            image.draw_line(xP1,yP1,xP2,yP2,color);
        }
    }
    return image;
}

EasyImage Functies::draw2DLines(const Lines2D &lines, Figures3D &figures3D, const int size, vector<double> backgroundColor)
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

    double DC_x = d * ((Xmin + Xmax) / 2);
    double DC_y = d * ((Ymin + Ymax) / 2);

    double d_x = (Image_x / 2) - DC_x;
    double d_y = (Image_y / 2) - DC_y;

    EasyImage image(lround(Image_x),lround(Image_y));

    // BACKGROUND-COLOR AANPASSEN
    for (unsigned int i = 0; i < lround(Image_y); i++)
    {
        for (unsigned int j = 0; j < ::lround(Image_x); j++)
        {
            image(j, i).red = lround(backgroundColor[0] * 255);
            image(j, i).green = lround(backgroundColor[1] * 255);
            image(j, i).blue = lround(backgroundColor[2] * 255);
        }
    }

    ZBuffer zBuffer(lround(Image_x),lround(Image_y));

    // TRIANGULATIE OP ELKE VLAK VAN ELKE FIGUUR
    for (Figure &figure : figures3D)
    {
        vector<Face> faces;
        for (Face &face : figure.faces)
        {
            vector<Face> triangles = zBuffer.triangulate(face);
            for (Face &f : triangles)
            {
                faces.push_back(f);
            }
        }
        figure.faces = faces;
    }

    for (Figure &figure : figures3D)
    {
        for (Face &face : figure.faces)
        {
            Color color(lround(figure.color.red*255), lround(figure.color.green*255), lround(figure.color.blue*255));
            zBuffer.draw_zbuf_triag(image,figure.points[face.point_indexes[0]],figure.points[face.point_indexes[1]],figure.points[face.point_indexes[2]],d,d_x,d_y,color);
        }
    }
    return image;
}

Lines2D Functies::drawLSystem(const LSystem2D &l_system, vector<double> color)
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

    tekenReplace(initiator,replacements, iterations);
    leesString(starting_angle,angle,lines,x,y,initiator, color);

    return lines;
}

void Functies::tekenReplace(string &initiator, vector<pair<char,string>> replacements, unsigned int iterations)
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
            else if (i == '&')
            {
                nieuw_initiator += i;
            }
            else if (i == '^')
            {
                nieuw_initiator += i;
            }
            else if (i == '/')
            {
                nieuw_initiator += i;
            }
            else if (i == '\\')
            {
                nieuw_initiator += i;
            }
            else if (i == '|')
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

void Functies::leesString(double starting_angle, double angle, Lines2D &lines, double &x, double &y, string string1, vector<double> color)
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
            Colour col(color[0],color[1],color[2]);
            Line2D line(point1,point2,col);
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
            for (int ind = 0; ind < face.point_indexes.size(); ind++)
            {
                Point2D p1;
                Point2D p2;
                double z1;
                double z2;
                if (ind == face.point_indexes.size() - 1)
                {
                    p1 = doProjection(figure.points[face.point_indexes[ind]],1);
                    p2 = doProjection(figure.points[face.point_indexes[0]],1);
                    z1 = figure.points[face.point_indexes[ind]].z;
                    z2 = figure.points[face.point_indexes[0]].z;
                }
                else
                {
                    p1 = doProjection(figure.points[face.point_indexes[ind]],1);
                    p2 = doProjection(figure.points[face.point_indexes[ind + 1]],1);
                    z1 = figure.points[face.point_indexes[ind]].z;
                    z2 = figure.points[face.point_indexes[ind + 1]].z;
                }
                Line2D line2D(p1,p2,figure.color,z1,z2);
                lines2D.push_back(line2D);
            }
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

Lines2D Functies::pasFigure(Figures3D &figures3D, const Configuration &configuration)
{
    int nrFigures = configuration["General"]["nrFigures"].as_int_or_die();
    string figure = "Figure";

    for (int cnt = 0; cnt < nrFigures; cnt++)
    {
        Figure figure1;
        figure += to_string(cnt);
        vector<double> colors = configuration[figure]["color"].as_double_tuple_or_die();
        string type = configuration[figure]["type"].as_string_or_die();

        if (type == "LineDrawing")
        {
            int nrPoints = configuration[figure]["nrPoints"].as_int_or_die();
            int nrLines = configuration[figure]["nrLines"].as_int_or_die();
            vector<Vector3D> points;
            vector<Face> pts_collection;
            string point = "point";
            for (int i = 0; i < nrPoints; i++)
            {
                point += to_string(i);
                vector<double> pts = configuration[figure][point].as_double_tuple_or_die();
                points.push_back(Vector3D::point(pts[0],pts[1],pts[2]));
                point = "point";
            }
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
        }
        else if (type == "BuckyBall" or type == "FractalBuckyBall")
        {
            figure1 = createBuckyBall();
        }
        else if (type == "MengerSponge")
        {
            figure1 = createMengerSponge(2);
        }
        else if (type == "Cube" or type == "FractalCube")
        {
            figure1 = createCube();
        }
        else if (type == "Tetrahedron" or type == "FractalTetrahedron")
        {
            figure1 = createTetrahedron();
        }
        else if (type == "Octahedron" or type == "FractalOctahedron")
        {
            figure1 = createOctahedron();
        }
        else if (type == "Icosahedron" or type == "FractalIcosahedron")
        {
            figure1 = createIcosahedron();
        }
        else if (type == "Dodecahedron" or type == "FractalDodecahedron")
        {
            figure1 = createDodecahedron();
        }
        else if (type == "Sphere" or type == "FractalSphere")
        {
            int n = configuration[figure]["n"].as_int_or_die();
            figure1 = createSphere(1.5,n);
        }
        else if (type == "Cone" or type == "FractalCone")
        {
            int n = configuration[figure]["n"].as_int_or_die();
            double h = configuration[figure]["height"].as_double_or_die();
            figure1 = createCone(n,h);
        }
        else if (type == "Cylinder" or type == "FractalCylinder")
        {
            int n = configuration[figure]["n"].as_int_or_die();
            double h = configuration[figure]["height"].as_double_or_die();
            figure1 = createCylinder(n,h);
        }
        else if (type == "Torus" or type == "FractalTorus")
        {
            double r = configuration[figure]["r"].as_double_or_die();
            double R = configuration[figure]["R"].as_double_or_die();
            int n = configuration[figure]["n"].as_int_or_die();
            int m = configuration[figure]["m"].as_int_or_die();
            figure1 = createTorus(r,R,n,m);
        }
        else if (type == "3DLSystem")
        {
            string input_file = configuration[figure]["inputfile"].as_string_or_die();

            LSystem3D lSystem3D;

            ifstream input_stream(input_file);
            input_stream >> lSystem3D;
            input_stream.close();

            figure1 = drawLSystem3D(lSystem3D);
        }
        else
        {
            figure1 = createCube();
        }

        figure1.color.red = colors[0];
        figure1.color.green = colors[1];
        figure1.color.blue = colors[2];

        // OMZETTEN VAN 3D NAAR 2D

        vector<double> center = configuration[figure]["center"].as_double_tuple_or_die();
        double scale = configuration[figure]["scale"].as_double_or_die();
        double rotX = configuration[figure]["rotateX"].as_double_or_die();
        double rotY = configuration[figure]["rotateY"].as_double_or_die();
        double rotZ = configuration[figure]["rotateZ"].as_double_or_die();

        double M_PI = 3.14159265358979323846;

        // 1) MATRIX OPSTELLEN

        Matrix S,M,T;

        S = scaleFigure(scale);
        M = rotateX(rotX*(M_PI/180)) * rotateY(rotY*(M_PI/180)) * rotateZ(rotZ*(M_PI/180));
        T = translate(Vector3D::point(center[0],center[1],center[2]));

        Matrix omzetMatrix = S * M * T;
        applyTransformation(figure1,omzetMatrix);

        // CREATE FRACTAL
        if (type.find("Fractal") != std::string::npos)
        {
            int nr_iterations = configuration[figure]["nrIterations"].as_int_or_die() ;
            double scale = configuration[figure]["fractalScale"].as_double_or_die();
            generateFractal(figure1, figures3D, nr_iterations, scale);
            figure = "Figure";
            continue;
        }

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

void Functies::generateFractal(Figure &figure, Figures3D &fractal, const int nr_iterations, const double scale)
{
    Figures3D new_figures = {figure};
    for(int i = 0; i < nr_iterations; i++)
    {
        Figures3D tempfigures;
        for (Figure &fig : new_figures)
        {
            for (int j = 0; j < fig.points.size(); ++j)
            {
                Figure new_figure = fig;
                applyTransformation(new_figure,scaleFigure(1/(scale)));
                Vector3D pointMixed = Vector3D::point(fig.points[j].x - new_figure.points[j].x,fig.points[j].y - new_figure.points[j].y,fig.points[j].z - new_figure.points[j].z);
                applyTransformation(new_figure,translate(pointMixed));
                tempfigures.push_back(new_figure);
            }
        }
        new_figures = tempfigures;
    }
    for (Figure figure1 : new_figures)
    {
        fractal.push_back(figure1);
    }
}

Figure Functies::createCube()
{
    Figure kubus;
    vector<vector<int>> pnt_collections =
    {
        {0,4,2,6},
        {4,1,7,2},
        {1,5,3,7},
        {5,0,6,3},
        {6,2,7,3},
        {0,5,1,4}
    };
    for (int i = 0; i < pnt_collections.size(); i++)
    {
        Face face;
        face.point_indexes = pnt_collections[i];
        kubus.faces.push_back(face);
    }
    Vector3D vector3D;
    kubus.points =
    {
        vector3D.point(1,-1,-1),
        vector3D.point(-1,1,-1),
        vector3D.point(1,1,1),
        vector3D.point(-1,-1,1),
        vector3D.point(1,1,-1),
        vector3D.point(-1,-1,-1),
        vector3D.point(1,-1,1),
        vector3D.point(-1,1,1)
    };
    return kubus;
}

Figure Functies::createTetrahedron()
{
    Figure tetrahedron;
    vector<vector<int>> pnt_collections =
    {
        {0,1,2},
        {1,3,2},
        {0,3,1},
        {0,2,3}
    };

    for (int i = 0; i < pnt_collections.size(); i++)
    {
        Face face;
        face.point_indexes = pnt_collections[i];
        tetrahedron.faces.push_back(face);
    }
    Vector3D vector3D;
    tetrahedron.points =
    {
        vector3D.point(1,-1,-1),
        vector3D.point(-1,1,-1),
        vector3D.point(1,1,1),
        vector3D.point(-1,-1,1)
    };
    return tetrahedron;
}

Figure Functies::createOctahedron()
{
    Figure octahedron;
    vector<vector<int>> pnt_collections =
    {
        {0,1,5},
        {1,2,5},
        {2,3,5},
        {3,0,5},
        {1,0,4},
        {2,1,4},
        {3,2,4},
        {0,3,4}
    };

    for (int i = 0; i < pnt_collections.size(); i++)
    {
        Face face;
        face.point_indexes = pnt_collections[i];
        octahedron.faces.push_back(face);
    }

    Vector3D vector3D;
    octahedron.points =
    {
        vector3D.point(1,0,0),
        vector3D.point(0,1,0),
        vector3D.point(-1,0,0),
        vector3D.point(0,-1,0),
        vector3D.point(0,0,-1),
        vector3D.point(0,0,1)
    };
    return octahedron;
}

Figure Functies::createIcosahedron()
{
    Figure icosahedron;
    vector<vector<int>> pnt_collections =
    {
        {0,1,2},
        {0,2,3},
        {0,3,4},
        {0,4,5},
        {0,5,1},
        {1,6,2},
        {2,6,7},
        {2,7,3},
        {3,7,8},
        {3,8,4},
        {4,8,9},
        {4,9,5},
        {5,9,10},
        {5,10,1},
        {1,10,6},
        {11,7,6},
        {11,8,7},
        {11,9,8},
        {11,10,9},
        {11,6,10}
    };

    for (int i = 0; i < pnt_collections.size(); i++)
    {
        Face face;
        face.point_indexes = pnt_collections[i];
        icosahedron.faces.push_back(face);
    }

    double M_PI = 3.14159265358979323846;
    Vector3D vector3D;
    icosahedron.points =
    {
        vector3D.point(0,0, sqrt(5)/2),
        vector3D.point(cos((2-2)*2*M_PI/5), sin((2-2)*2*M_PI/5), 0.5),
        vector3D.point(cos((3-2)*2*M_PI/5), sin((3-2)*2*M_PI/5), 0.5),
        vector3D.point(cos((4-2)*2*M_PI/5), sin((4-2)*2*M_PI/5), 0.5),
        vector3D.point(cos((5-2)*2*M_PI/5), sin((5-2)*2*M_PI/5), 0.5),
        vector3D.point(cos((6-2)*2*M_PI/5), sin((6-2)*2*M_PI/5), 0.5),
        vector3D.point(cos(M_PI/5 + (7-7)*2*M_PI/5),sin(M_PI/5 + (7-7)*2*M_PI/5),-0.5),
        vector3D.point(cos(M_PI/5 + (8-7)*2*M_PI/5),sin(M_PI/5 + (8-7)*2*M_PI/5),-0.5),
        vector3D.point(cos(M_PI/5 + (9-7)*2*M_PI/5),sin(M_PI/5 + (9-7)*2*M_PI/5),-0.5),
        vector3D.point(cos(M_PI/5 + (10-7)*2*M_PI/5),sin(M_PI/5 + (10-7)*2*M_PI/5),-0.5),
        vector3D.point(cos(M_PI/5 + (11-7)*2*M_PI/5),sin(M_PI/5 + (11-7)*2*M_PI/5),-0.5),
        vector3D.point(0,0,-sqrt(5)/2)
    };
    return icosahedron;
}

Figure Functies::createBuckyBall()
{
    Figure buckyBall;
    Figure icosahedron = createIcosahedron();
    int originalNumFaces = icosahedron.faces.size();
    int originalNumPoints = buckyBall.points.size();
    map<int,vector<Vector3D>> points;
    for (Face &face : icosahedron.faces)
    {
        int p1 = face.point_indexes[0];
        int p2 = face.point_indexes[1];
        int p3 = face.point_indexes[2];
        // ADD POINTS
        buckyBall.points.push_back(getPoint(icosahedron.points[p1], icosahedron.points[p2], 1.0/3.0));
        buckyBall.points.push_back(getPoint(icosahedron.points[p2], icosahedron.points[p3], 1.0/3.0));
        buckyBall.points.push_back(getPoint(icosahedron.points[p3], icosahedron.points[p1], 1.0/3.0));
        buckyBall.points.push_back(getPoint(icosahedron.points[p1], icosahedron.points[p2], 2.0/3.0));
        buckyBall.points.push_back(getPoint(icosahedron.points[p2], icosahedron.points[p3], 2.0/3.0));
        buckyBall.points.push_back(getPoint(icosahedron.points[p3], icosahedron.points[p1], 2.0/3.0));
        // INDEXES
        int p4 = originalNumPoints;
        int p5 = originalNumPoints + 1;
        int p6 = originalNumPoints + 2;
        int p7 = originalNumPoints + 3;
        int p8 = originalNumPoints + 4;
        int p9 = originalNumPoints + 5;
        // ADD FACE
        Face face1;
        face1.point_indexes = {p4, p9, p6, p8, p5, p7};
        buckyBall.faces.push_back(face1);
        originalNumPoints += 6;
        if(points.count(face.point_indexes[0]) < 1)
        {
            points[face.point_indexes[0]] = {getPoint(icosahedron.points[p1], icosahedron.points[p2], 1.0 / 3.0), getPoint(icosahedron.points[p3], icosahedron.points[p1], 2.0 / 3.0)};
        }
        else
        {
            points[face.point_indexes[0]].push_back(getPoint(icosahedron.points[p1], icosahedron.points[p2], 1.0 / 3.0));
            points[face.point_indexes[0]].push_back(getPoint(icosahedron.points[p3], icosahedron.points[p1], 2.0 / 3.0));
        }
        if(points.count(face.point_indexes[1]) < 1)
        {
            points[face.point_indexes[1]] = {getPoint(icosahedron.points[p1], icosahedron.points[p2], 2.0 / 3.0), getPoint(icosahedron.points[p2], icosahedron.points[p3], 1.0 / 3.0)};
        }
        else
        {
            points[face.point_indexes[1]].push_back(getPoint(icosahedron.points[p1], icosahedron.points[p2], 2.0 / 3.0));
            points[face.point_indexes[1]].push_back(getPoint(icosahedron.points[p2], icosahedron.points[p3], 1.0 / 3.0));
        }
        if(points.count(face.point_indexes[2]) < 1)
        {
            points[face.point_indexes[2]] = {getPoint(icosahedron.points[p2], icosahedron.points[p3], 2.0 / 3.0), getPoint(icosahedron.points[p3], icosahedron.points[p1], 1.0 / 3.0)};
        }
        else
        {
            points[face.point_indexes[2]].push_back(getPoint(icosahedron.points[p2], icosahedron.points[p3], 2.0 / 3.0));
            points[face.point_indexes[2]].push_back(getPoint(icosahedron.points[p3], icosahedron.points[p1], 1.0 / 3.0));
        }
    }
    for(pair<const int, vector<Vector3D>>& i: points)
    {
        vector<Vector3D> vijfPunten;
        for(int j = 0; j < i.second.size(); j++)
        {
            for(int k = j+1; k < i.second.size(); k++)
            {
                double lengthx = abs(i.second[j].x - i.second[k].x);
                double lengthy =  abs(i.second[j].y - i.second[k].y);
                double lengthz =  abs(i.second[j].z - i.second[k].z);
                if(lengthx < 0.00001 and lengthy < 0.00001 and lengthz < 0.00001)
                {
                    vijfPunten.push_back(i.second[j]);
                    break;
                }
            }
        }
        i.second = vijfPunten;
    }
    int length = buckyBall.points.size();
    for(pair<const int, vector<Vector3D>>& i: points)
    {
        vector<int> volgorde = {0};
        for(int j = 0; j < 4; j++)
        {
            pair<int,double> kort;
            for(int k = volgorde[j]; k < 5; k++)
            {
                if(k != 4)
                {
                    if(count(volgorde.begin(), volgorde.end(), k + 1) == 0)
                    {
                        kort.first = k + 1;
                        kort.second = sqrt(pow(i.second[volgorde[j]].x - i.second[k + 1].x, 2) + pow(i.second[volgorde[j]].y - i.second[k + 1].y, 2) + pow(i.second[volgorde[j]].z - i.second[k + 1].z, 2));
                        break;
                    }
                }
                else
                {
                    for(int f = 0; f < 4; f++)
                    {
                        if(count(volgorde.begin(), volgorde.end(), f) == 0)
                        {
                            kort.first = f;
                            kort.second = sqrt(pow(i.second[volgorde[j]].x - i.second[f].x, 2) + pow(i.second[volgorde[j]].y - i.second[f].y, 2) + pow(i.second[volgorde[j]].z - i.second[f].z, 2));

                        }
                    }
                }
            }
            for(int k = 0; k < i.second.size(); k++)
            {
                if(count(volgorde.begin(), volgorde.end(), k) == 0)
                {
                    double distance = sqrt(pow(i.second[volgorde[j]].x - i.second[k].x, 2) + pow(i.second[volgorde[j]].y - i.second[k].y, 2) + pow(i.second[volgorde[j]].z - i.second[k].z, 2));
                    if(distance < kort.second)
                    {
                        kort.first = k;
                        kort.second = distance;
                    }
                }
            }
            volgorde.push_back(kort.first);
        }
        Face face;
        face.point_indexes = {length, length + 1, length + 2, length + 3, length + 4};
        buckyBall.faces.push_back(face);
        for(int j : volgorde)
        {
            buckyBall.points.push_back(i.second[j]);
        }
        length += 5;
    }
    return buckyBall;
}

Figure Functies::createMengerSponge(int nrIterations)
{
    Figure mengerSponge = createCube();
    int originalNumPoints = mengerSponge.faces.size();
    for (int x = 0; x < 6; x++)
    {
        int p1 = mengerSponge.faces[x].point_indexes[0];
        int p2 = mengerSponge.faces[x].point_indexes[1];
        int p3 = mengerSponge.faces[x].point_indexes[2];
        int p4 = mengerSponge.faces[x].point_indexes[3];

        // ADD POINTS
        Vector3D AB1 = getPoint(mengerSponge.points[p1], mengerSponge.points[p2], 1.0/3.0);
        Vector3D AB2 = getPoint(mengerSponge.points[p1], mengerSponge.points[p2], 2.0/3.0);
        Vector3D BC1 = getPoint(mengerSponge.points[p2], mengerSponge.points[p3], 1.0/3.0);
        Vector3D BC2 = getPoint(mengerSponge.points[p2], mengerSponge.points[p3], 2.0/3.0);
        Vector3D CD1 = getPoint(mengerSponge.points[p3], mengerSponge.points[p4], 1.0/3.0);
        Vector3D CD2 = getPoint(mengerSponge.points[p3], mengerSponge.points[p4], 2.0/3.0);
        Vector3D DA1 = getPoint(mengerSponge.points[p4], mengerSponge.points[p1], 1.0/3.0);
        Vector3D DA2 = getPoint(mengerSponge.points[p4], mengerSponge.points[p1], 2.0/3.0);

        mengerSponge.points.push_back(AB1);
        mengerSponge.points.push_back(AB2);
        mengerSponge.points.push_back(BC1);
        mengerSponge.points.push_back(BC2);
        mengerSponge.points.push_back(CD1);
        mengerSponge.points.push_back(CD2);
        mengerSponge.points.push_back(DA1);
        mengerSponge.points.push_back(DA2);

        // INDEXES
        int ab1 = originalNumPoints;
        int ab2 = originalNumPoints + 1;
        int bc1 = originalNumPoints + 2;
        int bc2 = originalNumPoints + 3;
        int cd1 = originalNumPoints + 4;
        int cd2 = originalNumPoints + 5;
        int da1 = originalNumPoints + 6;
        int da2 = originalNumPoints + 7;
        // ADD FACE
        vector<vector<int>> points
        {
            {p1,p2,bc1,da2},
            {da2,bc1,bc2,da1},
            {p4,p3,bc2,da1}
        };
        for (vector<int> i : points)
        {
            Face vlak;
            vlak.point_indexes = i;
            mengerSponge.faces.push_back(vlak);
        }
        originalNumPoints += 8;
    }
    return mengerSponge;
}


Vector3D Functies::getPoint(Vector3D p1, Vector3D p2, double factor)
{
    double x = p1.x + factor * (p2.x - p1.x);
    double y = p1.y + factor * (p2.y - p1.y);
    double z = p1.z + factor * (p2.z - p1.z);

    Vector3D vector3D;
    return vector3D.point(x,y,z);
}

Figure Functies::createDodecahedron()
{
    Figure dodecahedron;
    Figure isocohedron = createIcosahedron();

    vector<vector<int>> pnt_collections =
    {
        {1-1,2-1,3-1,4-1,5-1},
        {1-1,6-1,7-1,8-1,2-1},
        {2-1,8-1,9-1,10-1,3-1},
        {3-1,10-1,11-1,12-1,4-1},
        {4-1,12-1,13-1,14-1,5-1},
        {5-1,14-1,15-1,6-1,1-1},
        {20-1,19-1,18-1,17-1,16-1},
        {20-1,15-1,14-1,13-1,19-1},
        {19-1,13-1,12-1,11-1,18-1},
        {18-1,11-1,10-1,9-1,17-1},
        {17-1,9-1,8-1,7-1,16-1},
        {16-1,7-1,6-1,15-1,20-1}
    };

    for (int i = 0; i < pnt_collections.size(); i++)
    {
        Face face;
        face.point_indexes = pnt_collections[i];
        dodecahedron.faces.push_back(face);
    }

    for (Face punt : isocohedron.faces)
    {
        Vector3D vector3D;
        vector3D.x = (isocohedron.points[punt.point_indexes[0]].x + isocohedron.points[punt.point_indexes[1]].x + isocohedron.points[punt.point_indexes[2]].x) / 3;
        vector3D.y = (isocohedron.points[punt.point_indexes[0]].y + isocohedron.points[punt.point_indexes[1]].y + isocohedron.points[punt.point_indexes[2]].y) / 3;
        vector3D.z = (isocohedron.points[punt.point_indexes[0]].z + isocohedron.points[punt.point_indexes[1]].z + isocohedron.points[punt.point_indexes[2]].z) / 3;

        dodecahedron.points.push_back(vector3D);
    }

    return dodecahedron;
}

Figure Functies::createSphere(const double radius, const int n)
{
    Figure sphere = createIcosahedron();

    vector<Face> faces = sphere.faces;
    vector<Vector3D> points = sphere.points;
    for (int i = 0; i < n; i++)
    {
        sphere.faces.clear();
        sphere.points.clear();

        for (Face &face : faces)
        {
            Vector3D A = points[face.point_indexes[0]];
            Vector3D B = points[face.point_indexes[1]];
            Vector3D C = points[face.point_indexes[2]];

            Vector3D D = Vector3D::point((A.x+B.x)/2,(A.y+B.y)/2,(A.z+B.z)/2);
            Vector3D E = Vector3D::point((A.x+C.x)/2,(A.y+C.y)/2,(A.z+C.z)/2);
            Vector3D F = Vector3D::point((B.x+C.x)/2,(B.y+C.y)/2,(B.z+C.z)/2);

            sphere.points.push_back(A);
            sphere.points.push_back(B);
            sphere.points.push_back(C);
            sphere.points.push_back(D);
            sphere.points.push_back(E);
            sphere.points.push_back(F);

            int a = sphere.points.size()-6;
            int b = sphere.points.size()-5;
            int c = sphere.points.size()-4;
            int d = sphere.points.size()-3;
            int e = sphere.points.size()-2;
            int f = sphere.points.size()-1;

            vector<vector<int>> pts_collections =
            {
                {a,d,e},
                {b,f,d},
                {c,e,f},
                {d,f,e}
            };

            for (vector<int> &vlak : pts_collections)
            {
                Face face;
                face.point_indexes = vlak;
                sphere.faces.push_back(face);
            }
        }
        faces = sphere.faces;
        points = sphere.points;
    }

    for (Vector3D &pnt : sphere.points)
    {
        pnt.normalise();
    }

    return sphere;
}

Figure Functies::createCone(const int n, const double h)
{
    Figure cone;
    double PI = 3.14159265358979323846;
    //                                                    PUNTEN
    for (int i = 0; i < n + 1; i++)
    {
        if (i == n)
        {
            // TOP KEGEL
            Vector3D top = Vector3D::point(0,0,h);
            cone.points.push_back(top);
        }
        else
        {
            // GRONDPUNT
            Vector3D grondpunt = Vector3D::point(cos((2*i*PI)/n),sin((2*i*PI)/n),0);
            cone.points.push_back(grondpunt);
        }
    }
    //                                                    VLAKKEN
    vector<vector<int>> vlakken;
    for (int j = 0; j < n + 1; j++)
    {
        if (j == n)
        {
            vector<int> grondvlak;
            for (int x = 0; x <= n; x++)
            {
                grondvlak.push_back(n-x);
            }
            vlakken.push_back(grondvlak);
        }
        else
        {
            vector<int> mantelvlak{j,(j+1)%n,n};
            vlakken.push_back(mantelvlak);
        }
    }
    for (int i = 0; i < vlakken.size(); i++)
    {
        Face face;
        face.point_indexes = vlakken[i];
        cone.faces.push_back(face);
    }
    return cone;
}

Figure Functies::createCylinder(const int n, const double h)
{
    Figure cylinder;
    double PI = 3.14159265358979323846;

    //                                                    PUNTEN

    cylinder.points.push_back(Vector3D::point(0, 0, h)); // top point
    cylinder.points.push_back(Vector3D::point(0, 0, 0)); // bottom point

    for (int i = 0; i < 2*n; i++)
    {
        // TOP KEGEL
        Vector3D toppunt = Vector3D::point(cos((2*i*PI)/n),sin((2*i*PI)/n),h);
        cylinder.points.push_back(toppunt);
        // GRONDPUNT
        Vector3D grondpunt = Vector3D::point(cos((2*i*PI)/n),sin((2*i*PI)/n),0);
        cylinder.points.push_back(grondpunt);
    }
    //                                                    VLAKKEN
    vector<vector<int>> vlakken;
    for (int j = 0; j < 2*n; j+=2)
    {
        vector<int> mantelvlak{j,(j+1)%(2*n),(j+3)%(2*n),(j+2)%(2*n)};
        vlakken.push_back(mantelvlak);
    }

    //                                                      TOP
    for (int i = 0; i < n; i++)
    {
        vector<int> topvlak{i*2, (i*2+2)%(2*n), 2*n};
        vlakken.push_back(topvlak);
    }

    //                                                     ONDER
    for (int i = 0; i < n; i++)
    {
        vector<int> bottomvlak{i*2+1, (i*2+3)%(2*n), 2*n+1};
        vlakken.push_back(bottomvlak);
    }

    for (int i = 0; i < vlakken.size(); i++)
    {
        Face face;
        face.point_indexes = vlakken[i];
        cylinder.faces.push_back(face);
    }

    return cylinder;
}

Figure Functies::createTorus(const double r, const double R, const int n, const int m)
{
    Figure torus;
    const double PI = 3.14159265358979323846;
    //                                                    PUNTEN
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            double u = (2*i*PI)/n;
            double v = (2*j*PI)/m;
            torus.points.push_back(Vector3D::point((R + r * cos(v))*cos(u),(R + r * cos(v))*sin(u),r*sin(v)));
        }
    }
    //                                                    VLAKKEN
    vector<vector<int>> vlakken;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            vlakken.push_back({i*m + j, ((i+1)%n)*m + j, ((i+1)%n)*m + (j+1)%m, i*m + (j+1)%m});
        }
    }
    for (int i = 0; i < vlakken.size(); i++)
    {
        Face face;
        face.point_indexes = vlakken[i];
        torus.faces.push_back(face);
    }
    return torus;
}



Figure Functies::drawLSystem3D(const LSystem3D &lSystem3D)
{
    Figure figure;
    const double PI = 3.14159265358979323846;
    // DE COMPONENTEN
    set<char> alfabet = lSystem3D.get_alphabet();
    string initiator = lSystem3D.get_initiator();
    unsigned int iterations = lSystem3D.get_nr_iterations();
    double angle = lSystem3D.get_angle();
    //REPLACEMENTS
    vector<pair<char,string>> replacements;
    for (char i : alfabet)
    {
        string str = lSystem3D.get_replacement(i);
        pair<char,string> letter(i,str);
        replacements.push_back(letter);
    }
    angle *= (PI/180);

    tekenReplace(initiator,replacements,iterations);
    leesString(angle,initiator,figure,lSystem3D);

    return figure;
}

void Functies::leesString(const double angle, string string, Figure& figure, const LSystem3D &lSystem3D)
{
    stack<vector<Vector3D>> punten_stack;

    Vector3D H = Vector3D::point(1,0,0);
    Vector3D L = Vector3D::point(0,1,0);
    Vector3D U = Vector3D::point(0,0,1);
    Vector3D A = Vector3D::point(0,0,0);
    Vector3D old;

    figure.points.push_back(A);
    int count = 1;
    for (const char k : string)
    {
        if (k == '+')
        {
            old = H;
            H = H * cos(angle) + L * sin(angle);
            L = -old * sin(angle) + L * cos(angle);
        }
        else if (k == '-')
        {
            old = H;
            H = H * cos(-angle) + L * sin(-angle);
            L = -old * sin(-angle) + L * cos(-angle);
        }
        else if (k == '^')
        {
            old = H;
            H = H * cos(angle) + U * sin(angle);
            U = -old * sin(angle) + U * cos(angle);
        }
        else if (k == '&')
        {
            old = H;
            H = H * cos(-angle) + U * sin(-angle);
            U = -old * sin(-angle) + U * cos(-angle);
        }
        else if (k == '/')
        {
            old = L;
            L = L * cos(-angle) - U * sin (-angle);
            U = old * sin(-angle) + U * cos(-angle);
        }
        else if (k == '\\')
        {
            old = L;
            L = L * cos(angle) - U * sin (angle);
            U = old * sin(angle) + U * cos(angle);
        }
        else if (k == '|')
        {
            H = -H;
            L = -L;
        }
        else if (k == '(')
        {
            punten_stack.push({H,L,U,A});
        }
        else if (k == ')')
        {
            vector<Vector3D> top = punten_stack.top();

            H = top[0];
            L = top[1];
            U = top[2];
            A = top[3];

            figure.points.push_back(A);
            punten_stack.pop();
            count++;
        }
        else
        {
            A += H;
            figure.points.push_back(A);
            if (lSystem3D.draw(k))
            {
                Face face;
                face.point_indexes = {count - 1, count};
                figure.faces.push_back(face);
            }
            count++;
        }
    }
}