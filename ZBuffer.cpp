//
// Created by Abobaker Rahman on 12/04/2023.
//

#include <sstream>
#include "ZBuffer.h"
#include "limits"
#include "math.h"
#include "Point2D.h"
#include "Line2D.h"

ZBuffer::ZBuffer(const int width, const int height) : width(width), height(height)
{
    for (int i = 0; i < height; i++)
    {
        vector<double> rij;
        for (int j = 0; j < width; j++)
        {
            double infinity = numeric_limits<double>::infinity();
            rij.push_back(infinity);
        }
        this->push_back(rij);
    }
}

void ZBuffer::draw_zbuf_line(EasyImage &image, unsigned int x0, unsigned int y0, double z0, unsigned int x1, unsigned int y1, double z1,const Color &color)
{
    if (x0 >= image.get_width() || y0 >= image.get_height() || x1 >= image.get_width() || y1 > image.get_height())
    {
        std::stringstream ss;
        ss << "Drawing line from (" << x0 << "," << y0 << ") to (" << x1 << "," << y1 << ") in image of width "
           << this->width << " and height " << this->height;
        throw std::runtime_error(ss.str());
    }
    if (x0 == x1 and y0 == y1)
    {
        double max_z = max(z0,z1);
        double inv_z = 1/max_z;
        if (inv_z < (*this)[y0][x0])
        {
            (*this)[y0][x0] = inv_z;
            (image)(x0, y0) = color;
        }
    }
    else if (x0 == x1)
    {
        if (y0 > y1)
        {
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1);
        }
        //special case for x0 == x1
        double a = max(y0, y1) - min(y0, y1);
        int n = 0;
        for (unsigned int i = min(y0, y1); i <= max(y0, y1); i++)
        {
            double p = 1-(n/a);
            double inv_z = (p/z0) + ((1-p)/z1);
            if (inv_z < (*this)[i][x0])
            {
                (*this)[i][x0] = inv_z;
                (image)(x0, i) = color;
            }
            n++;
        }
    }
    else if (y0 == y1)
    {
        if (x0 > x1)
        {
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1);
        }
        //special case for y0 == y1
        double a = max(x0, x1) - min(x0, x1);
        int n = 0;
        for (unsigned int i = min(x0, x1); i <= max(x0, x1); i++)
        {
            double p = 1-(n/a);
            double inv_z = (p/z0) + ((1-p)/z1);
            if (inv_z < (*this)[y0][i])
            {
                (*this)[y0][i] = inv_z;
                (image)(i, y0) = color;
            }
            n++;
        }
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            double a = x1 - x0;
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                double p = 1-(i/a);
                double inv_z = (p/z0) + ((1-p)/z1);
                if (inv_z < (*this)[(unsigned int) round(y0 + m * i)][x0 + i])
                {
                    (*this)[(unsigned int) round(y0 + m * i)][x0 + i] = inv_z;
                    (image)(x0 + i, (unsigned int) round(y0 + m * i)) = color;
                }
            }
        }
        else if (m > 1.0)
        {
            double a = y1 - y0;
            for (unsigned int i = 0; i <= (y1 - y0); i++)
            {
                double p = 1-(i/a);
                double inv_z = (p/z0) + ((1-p)/z1);
                if (inv_z < (*this)[y0 + i][(unsigned int) round(x0 + (i / m))])
                {
                    (*this)[y0 + i][(unsigned int) round(x0 + (i / m))] = inv_z;
                    (image)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
                }
            }
        }
        else if (m < -1.0)
        {
            double a = y0-y1;
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                double p = 1-(i/a);
                double inv_z = (p/z0) + ((1-p)/z1);
                if (inv_z < (*this)[y0 - i][(unsigned int) round(x0 - (i / m))])
                {
                    (*this)[y0 - i][(unsigned int) round(x0 - (i / m))] = inv_z;
                    (image)((unsigned int) round(x0 - (i / m)), y0 - i) = color;
                }
            }
        }
    }
}

vector<Face> ZBuffer::triangulate(const Face &face)
{
    vector<Face> faces;
    for (int i = 1; i < face.point_indexes.size() - 1; i++)
    {
        Face face1;
        face1.point_indexes = {face.point_indexes[0],face.point_indexes[i],face.point_indexes[i+1]};
        faces.push_back(face1);
    }
    return faces;
}

void ZBuffer::draw_zbuf_triag(EasyImage &image, const Vector3D &A, const Vector3D &B, const Vector3D &C, double d, double dx, double dy, Color color)
{
    Point2D A_(((d*A.x)/(-A.z))+dx,((d*A.y)/(-A.z))+dy);
    Point2D B_(((d*B.x)/(-B.z))+dx,((d*B.y)/(-B.z))+dy);
    Point2D C_(((d*C.x)/(-C.z))+dx,((d*C.y)/(-C.z))+dy);

    int yMin = lround(min(min(A_.y,B_.y),C_.y) + 0.5);
    int yMax = lround(max(max(A_.y,B_.y),C_.y) - 0.5);

    double xG = (A_.x + B_.x + C_.x)/3;
    double yG = (A_.y + B_.y + C_.y)/3;
    double inv_zG = (1/(3*A.z)) + (1/(3*B.z)) + (1/(3*C.z));

    Vector3D u = B - A;
    Vector3D v = C - A;
    double w1 = (u.y * v.z) - (u.z * v.y);
    double w2 = (u.z * v.x) - (u.x * v.z);
    double w3 = (u.x * v.y) - (u.y * v.x);
    double k = w1 * A.x + w2 * A.y + w3 * A.z;

    double dzdx = w1/(-(d*k));
    double dzdy = w2/(-(d*k));

    for (int y = yMin; y <= yMax; y++)
    {
        double xL_AB = numeric_limits<double>::infinity();
        double xL_AC = numeric_limits<double>::infinity();
        double xL_BC = numeric_limits<double>::infinity();

        double xR_AB = -numeric_limits<double>::infinity();
        double xR_AC = -numeric_limits<double>::infinity();
        double xR_BC = -numeric_limits<double>::infinity();

        // HANDMATIGE FOR LOOP DOOR DE DRIE LIJNEN
        // ---------------------------------------
        // LIJNSTUK AB
        if (((y-A.y)*(y-B.y) <= 0 and A.y != B.y))
        {
            double xI = B.x + (((A.x - B.x)*(y-B.y))/(A.y-B.y));
            xL_AB = xR_AB = xI;
        }
        // LIJNSTUK AC
        else if (((y-A.y)*(y-C.y) <= 0 and A.y != C.y))
        {
            double xI = C.x + (((A.x - C.x)*(y-C.y))/(A.y-C.y));
            xL_AC = xR_AC = xI;
        }
        // LIJNSTUK BC
        else if (((y-B.y)*(y-C.y) <= 0 and B.y != C.y))
        {
            double xI = C.x + (((B.x - C.x)*(y-C.y))/(B.y-C.y));
            xL_BC = xR_BC = xI;
        }
        double xL = lround(min(min(xL_AB,xL_AC),xL_BC) + 0.5);
        double xR = lround(max(max(xR_AB,xR_AC),xR_BC) - 0.5);

        for (int x = xL; x <= xR; x++)
        {
            double inv_z = 1.0001 * inv_zG + (x - xG) * dzdx + (y - yG) * dzdy;
            if (inv_z < (*this)[y][x])
            {
                (*this)[y][x] = inv_z;
                (image)(x,y) = color;
            }
        }
    }
}
