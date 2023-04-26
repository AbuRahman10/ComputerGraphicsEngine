//
// Created by Abobaker Rahman on 12/04/2023.
//

#include <sstream>
#include "ZBuffer.h"
#include "limits"
#include "math.h"
#include <assert.h>

inline int roundToInt(double d) { return static_cast<int>(round(d)); }

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
