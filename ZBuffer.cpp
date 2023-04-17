//
// Created by Abobaker Rahman on 12/04/2023.
//

#include "ZBuffer.h"

ZBuffer::ZBuffer(const int width, const int height) : width(width), height(height) {}

void ZBuffer::draw_zbuf_line(ZBuffer &zBuffer, EasyImage &image, const unsigned int x0, const unsigned int y0,const double z0, const unsigned int x1, const unsigned int y1, const double z1,const Color &color)
{
    image.draw_line(x0,y0,x1,y1,color);
}
