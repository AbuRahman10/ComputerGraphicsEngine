//
// Created by Abobaker Rahman on 12/04/2023.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H

#include "iostream"
#include "vector"
#include "easy_image.h"

using namespace std;
using namespace img;

class ZBuffer: public vector<vector<double>>
{
public:

    const int width;
    const int height;

    ZBuffer(const int width, const int height);

    void draw_zbuf_line(ZBuffer &zBuffer,EasyImage &image,const unsigned int x0, const unsigned int y0, const double z0, const unsigned int x1, const unsigned int y1, const double z1, const Color &color);
};


#endif //ENGINE_ZBUFFER_H
