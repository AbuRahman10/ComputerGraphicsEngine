//
// Created by Abobaker Rahman on 12/04/2023.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H

#include "iostream"
#include "vector"
#include "easy_image.h"
#include "Face.h"
#include "vector3d.h"

using namespace std;
using namespace img;

class ZBuffer: public vector<vector<double>>
{
public:

    const int width;
    const int height;

    ZBuffer(const int width, const int height);

    void draw_zbuf_line(EasyImage &image,unsigned int x0, unsigned int y0, double z0, unsigned int x1, unsigned int y1, double z1, const Color &color);

    vector<Face> triangulate(const Face &face);

    void draw_zbuf_triag(EasyImage &image, Vector3D const &A, Vector3D const &B, Vector3D const &C, double d, double dx, double dy, Color color);
};


#endif //ENGINE_ZBUFFER_H
