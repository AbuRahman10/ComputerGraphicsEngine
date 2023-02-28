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

using namespace std;
using namespace img;

class Functies
{
public:
    static EasyImage draw2DLines(const Lines2D &lines, const int size);
};

#endif //ENGINE_FUNCTIES_H
