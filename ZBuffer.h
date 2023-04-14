//
// Created by Abobaker Rahman on 12/04/2023.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H

#include "iostream"
#include "vector"

using namespace std;

class ZBuffer: public vector<vector<double>>
{
public:

    const int width;
    const int height;

    ZBuffer(const int width, const int height);
};


#endif //ENGINE_ZBUFFER_H
