#ifndef TOOLS_H
#define TOOLS_H

#include <glm.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

using namespace glm;
using namespace std;

int toint(string input)
{
    int ret;
    istringstream ss(input);
    ss >> ret;
    return ret;
}
template <typename T>
string tostring(T in)
{
    stringstream ss;
    ss << in;
    string ret;
    ss >> ret;
    return ret;
}
float tofloat(string input)
{
    float ret;
    istringstream ss(input);
    ss >> ret;
    return ret;
}
string print_vec3(vec3 in)
{
    stringstream ss;
    ss << '(' << in.x << ',' << in.y << ',' << in.z << ')';
    string ret;
    ss >> ret;
    return ret;
}
string print_mat3x4(imat3x4 in)
{
    stringstream ss;
    string ret;
    ss << '(';
    for (int x = 0; x < 3; ++x)
    {
        ss << in[x][0] << ',' << in[x][1] << ',' << in[x][2] << ',' << in[x][3] << ';';
    }
    ss << ')';
    ss >> ret;
    return ret;
}
string print_mat3x3(mat3x3 in)
{
    stringstream ss;
    string ret;
    ss << '(';
    for (int x = 0; x < 3; ++x)
    {
        ss << in[x][0] << ',' << in[x][1] << ',' << in[x][2] << ';';
    }
    ss << ')';
    ss >> ret;
    return ret;
}

template <typename T>
inline const T &max_(const T &__a, const T &__b)
{
    return __a < __b ? __b : __a;
}
template <typename T>
inline const T &min_(const T &__a, const T &__b)
{
    return __a < __b ? __a : __b;
}
template <typename T>
inline T max3(T a, T b, T c)
{
    return max_(max_(a, b), c);
}
template <typename T>
inline T min3(T a, T b, T c)
{
    return min_(min_(a, b), c);
}
template <typename T>
inline T between(T xmin, T xmax, T x)
{
    return min_(max_(xmin, x), xmax);
}
inline uint32_t toRGB(vec3 a)
{
    int R = between(0, 255, int(a.x * 255.0));
    int G = between(0, 255, int(a.y * 255.0));
    int B = between(0, 255, int(a.z * 255.0));
    return 0xff000000 | B << 16 | G << 8 | R;
}

#endif
