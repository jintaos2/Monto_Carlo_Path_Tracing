#ifndef TOOLS_H
#define TOOLS_H

#include <glm.hpp>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

int toint(std::string input)
{
    int ret;
    std::istringstream ss(input);
    ss >> ret;
    return ret;
}
template <typename T>
std::string tostring(T in)
{
    std::stringstream ss;
    ss << in;
    std::string ret;
    ss >> ret;
    return ret;
}
float tofloat(std::string input)
{
    float ret;
    std::istringstream ss(input);
    ss >> ret;
    return ret;
}
std::string print_vec3(glm::dvec3 in)
{
    std::stringstream ss;
    ss << '(' << in.x << ',' << in.y << ',' << in.z << ')';
    std::string ret;
    ss >> ret;
    return ret;
}
std::string print_mat3x4(glm::imat3x4 in)
{
    std::stringstream ss;
    std::string ret;
    ss << '(';
    for (int x = 0; x < 3; ++x)
    {
        ss << in[x][0] << ',' << in[x][1] << ',' << in[x][2] << ',' << in[x][3] << ';';
    }
    ss << ')';
    ss >> ret;
    return ret;
}
std::string print_mat3x3(glm::dmat3x3 in)
{
    std::stringstream ss;
    std::string ret;
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
inline const T &max3(const T &a, const T &b, const T &c)
{
    return max_(max_(a, b), c);
}
template <typename T>
inline const T &min3(const T &a, const T &b, const T &c)
{
    return min_(min_(a, b), c);
}
template <typename T>
inline const T &between(const T &xmin, const T &xmax, const T &x)
{
    return min_(max_(xmin, x), xmax);
}
// TODO hdr to ldr
inline uint32_t toRGB(glm::vec3 a)
{
    int R = between(0, 255, int(a.x * 255.0));
    int G = between(0, 255, int(a.y * 255.0));
    int B = between(0, 255, int(a.z * 255.0));
    return 0xff000000 | B << 16 | G << 8 | R;
}

#endif
