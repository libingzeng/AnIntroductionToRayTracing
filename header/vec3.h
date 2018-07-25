#ifndef VEC3_H
#define VEC3_H

#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include "log.h"
#include <iomanip>

using namespace std;

class vec3
{
    public:
        vec3() {}
        vec3(float e0, float e1, float e2) {e[0] = e0; e[1] = e1; e[2] = e2; }
        inline float x() const { return e[0]; }
        inline float y() const { return e[1]; }
        inline float z() const { return e[2]; }
        inline float r() const { return e[0]; }
        inline float g() const { return e[1]; }
        inline float b() const { return e[2]; }

        inline const vec3& operator+() const { return *this; }
        inline vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
        inline float operator[](int i) const { return e[i]; }
        inline float& operator[](int i){ return e[i]; };

        inline vec3& operator+=(const vec3 &v2);
        inline vec3& operator-=(const vec3 &v2);
        inline vec3& operator*=(const vec3 &v2);
        inline vec3& operator/=(const vec3 &v2);
        inline vec3& operator*=(const float t);
        inline vec3& operator/=(const float t);
        inline float length() const
        {
            return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
        }
        inline float squared_length() const
        {
            return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
        }
        inline void make_unit_vector();
        inline vec3 unit_vector(vec3 v);

        float e[3];
};

inline std::istream& operator>>(std::istream &is, vec3 &t)
{
    is >> t.e[0] >> t.e[1] >> t.e[2];
    return is;
}

inline std::ostream& operator<<(std::ostream &os, const vec3 &t)
{
    os << t.e[0] << " " << t.e[1] << " " << t.e[2];
    return os;
}

inline void vec3::make_unit_vector()
{
    float k = 1.0 / sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
    e[0] *= k;
    e[1] *= k;
    e[2] *= k;
}

inline vec3 operator+(const vec3 &v1, const vec3 &v2)
{
    return vec3(v1.e[0] + v2.e[0], v1.e[1] + v2.e[1], v1.e[2] + v2.e[2]);
}

inline vec3 operator-(const vec3 &v1, const vec3 &v2)
{
    return vec3(v1.e[0] - v2.e[0], v1.e[1] - v2.e[1], v1.e[2] - v2.e[2]);
}

inline vec3 operator*(const vec3 &v1, const vec3 &v2)
{
    return vec3(v1.e[0] * v2.e[0], v1.e[1] * v2.e[1], v1.e[2] * v2.e[2]);
}

inline vec3 operator/(const vec3 &v1, const vec3 &v2)
{
    return vec3(v1.e[0] / v2.e[0], v1.e[1] / v2.e[1], v1.e[2] / v2.e[2]);
}

inline vec3 operator*(float t, const vec3 &v)
{
    return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline vec3 operator/(vec3 v, float t)
{
    return vec3(v.e[0]/t, v.e[1]/t, v.e[2]/t);
}

inline vec3 operator*(const vec3 &v, float t)
{
    return vec3(t*v.e[0], t*v.e[1], t*v.e[2]);
}

inline float dot(const vec3 &v1, const vec3 &v2)
{
    return v1.e[0] * v2.e[0] + v1.e[1] * v2.e[1] + v1.e[2] * v2.e[2];
}

inline vec3 cross(const vec3 &v1, const vec3 &v2)
{
    return vec3( (v1.e[1]*v2.e[2] - v1.e[2]*v2.e[1]),
                 (v1.e[2]*v2.e[0] - v1.e[0]*v2.e[2]),
                 (v1.e[0]*v2.e[1] - v1.e[1]*v2.e[0]));
}

inline bool vector_equ(const vec3 &v1, const vec3 &v2)
{
    return ((v1.e[0]==v2.e[0]) && (v1.e[1]==v2.e[1]) && (v1.e[2]==v2.e[2])) ? true : false;
}

inline vec3& vec3::operator+=(const vec3 &v)
{
    e[0] += v.e[0];
    e[1] += v.e[1];
    e[2] += v.e[2];

    return *this;
}

inline vec3& vec3::operator*=(const vec3 &v)
{
    e[0] *= v.e[0];
    e[1] *= v.e[1];
    e[2] *= v.e[2];

    return *this;
}

inline vec3& vec3::operator-=(const vec3 &v)
{
    e[0] -= v.e[0];
    e[1] -= v.e[1];
    e[2] -= v.e[2];

    return *this;
}

inline vec3& vec3::operator/=(const vec3 &v)
{
    e[0] /= v.e[0];
    e[1] /= v.e[1];
    e[2] /= v.e[2];

    return *this;
}

inline vec3& vec3::operator*=(const float t)
{
    e[0] *= t;
    e[1] *= t;
    e[2] *= t;

    return *this;
}

inline vec3& vec3::operator/=(const float t)
{
    float k = 1.0 / t;
    e[0] *= k;
    e[1] *= k;
    e[2] *= k;

    return *this;
}

inline vec3 unit_vector(vec3 v)
{
    return v / v.length();
}

#define EPS 1e-10
vec3 get_vector_v(const vec3& vector_u, float angle);
vec3 vector_trans(const vec3& v1, const vec3& u, const vec3& v, const vec3& w);
vec3 vector_trans_back(const vec3& v1, const vec3& u, const vec3& v, const vec3& w);
bool get_vector_vw(vec3& u, float angle, vec3& v, vec3& w);

//bool roots_quadratic_equation2(float a, float b, float c, float (&roots)[3]);
float* roots_quadratic_equation(double a, double b, double c);
float* roots_cubic_equation(double a, double b, double c, double d);
//float* roots_quartic_equation(float a, float b, float c, float d, float e);
bool roots_quartic_equation2(double a, double b, double c, double d, double e, float (&roots)[5]);

double* roots_quadratic_equation_rain(double a, double b, double c);
double* roots_cubic_equation_rain(double a, double b, double c, double d);
bool roots_quartic_equation2_rain(double a, double b, double c, double d, double e, double (&roots)[5]);

bool get_matrix_3_1(const vec3 a, float (&m)[3][1]);
bool get_matrix_3_3(const vec3 a, const vec3 b, const vec3 c, float (&m)[3][3]);
bool get_matrix_inverse_3_3(const float m[3][3], float (&inverse)[3][3]);
bool matrix_3_3_multiply_3_1(const float m1[3][3], const float m2[3][1], float (&m)[3][1]);
bool matrix_3_1_minus_3_1(const float m1[3][1], const float m2[3][1], float (&m)[3][1]);

bool get_bezier_matrix_xyz_4_4(const vec3 point[4][4], float (&bezier_x)[4][4], float (&bezier_y)[4][4], float (&bezier_z)[4][4]);
bool matrix_4_4_multiply_4_4(const float matrix1[4][4], const float matrix2[4][4], float (&result)[4][4]);
bool matrix_4_4_multiply_4_1(const float matrix1[4][4], const float matrix2[4][1], float (&result)[4][1]);
bool get_matrix_transpose_4_4(const float matrix[4][4], float (&result)[4][4]);

bool get_teapot_data(int (&patches)[32][16], float (&vertices)[306][3]);

bool roots_num_equation_6th(float ee6[7], double a, double b, int &num);
bool roots_equation_6th(float ee6[7], float a, float b, float tol, float (&roots)[7]);
bool roots_num_equation_10th(float ee10[11], double a, double b, int &num);
bool roots_equation_10th(float ee10[11], float a, float b, float tol, float (&roots)[11]);
bool roots_equation_10th_bt(float ee10[11], float a, float b, float tol, float (&roots)[11]);

#endif // VEC3_H
