#ifndef PARAMETRIC_SURFACE_H
#define PARAMETRIC_SURFACE_H

#include "hitable.h"
#include "material.h"
#include "log.h"

#define TYPE 3
/*
TYPE=1: sphere;
TYPE=2: horn;
TYPE=3: bezier3;
*/
#define CURVE 2
/*
CURVE=1: bezier
CURVE=2: b-spline
*/

class parametric_surface : public hitable
{
    public:
        parametric_surface() {}
#if ((TYPE == 1) || (TYPE ==2))
        parametric_surface(vec3 cen, float a, float b, float c, float hy, material *m, float r, float t, float n) : center(cen), intercept_x(a), intercept_y(b), intercept_z(c), height_half_y(hy), ma(m), rho(r), theta(t), num(n) {}
#endif // TYPE
#if TYPE == 3
        parametric_surface(material *m, float r, float t, float n, vec3 *cp) : ma(m), rho(r), theta(t), num(n) {
            vec3 ctrl_points[4][4];
            for (int i=0; i<16; i++) {
                ctrl_points[(i/4)][(i%4)] = cp[i];
            }
#if CURVE == 1//bezier
            float matrix_t[4][4] = {{ 1,  0,  0, 0},
                                    {-3,  3,  0, 0},
                                    { 3, -6,  3, 0},
                                    {-1,  3, -3, 1}};
            float matrix[4][4] = {{1, -3,  3, -1},
                                  {0,  3, -6,  3},
                                  {0,  0,  3, -3},
                                  {0,  0,  0,  1}};
#endif // CURVE
#if CURVE == 2//b-spline
            float matrix_t_6[4][4] = {{ 1,  4,  1, 0},
                                      {-3,  0,  3, 0},
                                      { 3, -6,  3, 0},
                                      {-1,  3, -3, 1}};
            float matrix_t[4][4], matrix[4][4];
            for (int i=0; i<4; i++) {
                for (int j=0; j<4; j++) {
                    matrix_t[i][j] = matrix_t_6[i][j] / 6.0;
                    matrix[j][i] = matrix_t[i][j];
                }
            }
#endif // CURVE
            float points_x[4][4], points_y[4][4], points_z[4][4], points_x_t[4][4], points_y_t[4][4], points_z_t[4][4];
            get_bezier_matrix_xyz_4_4(ctrl_points, points_x, points_y, points_z);
            matrix_4_4_multiply_4_4(matrix_t, points_x, points_x_t);
            matrix_4_4_multiply_4_4(points_x_t, matrix, matrix_c_x);
            matrix_4_4_multiply_4_4(matrix_t, points_y, points_y_t);
            matrix_4_4_multiply_4_4(points_y_t, matrix, matrix_c_y);
            matrix_4_4_multiply_4_4(matrix_t, points_z, points_z_t);
            matrix_4_4_multiply_4_4(points_z_t, matrix, matrix_c_z);
            float min_x = points_x[0][0];
            float max_x = points_x[0][0];
            float min_y = points_y[0][0];
            float max_y = points_y[0][0];
            float min_z = points_z[0][0];
            float max_z = points_z[0][0];
            for (int i=0; i<4; i++) {
                for (int j=0; j<4; j++) {
                    if (min_x > points_x[i][j]) {
                        min_x = points_x[i][j];
                    }
                    if (max_x < points_x[i][j]) {
                        max_x = points_x[i][j];
                    }
                    if (min_y > points_y[i][j]) {
                        min_y = points_y[i][j];
                    }
                    if (max_y < points_y[i][j]) {
                        max_y = points_y[i][j];
                    }
                    if (min_z > points_z[i][j]) {
                        min_z = points_z[i][j];
                    }
                    if (max_z < points_z[i][j]) {
                        max_z = points_z[i][j];
                    }
                }
            }
            bezier3_bl = vec3(min_x, min_y, min_z);
            bezier3_bh = vec3(max_x, max_y, max_z);
        }
#endif // TYPE
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
#if ((TYPE == 1) || (TYPE ==2))
        vec3 center;
        float intercept_x;
        float intercept_y;
        float intercept_z;
        float height_half_y;
#endif // TYPE
        material *ma;
        float rho, theta, num;
#if TYPE == 3
        float matrix_c_x[4][4], matrix_c_y[4][4], matrix_c_z[4][4];
        vec3 bezier3_bl, bezier3_bh;
#endif // TYPE
};

#endif // PARAMETRIC_SURFACE_H
