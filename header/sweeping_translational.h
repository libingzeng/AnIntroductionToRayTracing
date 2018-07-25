#ifndef SWEEPING_TRANSLATIONAL_H
#define SWEEPING_TRANSLATIONAL_H

#include "hitable.h"
#include "log.h"

class sweeping_translational : public hitable
{
    public:
        sweeping_translational() {}
        sweeping_translational(vec3 *cp, float by, float cy, material *m, int ty){
/*
ty=1: translational sweeping
ty=2: conic sweeping
*/
            base_y = by;
            cap_y = cy;
            ma = m;
            type = ty;
            float matrix_t_6[4][4] = {{ 1,  4,  1, 0},
                                      {-3,  0,  3, 0},
                                      { 3, -6,  3, 0},
                                      {-1,  3, -3, 1}};
            float matrix_t[4][4];
            for (int i=0; i<4; i++) {
                for (int j=0; j<4; j++) {
                    matrix_t[i][j] = matrix_t_6[i][j] / 6.0;
                }
            }
            float points_x[17], points_z[17], b_points_x[4][1], b_points_z[4][1], b_result_x[4][1], b_result_z[4][1];
            int b_num = 0;
            float min_x, max_x, min_z, max_z;
            for (int i=0; i<17; i++) {
                points_x[i] = cp[i].x();
                points_z[i] = cp[i].z();
            }
            min_x = points_x[0];
            max_x = points_x[0];
            min_z = points_z[0];
            max_z = points_z[0];
            for (int i=0; i<17; i++) {
                if (min_x > points_x[i]) {
                    min_x = points_x[i];
                }
                if (max_x < points_x[i]) {
                    max_x = points_x[i];
                }
                if (min_z > points_z[i]) {
                    min_z = points_z[i];
                }
                if (max_z < points_z[i]) {
                    max_z = points_z[i];
                }
            }
            sweeping_bl = vec3(min_x, base_y, min_z);
            sweeping_bh = vec3(max_x, cap_y, max_z);

            for (int i=0; i<14; i=i+1) {
                b_points_x[0][0] = points_x[i];
                b_points_x[1][0] = points_x[i+1];
                b_points_x[2][0] = points_x[i+2];
                b_points_x[3][0] = points_x[i+3];
                b_points_z[0][0] = points_z[i];
                b_points_z[1][0] = points_z[i+1];
                b_points_z[2][0] = points_z[i+2];
                b_points_z[3][0] = points_z[i+3];
                matrix_4_4_multiply_4_1(matrix_t, b_points_x, b_result_x);
                matrix_4_4_multiply_4_1(matrix_t, b_points_z, b_result_z);
                for (int j=0; j<4; j++) {
                    matrix_c_x[j][b_num] = ((fabs(b_result_x[j][0])<1e-6)? (0.0):(b_result_x[j][0]));
                    matrix_c_z[j][b_num] = ((fabs(b_result_z[j][0])<1e-6)? (0.0):(b_result_z[j][0]));
                }
                b_num ++;
            }

        }
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        float matrix_c_x[4][14], matrix_c_z[4][14];
        vec3 sweeping_bl, sweeping_bh;
        float base_y, cap_y;
        material *ma;
        int type;
};

#endif // SWEEPING_TRANSLATIONAL_H
