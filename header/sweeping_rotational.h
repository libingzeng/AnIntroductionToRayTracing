#ifndef SWEEPING_ROTATIONAL_H
#define SWEEPING_ROTATIONAL_H

#include "hitable.h"
#include "log.h"

class sweeping_rotational : public hitable
{
    public:
        sweeping_rotational() {}
        sweeping_rotational(vec3 cen, vec3 *cp, material *m, bool c){
            center = cen;
            ma = m;
            cut = c;
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
            float points_u[8], points_v[8], b_points_u[4][1], b_points_v[4][1], b_result_u[4][1], b_result_v[4][1];
            int b_num = 0;
            float max_u, min_v, max_v;
            for (int i=0; i<8; i++) {//8
                points_u[i] = cp[i].x();
                points_v[i] = cp[i].y();
            }
            max_u = points_u[0];
            min_v = points_v[0];
            max_v = points_v[0];
            for (int i=0; i<8; i++) {//8
                if (max_u < points_u[i]) {
                    max_u = points_u[i];
                }
                if (min_v > points_v[i]) {
                    min_v = points_v[i];
                }
                if (max_v < points_v[i]) {
                    max_v = points_v[i];
                }
            }
            sweeping_bl = vec3(-max_u+center.x(), min_v+center.y(), -max_u+center.z());
            sweeping_bh = vec3( max_u+center.x(), max_v+center.y(),  max_u+center.z());

            for (int i=0; i<5; i=i+1) {//5
                b_points_u[0][0] = points_u[i];
                b_points_u[1][0] = points_u[i+1];
                b_points_u[2][0] = points_u[i+2];
                b_points_u[3][0] = points_u[i+3];
                b_points_v[0][0] = points_v[i];
                b_points_v[1][0] = points_v[i+1];
                b_points_v[2][0] = points_v[i+2];
                b_points_v[3][0] = points_v[i+3];
                matrix_4_4_multiply_4_1(matrix_t, b_points_u, b_result_u);
                matrix_4_4_multiply_4_1(matrix_t, b_points_v, b_result_v);
                for (int j=0; j<4; j++) {
                    matrix_c_u[j][b_num] = ((fabs(b_result_u[j][0])<1e-6)? (0.0):(b_result_u[j][0]));
                    matrix_c_v[j][b_num] = ((fabs(b_result_v[j][0])<1e-6)? (0.0):(b_result_v[j][0]));
                }
                b_num ++;
            }

        }
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        float matrix_c_u[4][5], matrix_c_v[4][5];
        vec3 sweeping_bl, sweeping_bh, center;
        material *ma;
        bool cut;
};

#endif // SWEEPING_ROTATIONAL_H
