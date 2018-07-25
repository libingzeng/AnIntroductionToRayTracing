#ifndef SWEEPING_SPHERE_H
#define SWEEPING_SPHERE_H

#include "hitable.h"
#include "log.h"

class sweeping_sphere : public hitable
{
    public:
        sweeping_sphere() {}
        sweeping_sphere(vec3 *ccp, vec3 *rcp, material *m){
            ma = m;
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
            float points_cx[8], points_cy[8], points_cz[8], points_r[8];
            float b_points_cx[4][1], b_points_cy[4][1], b_points_cz[4][1], b_points_r[4][1];
            float b_result_cx[4][1], b_result_cy[4][1], b_result_cz[4][1], b_result_r[4][1];
            int b_num = 0;
            float min_cx, max_cx, min_cy, max_cy, min_cz, max_cz, max_r;
            for (int i=0; i<8; i++) {//8
                points_cx[i] = ccp[i].x();
                points_cy[i] = ccp[i].y();
                points_cz[i] = ccp[i].z();
                points_r[i]  = rcp[i].x();
            }
            min_cx = points_cx[0];
            max_cx = points_cx[0];
            min_cy = points_cy[0];
            max_cy = points_cy[0];
            min_cz = points_cz[0];
            max_cz = points_cz[0];
            max_r  = points_r[0];
            for (int i=0; i<8; i++) {//8
                if (min_cx > points_cx[i]) {
                    min_cx = points_cx[i];
                }
                if (max_cx < points_cx[i]) {
                    max_cx = points_cx[i];
                }
                if (min_cy > points_cy[i]) {
                    min_cy = points_cy[i];
                }
                if (max_cy < points_cy[i]) {
                    max_cy = points_cy[i];
                }
                if (min_cz > points_cz[i]) {
                    min_cz = points_cz[i];
                }
                if (max_cz < points_cz[i]) {
                    max_cz = points_cz[i];
                }
               if (max_r < points_r[i]) {
                    max_r = points_r[i];
                }
            }
            sweeping_bl = vec3(min_cx-max_r, min_cy-max_r, min_cz-max_r);
            sweeping_bh = vec3(max_cx+max_r, max_cy+max_r, max_cz+max_r);

            for (int i=0; i<5; i=i+1) {//5
                b_points_cx[0][0] = points_cx[i];
                b_points_cx[1][0] = points_cx[i+1];
                b_points_cx[2][0] = points_cx[i+2];
                b_points_cx[3][0] = points_cx[i+3];
                b_points_cy[0][0] = points_cy[i];
                b_points_cy[1][0] = points_cy[i+1];
                b_points_cy[2][0] = points_cy[i+2];
                b_points_cy[3][0] = points_cy[i+3];
                b_points_cz[0][0] = points_cz[i];
                b_points_cz[1][0] = points_cz[i+1];
                b_points_cz[2][0] = points_cz[i+2];
                b_points_cz[3][0] = points_cz[i+3];
                b_points_r[0][0]  = points_r[i];
                b_points_r[1][0]  = points_r[i+1];
                b_points_r[2][0]  = points_r[i+2];
                b_points_r[3][0]  = points_r[i+3];
                matrix_4_4_multiply_4_1(matrix_t, b_points_cx, b_result_cx);
                matrix_4_4_multiply_4_1(matrix_t, b_points_cy, b_result_cy);
                matrix_4_4_multiply_4_1(matrix_t, b_points_cz, b_result_cz);
                matrix_4_4_multiply_4_1(matrix_t, b_points_r, b_result_r);
                for (int j=0; j<4; j++) {
                    matrix_c_cx[j][b_num] = ((fabs(b_result_cx[j][0])<1e-6)? (0.0):(b_result_cx[j][0]));
                    matrix_c_cy[j][b_num] = ((fabs(b_result_cy[j][0])<1e-6)? (0.0):(b_result_cy[j][0]));
                    matrix_c_cz[j][b_num] = ((fabs(b_result_cz[j][0])<1e-6)? (0.0):(b_result_cz[j][0]));
                    matrix_c_r[j][b_num]  = ((fabs(b_result_r[j][0])<1e-6)? (0.0):(b_result_r[j][0]));
                }
                b_num ++;
            }

        }
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        float matrix_c_cx[4][5], matrix_c_cy[4][5], matrix_c_cz[4][5], matrix_c_r[4][5];
        vec3 sweeping_bl, sweeping_bh;
        material *ma;
};

#endif // SWEEPING_SPHERE_H
