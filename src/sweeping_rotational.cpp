#include "sweeping_rotational.h"

#include <iostream>
#include <limits>
#include "float.h"

bool ray_hit_box_sr(const ray& r, const vec3& vertex_l, const vec3& vertex_h, float& t_near, float& t_far) {
        t_near = (numeric_limits<float>::min)();
        t_far = (numeric_limits<float>::max)();
        vec3 direction = r.direction();
        vec3 origin = r.origin();
        vec3 bl = vertex_l;
        vec3 bh = vertex_h;
        float array1[6];

        if(direction.x() == 0) {
            if((origin.x() < bl.x()) || (origin.x() > bh.x())) {
                return false;
            }
            array1[0] = (numeric_limits<float>::min)();
            array1[1] = (numeric_limits<float>::max)();
        }
        if(direction.y() == 0) {
            if((origin.y() < bl.y()) || (origin.y() > bh.y())) {
                return false;
            }
            array1[2] = (numeric_limits<float>::min)();
            array1[3] = (numeric_limits<float>::max)();
        }
        if(direction.z() == 0) {
            if((origin.z() < bl.z()) || (origin.z() > bh.z())) {
                return false;
            }
            array1[4] = (numeric_limits<float>::min)();
            array1[5] = (numeric_limits<float>::max)();
        }

        if((direction.x() != 0) && (direction.y() != 0) && (direction.z() != 0)) {
            array1[0] = (bl.x()-origin.x())/direction.x();
            array1[1] = (bh.x()-origin.x())/direction.x();
            array1[2] = (bl.y()-origin.y())/direction.y();
            array1[3] = (bh.y()-origin.y())/direction.y();
            array1[4] = (bl.z()-origin.z())/direction.z();
            array1[5] = (bh.z()-origin.z())/direction.z();
        }

        for (int i=0; i<6; i=i+2){
            if(array1[i] > array1[i+1]) {
                float t = array1[i];
                array1[i] = array1[i+1];
                array1[i+1] = t;
            }
            if(array1[i] >= t_near) {t_near = array1[i];}
            if(array1[i+1] <= t_far) {t_far = array1[i+1];}
            if(t_near > t_far) {
                return false;
            }
            if(t_far < 0) {
                return false;
            }
        }
        if (t_near != t_near) {
            t_near = t_near * 1;
        }
        return true;
}

bool sweeping_rotational::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if SWEEPING_ROTATIONAL_LOG == 1
        std::cout << "-------------sweeping_rotational::hit---1-------------" << endl;
#endif // SWEEPING_ROTATIONAL_LOG
        float t_near, t_far;
        if (ray_hit_box_sr(r, sweeping_bl, sweeping_bh, t_near, t_far)) {
            float xx0 = r.origin().x() - center.x();
            float yy0 = r.origin().y() - center.y();
            float zz0 = r.origin().z() - center.z();
            float xxd = r.direction().x();
            float yyd = r.direction().y();
            float zzd = r.direction().z();
            if (fabs(xxd) < 1e-6) {
                xxd = 1e-6;
            }
            if (fabs(yyd) < 1e-6) {
                yyd = 1e-6;
            }
            if (fabs(zzd) < 1e-6) {
                zzd = 1e-6;
            }
            float aaa = xxd*xxd + zzd*zzd;
            float bbb = 2 * ((xx0*xxd+zz0*zzd)*yyd - aaa*yy0);
            float ccc = (xx0*yyd-xxd*yy0)*(xx0*yyd-xxd*yy0) + (zz0*yyd-zzd*yy0)*(zz0*yyd-zzd*yy0);
            float ddd = yyd*yyd;
            float ss6[7], roots[7], roots_t[19][3], yyv, temp0, temp1, temp2;
            float tol = 1e-6;
            int num_roots_t = 0;

            for (int i=0; i<5; i++) {//5
                ss6[0] = aaa*matrix_c_v[3][i]*matrix_c_v[3][i] - ddd*matrix_c_u[3][i]*matrix_c_u[3][i];
                ss6[1] = aaa*2*matrix_c_v[2][i]*matrix_c_v[3][i] - ddd*2*matrix_c_u[2][i]*matrix_c_u[3][i];
                ss6[2] = aaa*(matrix_c_v[2][i]*matrix_c_v[2][i]+2*matrix_c_v[1][i]*matrix_c_v[3][i]) -
                      ddd*(matrix_c_u[2][i]*matrix_c_u[2][i]+2*matrix_c_u[1][i]*matrix_c_u[3][i]);
                ss6[3] = aaa*2*(matrix_c_v[0][i]*matrix_c_v[3][i]+matrix_c_v[1][i]*matrix_c_v[2][i]) + bbb*matrix_c_v[3][i] -
                      ddd*2*(matrix_c_u[0][i]*matrix_c_u[3][i]+matrix_c_u[1][i]*matrix_c_u[2][i]);
                ss6[4] = aaa*(matrix_c_v[1][i]*matrix_c_v[1][i]+2*matrix_c_v[0][i]*matrix_c_v[2][i]) + bbb*matrix_c_v[2][i] -
                      ddd*(matrix_c_u[1][i]*matrix_c_u[1][i]+2*matrix_c_u[0][i]*matrix_c_u[2][i]);
                ss6[5] = aaa*2*matrix_c_v[0][i]*matrix_c_v[1][i] + bbb*matrix_c_v[1][i]  - ddd*2*matrix_c_u[0][i]*matrix_c_u[1][i];
                ss6[6] = aaa*matrix_c_v[0][i]*matrix_c_v[0][i] + bbb*matrix_c_v[0][i]  - ddd*matrix_c_u[0][i]*matrix_c_u[0][i] + ccc;
                roots_equation_6th(ss6, 0, 1, tol, roots);
                for (int j=1; j<(int(roots[0])+1); j++) {
                    yyv = matrix_c_v[0][i]+matrix_c_v[1][i]*roots[j]+matrix_c_v[2][i]*roots[j]*roots[j]+matrix_c_v[3][i]*roots[j]*roots[j]*roots[j];
                    roots_t[num_roots_t+1][0] = (yyv-yy0)/yyd;
                    if (cut) {
                        rec.t = roots_t[num_roots_t+1][0];
                        rec.p = r.point_at_parameter(rec.t);
                        if (fabs(rec.p.z())<0.5) {//this restrict can help cut strips from revolution
                            roots_t[num_roots_t+1][1] = i;
                            roots_t[num_roots_t+1][2] = roots[j];
                            num_roots_t ++;
                        }
                    }
                    else {
                        roots_t[num_roots_t+1][1] = i;
                        roots_t[num_roots_t+1][2] = roots[j];
                        num_roots_t ++;
                    }
                }
            }
            roots_t[0][0] = float(num_roots_t);

            if (roots_t[0][0] != 0.0) {
                for (int i=1; i<int(roots_t[0][0]); i++) {
                    for (int j=i+1; j<int(roots_t[0][0])+1; j++) {
                        if (roots_t[i][0] > roots_t[j][0]) {
                            temp0 = roots_t[i][0];
                            roots_t[i][0] = roots_t[j][0];
                            roots_t[j][0] = temp0;
                            temp1 = roots_t[i][1];
                            roots_t[i][1] = roots_t[j][1];
                            roots_t[j][1] = temp1;
                            temp2 = roots_t[i][2];
                            roots_t[i][2] = roots_t[j][2];
                            roots_t[j][2] = temp2;
                        }
                    }
                }
                if ((roots_t[1][0] < t_max) && (roots_t[1][0] > t_min)) {
                    rec.t = roots_t[1][0];
                    rec.p = r.point_at_parameter(rec.t);
                    float nx = rec.p.x() - center.x();
                    float nz = rec.p.z() - center.z();
                    int num = int(roots_t[1][1]);
                    float sss = roots_t[1][2];
                    float nu = -(matrix_c_v[1][num]+2.0*matrix_c_v[2][num]*sss+3.0*matrix_c_v[3][num]*sss*sss);
                    float nv =  (matrix_c_u[1][num]+2.0*matrix_c_u[2][num]*sss+3.0*matrix_c_u[3][num]*sss*sss);
                    float ny = sqrt(nx*nx+nz*nz)*nv/nu - center.y();
                    rec.normal = unit_vector(vec3(nx, ny, nz));
                    if(dot(r.direction(), rec.normal) > 0) {
                        rec.normal = - rec.normal;
                    }
                    rec.mat_ptr = ma;
                    rec.u = -1.0;
                    rec.v = -1.0;
                    return true;
                }
            }
        }
        return false;
}
