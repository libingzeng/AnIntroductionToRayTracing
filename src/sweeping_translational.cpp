#include "sweeping_translational.h"

#include <iostream>
#include <limits>
#include "float.h"

bool ray_hit_box_st(const ray& r, const vec3& vertex_l, const vec3& vertex_h, float& t_near, float& t_far) {
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

bool sweeping_translational::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if SWEEPING_TRANSLATIONAL_LOG == 1
        std::cout << "-------------sweeping_translational::hit---1-------------" << endl;
#endif // SWEEPING_TRANSLATIONAL_LOG
        float t_near, t_far;
        if (ray_hit_box_st(r, sweeping_bl, sweeping_bh, t_near, t_far)) {
            float x0 = r.origin().x();
            float y0 = r.origin().y();
            float z0 = r.origin().z();
            float xd = r.direction().x();
            float yd = r.direction().y();
            float zd = r.direction().z();
            if (fabs(xd) < 1e-6) {
                xd = 1e-6;
            }
            if (fabs(yd) < 1e-6) {
                yd = 1e-6;
            }
            if (fabs(zd) < 1e-6) {
                zd = 1e-6;
            }
            float cos_theta = dot(r.direction(), vec3(0, -1, 0))/(r.direction().length());
            float t0, t1, t0_t1_smaller, t0_t1_bigger, a3, a2, a1, a0, *roots3_base, roots_base[21][3], temp0, temp1, temp2, root[5], u, xu, rbt, rbt_s, rbt_b, check_y;
            roots_base[0][0] = 0.0;
            root[0] = 0.0;
            root[2] = 0.0;
            root[3] = -1.0;
            root[4] = -1.0;
            int num_roots_base = 0;
            int num_t0_t1_smaller = 0;
            int num_t0_t1_bigger = 0;
            int num_roots_base_valid = 0;
            t0 = (base_y - y0) / yd;
            t1 = (cap_y - y0) / yd;
            t0_t1_smaller = (t0>t1)? t1:t0;
            t0_t1_bigger = (t0>t1)? t0:t1;
            for (int i=0; i<14; i++) {
                if (type == 1) {
                    a3 = zd*matrix_c_x[3][i] - xd*matrix_c_z[3][i];
                    a2 = zd*matrix_c_x[2][i] - xd*matrix_c_z[2][i];
                    a1 = zd*matrix_c_x[1][i] - xd*matrix_c_z[1][i];
                    a0 = zd*matrix_c_x[0][i] - xd*matrix_c_z[0][i] + xd*z0 - zd*x0;
                }
                if (type == 2) {
                    a3 = ((y0-base_y)*zd-z0*yd)*matrix_c_x[3][i] - ((y0-base_y)*xd-x0*yd)*matrix_c_z[3][i];
                    a2 = ((y0-base_y)*zd-z0*yd)*matrix_c_x[2][i] - ((y0-base_y)*xd-x0*yd)*matrix_c_z[2][i];
                    a1 = ((y0-base_y)*zd-z0*yd)*matrix_c_x[1][i] - ((y0-base_y)*xd-x0*yd)*matrix_c_z[1][i];
                    a0 = ((y0-base_y)*zd-z0*yd)*matrix_c_x[0][i] - ((y0-base_y)*xd-x0*yd)*matrix_c_z[0][i] + (cap_y-base_y)*(xd*z0 - zd*x0);
                }
                roots3_base = roots_cubic_equation(a3, a2, a1, a0);
                if (roots3_base[0] != 0.0) {
                    for (int j=1; j<(int(roots3_base[0])+1); j++) {
                        u = roots3_base[j];
                        if ((u>=0.0) && (u<=1.0)) {
                            xu = matrix_c_x[0][i]+matrix_c_x[1][i]*u+matrix_c_x[2][i]*u*u+matrix_c_x[3][i]*u*u*u;
                            if (type == 1) {
                                roots_base[num_roots_base+1][0] = (xu-x0)/xd;
                                roots_base[num_roots_base+1][1] = i;
                                roots_base[num_roots_base+1][2] = u;
                                num_roots_base ++;
                            }
                            if (type == 2) {
//                                roots_base[num_roots_base+1][0] = (xu-x0)/xd;
                                rbt = ((y0-base_y)*xu-(cap_y-base_y)*x0)/((cap_y-base_y)*xd-yd*xu);
                                if (yd > 0) {
                                    rbt_s = (base_y-y0)/yd;
                                    rbt_b = (cap_y-y0)/yd;
                                }
                                if (yd < 0) {
                                    rbt_s = (cap_y-y0)/yd;
                                    rbt_b = (base_y-y0)/yd;
                                }
                                if ((rbt>rbt_s) && (rbt<rbt_b)) {
//                                if (1) {
                                    roots_base[num_roots_base+1][0] = rbt;
                                    roots_base[num_roots_base+1][1] = i;
                                    roots_base[num_roots_base+1][2] = u;
                                    num_roots_base ++;
                                }
                            }
                        }
                    }
                    roots_base[0][0] = float(num_roots_base);
                }
                delete [] roots3_base;
            }
            if (roots_base[0][0] != 0.0) {
                for (int i=1; i<int(roots_base[0][0]); i++) {
                    for (int j=i+1; j<int(roots_base[0][0])+1; j++) {
                        if (roots_base[i][0] > roots_base[j][0]) {
                            temp0 = roots_base[i][0];
                            roots_base[i][0] = roots_base[j][0];
                            roots_base[j][0] = temp0;
                            temp1 = roots_base[i][1];
                            roots_base[i][1] = roots_base[j][1];
                            roots_base[j][1] = temp1;
                            temp2 = roots_base[i][2];
                            roots_base[i][2] = roots_base[j][2];
                            roots_base[j][2] = temp2;
                        }
                    }
                }
                for (int i=1; i<(int(roots_base[0][0])+1); i++) {
                    if (roots_base[i][0] > (t0_t1_smaller)) {
                        num_t0_t1_smaller = int(roots_base[0][0])+1-i;
                        break;
                    }
                }
                if ((num_t0_t1_smaller%2) != 0) {
                    root[1] = t0_t1_smaller;
                    root[0] = 1.0;
                    root[2] = 1.0;
                }
                else {
                    for (int i=1; i<(int(roots_base[0][0])+1); i++) {
                        if (roots_base[i][0] > (t0_t1_bigger)) {
                            num_t0_t1_bigger = int(roots_base[0][0])+1-i;
                            break;
                        }
                    }
                    if ((num_t0_t1_bigger%2) != 0) {
                        root[1] = roots_base[1][0];//fabs(cos_theta);
                        root[0] = 1.0;
                        root[2] = 2.0;
                        root[3] = roots_base[1][1];
                        root[4] = roots_base[1][2];
                    }
                    else {
                        for (int i=1; i<(int(roots_base[0][0])+1); i++) {
                            if (fabs(roots_base[i][0]-t0_t1_bigger)*fabs(dot(r.direction(), vec3(0, -1, 0)))<fabs(cap_y-base_y)) {
//                            if (roots_base[i][0] > t0_t1_smaller) {
                                num_roots_base_valid = i;
                                break;
                            }
                        }
                        if (num_roots_base_valid != 0) {
                            root[1] = roots_base[num_roots_base_valid][0];//fabs(cos_theta);
                            check_y = y0 + yd * root[1];
                            if ((check_y > base_y) && (check_y < cap_y)) {
                                root[0] = 1.0;
                                root[2] = 3.0;
                                root[3] = roots_base[num_roots_base_valid][1];
                                root[4] = roots_base[num_roots_base_valid][2];
                            }
                        }
                    }
                }
                if (root[0] != 0.0) {
                    if ((root[1] < t_max) && (root[1] > t_min)) {
                        rec.t = root[1];
                        rec.p = r.point_at_parameter(rec.t);
                        if (root[2] == 1.0) {
                            rec.normal = vec3(0, 1, 0);
                        }
                        else {
                            float nx, nz, nu;
                            int num = int(root[3]);
                            nu = root[4];
                            nx = matrix_c_x[1][num]+2.0*matrix_c_x[2][num]*nu+3.0*matrix_c_x[3][num]*nu*nu;
                            nz = matrix_c_z[1][num]+2.0*matrix_c_z[2][num]*nu+3.0*matrix_c_z[3][num]*nu*nu;
                            if (type == 1) {
                                rec.normal = unit_vector(vec3(-nx, 0, nz));
                            }
                            if (type == 2) {
                                vec3 va = vec3(nx, 0, nz);
                                float bx = matrix_c_x[0][num]+matrix_c_x[1][num]*nu+matrix_c_x[2][num]*nu*nu+matrix_c_x[3][num]*nu*nu*nu;
                                float bz = matrix_c_z[0][num]+matrix_c_z[1][num]*nu+matrix_c_z[2][num]*nu*nu+matrix_c_z[3][num]*nu*nu*nu;
                                vec3 vb = vec3(bx, (cap_y-base_y), bz);
                                rec.normal = unit_vector(cross(va, vb));
                            }
                        }
                        if(dot(r.direction(), rec.normal) > 0) {
                            rec.normal = - rec.normal;
                        }
                        rec.mat_ptr = ma;
                        rec.u = -1.0;
                        rec.v = -1.0;
                        return true;
                    }
                    else {
                        return false;
                    }
                }
                else {
                    return false;
                }
            }
            else {
                return false;
            }
        }
        else {
            return false;
        }
}
