#include "quadratic.h"

#include <iostream>
using namespace std;

bool quadratic::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if QUADRATIC_LOG == 1
        std::cout << "-------------quadratic::hit----------------" << endl;
#endif // QUADRATIC_LOG
        float ab_square = intercept_x*intercept_x*intercept_y*intercept_y;
        float bc_square = intercept_y*intercept_y*intercept_z*intercept_z;
        float ac_square = sign1*intercept_x*intercept_x*intercept_z*intercept_z;
        float abc_square = sign2*intercept_x*intercept_x*intercept_y*intercept_y*intercept_z*intercept_z;

        vec3 inter_square = vec3(bc_square, ac_square, ab_square);
        vec3 rd_square = vec3(r.direction().x()*r.direction().x(),
                              r.direction().y()*r.direction().y(),
                              r.direction().z()*r.direction().z());
        float A = dot(inter_square, rd_square);
        vec3 r0_c = r.origin() - center;
        vec3 r0_c_rd = vec3(r0_c.x()*r.direction().x(),
                            r0_c.y()*r.direction().y(),
                            r0_c.z()*r.direction().z());
        float B = 2*dot(r0_c_rd, inter_square);
        vec3 r0_c_square = vec3(r0_c.x()*r0_c.x(),
                                r0_c.y()*r0_c.y(),
                                r0_c.z()*r0_c.z());
        float C = dot(r0_c_square, inter_square) - abc_square;
        float temp, temp1, temp2;
        vec3 pc;

        if(A == 0) {
            if (B == 0) {
                return false;
            }
            else {
                temp = -C/B;
                if (temp < t_max && temp > t_min) {
                    rec.t = temp;
                    rec.p = r.point_at_parameter(rec.t);
                    if (((rec.p.y()-center.y()) > -height_half_y) && ((rec.p.y()-center.y()) < height_half_y)) {
                        pc = rec.p - center;
                        rec.normal = unit_vector(vec3(2*bc_square*pc.x(), 2*ac_square*pc.y(), 2*ab_square*pc.z()));
                        if (dot(rec.normal, r.direction()) > 0) {
                            rec.normal = -rec.normal;
                        }
                        rec.mat_ptr = ma;
                        if ((sign1 == 0) && (sign2 == 1) && (intercept_x == intercept_y) && (intercept_y == intercept_z)) {//0/1:cylinder
                            rec.v = (rec.p.y() - (center.y() - height_half_y)) / (2*height_half_y);
                            vec3 pc1 = vec3(rec.p.x()-center.x(), 0, rec.p.z()-center.z());
                            vec3 vx = vec3(0, 0, 1);
                            float u = acos(dot(pc1, vx) / (pc1.length()*vx.length())) / (2*M_PI);
//                            float u = acos((rec.p.x() - center.x()) / intercept_x) / (2*M_PI);
                            if ((rec.p.x() - center.x()) < 0) {
                                rec.u = 1-u;
                            }
                            else {
                                rec.u = u;
                            }
                            return true;
                        }
                        else if ((sign1 == -1) && (sign2 == 0) && (intercept_x == intercept_z)) {//-1/0:cone
                            if (((rec.p.y()-center.y()) > -height_half_y) && ((rec.p.y()-center.y()) < 0)) {
                                pc = rec.p - center;
                                rec.normal = unit_vector(vec3(2*bc_square*pc.x(), 2*ac_square*pc.y(), 2*ab_square*pc.z()));
                                if (dot(rec.normal, r.direction()) > 0) {
                                    rec.normal = -rec.normal;
                                }
                                rec.mat_ptr = ma;

                                rec.v = (rec.p.y() - (center.y() - height_half_y)) / (height_half_y);
//                                float rate = intercept_y / intercept_x;
                                vec3 pc1 = vec3(rec.p.x()-center.x(), 0, rec.p.z()-center.z());
                                vec3 vx = vec3(0, 0, 1);
                                float u = acos(dot(pc1, vx) / (pc1.length()*vx.length())) / (2*M_PI);
        //                        float u = acos((rec.p.x() - center.x()) / ((center.y() - rec.p.y()) / rate)) / (2*M_PI);
                                if ((rec.p.x() - center.x()) < 0) {
                                   rec.u = 1-u;
                                }
                                else {
                                    rec.u = u;
                                }
                                return true;
                            }
                            else {
                                rec.u = -1.0;
                                rec.v = -1.0;
                            }
                       }
                        else {
                            rec.u = -1.0;
                            rec.v = -1.0;
                            return true;
                        }
                    }
                    else {
                        return false;
                    }
                }
            }
        }
        else {
            float discriminant = B*B - 4*A*C;
            if (discriminant >= 0) {
                temp1 = (-B - sqrt(discriminant)) / (2.0*A);
                temp2 = (-B + sqrt(discriminant)) / (2.0*A);
                if (temp1 > temp2) {//make sure that temp1 is smaller than temp2
                    temp = temp1;
                    temp1 = temp2;
                    temp2 = temp;
                }
                if (temp1 < t_max && temp1 > t_min) {
                    rec.t = temp1;
                    rec.p = r.point_at_parameter(rec.t);
                    if (((rec.p.y()-center.y()) > -height_half_y) && ((rec.p.y()-center.y()) < height_half_y)) {
                        pc = rec.p - center;
                        rec.normal = unit_vector(vec3(2*bc_square*pc.x(), 2*ac_square*pc.y(), 2*ab_square*pc.z()));
                        if (dot(rec.normal, r.direction()) > 0) {
                            rec.normal = -rec.normal;
                        }
                        rec.mat_ptr = ma;
                        if ((sign1 == 0) && (sign2 == 1) && (intercept_x == intercept_y) && (intercept_y == intercept_z)) {//cylinder
                            rec.v = (rec.p.y() - (center.y() - height_half_y)) / (2*height_half_y);
                            vec3 pc1 = vec3(rec.p.x()-center.x(), 0, rec.p.z()-center.z());
                            vec3 vx = vec3(0, 0, 1);
                            float u = acos(dot(pc1, vx) / (pc1.length()*vx.length())) / (2*M_PI);
//                            float u = acos((rec.p.x() - center.x()) / intercept_x) / (2*M_PI);
                            if ((rec.p.x() - center.x()) < 0) {
                                rec.u = 1-u;
                            }
                            else {
                                rec.u = u;
                            }
                            return true;
                        }
                        else if ((sign1 == -1) && (sign2 == 0) && (intercept_x == intercept_z)) {//cone
                            if (((rec.p.y()-center.y()) > -height_half_y) && ((rec.p.y()-center.y()) < 0)) {
                                pc = rec.p - center;
                                rec.normal = unit_vector(vec3(2*bc_square*pc.x(), 2*ac_square*pc.y(), 2*ab_square*pc.z()));
                                if (dot(rec.normal, r.direction()) > 0) {
                                    rec.normal = -rec.normal;
                                }
                                rec.mat_ptr = ma;

                                rec.v = (rec.p.y() - (center.y() - height_half_y)) / (height_half_y);
//                                float rate = intercept_y / intercept_x;
                                vec3 pc1 = vec3(rec.p.x()-center.x(), 0, rec.p.z()-center.z());
                                vec3 vx = vec3(0, 0, 1);
                                float u = acos(dot(pc1, vx) / (pc1.length()*vx.length())) / (2*M_PI);
        //                        float u = acos((rec.p.x() - center.x()) / ((center.y() - rec.p.y()) / rate)) / (2*M_PI);
                                if ((rec.p.x() - center.x()) < 0) {
                                    rec.u = 1-u;
                                }
                                else {
                                    rec.u = u;
                                }
                                return true;
                            }
                            else {
                                rec.u = -1.0;
                                rec.v = -1.0;
                            }
                       }
                        else {
                            rec.u = -1.0;
                            rec.v = -1.0;
                            return true;
                        }
                    }
                    else {
//                        return false;
                    }
                }
                if (temp2 < t_max && temp2 > t_min) {
                    rec.t = temp2;
                    rec.p = r.point_at_parameter(rec.t);
                    if (((rec.p.y()-center.y()) > -height_half_y) && ((rec.p.y()-center.y()) < height_half_y)) {
                        pc = rec.p - center;
                        rec.normal = unit_vector(vec3(2*bc_square*pc.x(), 2*ac_square*pc.y(), 2*ab_square*pc.z()));
                        if (dot(rec.normal, r.direction()) > 0) {
                            rec.normal = -rec.normal;
                        }
                        rec.mat_ptr = ma;
                        if ((sign1 == 0) && (sign2 == 1) && (intercept_x == intercept_y) && (intercept_y == intercept_z)) {//cylinder
                            rec.v = (rec.p.y() - (center.y() - height_half_y)) / (2*height_half_y);
                            vec3 pc1 = vec3(rec.p.x()-center.x(), 0, rec.p.z()-center.z());
                            vec3 vx = vec3(0, 0, 1);
                            float u = acos(dot(pc1, vx) / (pc1.length()*vx.length())) / (2*M_PI);
//                            float u = acos((rec.p.x() - center.x()) / intercept_x) / (2*M_PI);
                            if ((rec.p.x() - center.x()) < 0) {
                                rec.u = 1-u;
                            }
                            else {
                                rec.u = u;
                            }
                            return true;
                        }
                        else if ((sign1 == -1) && (sign2 == 0) && (intercept_x == intercept_z)) {//cone
                            if (((rec.p.y()-center.y()) > -height_half_y) && ((rec.p.y()-center.y()) < 0)) {
                                pc = rec.p - center;
                                rec.normal = unit_vector(vec3(2*bc_square*pc.x(), 2*ac_square*pc.y(), 2*ab_square*pc.z()));
                                if (dot(rec.normal, r.direction()) > 0) {
                                    rec.normal = -rec.normal;
                                }
                                rec.mat_ptr = ma;

                                rec.v = (rec.p.y() - (center.y() - height_half_y)) / (height_half_y);
//                                float rate = intercept_y / intercept_x;
                                vec3 pc1 = vec3(rec.p.x()-center.x(), 0, rec.p.z()-center.z());
                                vec3 vx = vec3(0, 0, 1);
                                float u = acos(dot(pc1, vx) / (pc1.length()*vx.length())) / (2*M_PI);
        //                        float u = acos((rec.p.x() - center.x()) / ((center.y() - rec.p.y()) / rate)) / (2*M_PI);
                                if ((rec.p.x() - center.x()) < 0) {
                                    rec.u = 1-u;
                                }
                                else {
                                    rec.u = u;
                                }
                                return true;
                            }
                            else {
                                rec.u = -1.0;
                                rec.v = -1.0;
                            }
                       }
                        else {
                            rec.u = -1.0;
                            rec.v = -1.0;
                            return true;
                        }
                    }
                    else {
//                        return false;
                    }
                }
            }
            return false;
        }
    return false;
}
