#include "quadratic_cylinder_all.h"

#include <iostream>
using namespace std;

bool quadratic_cylinder_all::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if QUADRATIC_CYLINDER_ALL_LOG == 1
        std::cout << "-------------quadratic_cylinder_all::hit----------------" << endl;
#endif // QUADRATIC_CYLINDER_ALL_LOG
        vec3 direction = vector_trans(r.direction(), vector_v, vector_u, vector_w);
        vec3 origin = vector_trans(r.origin(), vector_v, vector_u, vector_w);
        vec3 center_vuw = vector_trans(center, vector_v, vector_u, vector_w);

        float ab_square = intercept_x*intercept_x*intercept_y*intercept_y;
        float bc_square = intercept_y*intercept_y*intercept_z*intercept_z;
        float ac_square = 0.0*intercept_x*intercept_x*intercept_z*intercept_z;
        float abc_square = 1.0*intercept_x*intercept_x*intercept_y*intercept_y*intercept_z*intercept_z;

        vec3 inter_square = vec3(bc_square, ac_square, ab_square);
        vec3 rd_square = vec3(direction.x()*direction.x(),
                              direction.y()*direction.y(),
                              direction.z()*direction.z());
        float A = dot(inter_square, rd_square);
        vec3 r0_c = origin - center_vuw;
        vec3 r0_c_rd = vec3(r0_c.x()*direction.x(),
                            r0_c.y()*direction.y(),
                            r0_c.z()*direction.z());
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
                    rec.p = vector_trans(r.point_at_parameter(rec.t), vector_v, vector_u, vector_w);
                    if (((rec.p.y()-center_vuw.y()) > -height_half_y) && ((rec.p.y()-center_vuw.y()) < height_half_y)) {
                        pc = rec.p - center_vuw;
                        rec.normal = unit_vector(vec3(2*bc_square*pc.x(), 2*ac_square*pc.y(), 2*ab_square*pc.z()));
                        if (dot(rec.normal, direction) > 0) {
                            rec.normal = -rec.normal;
                        }
                        rec.normal = vector_trans_back(rec.normal, vector_v, vector_u, vector_w);
                        rec.mat_ptr = ma;
                        if ((intercept_x == intercept_y) && (intercept_y == intercept_z)) {//0/1:cylinder
                            rec.v = (rec.p.y() - (center_vuw.y() - height_half_y)) / (2*height_half_y);
                            vec3 pc1 = vec3(rec.p.x()-center_vuw.x(), 0, rec.p.z()-center_vuw.z());
                            vec3 vx = vec3(0, 0, 1);
                            float u = acos(dot(pc1, vx) / (pc1.length()*vx.length())) / (2*M_PI);
//                            float u = acos((rec.p.x() - center_vuw.x()) / intercept_x) / (2*M_PI);
                            if ((rec.p.x() - center_vuw.x()) < 0) {
                                rec.u = 1-u;
                            }
                            else {
                                rec.u = u;
                            }
                            rec.p = vector_trans_back(rec.p, vector_v, vector_u, vector_w);
                            return true;
                        }
                        else {
                            rec.u = -1.0;
                            rec.v = -1.0;
                            rec.p = vector_trans_back(rec.p, vector_v, vector_u, vector_w);
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
                    rec.p = vector_trans(r.point_at_parameter(rec.t), vector_v, vector_u, vector_w);
                    if (((rec.p.y()-center_vuw.y()) > -height_half_y) && ((rec.p.y()-center_vuw.y()) < height_half_y)) {
                        pc = rec.p - center_vuw;
                        rec.normal = unit_vector(vec3(2*bc_square*pc.x(), 2*ac_square*pc.y(), 2*ab_square*pc.z()));
                        if (dot(rec.normal, direction) > 0) {
                            rec.normal = -rec.normal;
                        }
                        rec.normal = vector_trans_back(rec.normal, vector_v, vector_u, vector_w);
                        rec.mat_ptr = ma;
                        if ((intercept_x == intercept_y) && (intercept_y == intercept_z)) {//cylinder
                            rec.v = (rec.p.y() - (center_vuw.y() - height_half_y)) / (2*height_half_y);
                            vec3 pc1 = vec3(rec.p.x()-center_vuw.x(), 0, rec.p.z()-center_vuw.z());
                            vec3 vx = vec3(0, 0, 1);
                            float u = acos(dot(pc1, vx) / (pc1.length()*vx.length())) / (2*M_PI);
//                            float u = acos((rec.p.x() - center_vuw.x()) / intercept_x) / (2*M_PI);
                            if ((rec.p.x() - center_vuw.x()) < 0) {
                                rec.u = 1-u;
                            }
                            else {
                                rec.u = u;
                            }
                            rec.p = vector_trans_back(rec.p, vector_v, vector_u, vector_w);
                            return true;
                        }
                        else {
                            rec.u = -1.0;
                            rec.v = -1.0;
                            rec.p = vector_trans_back(rec.p, vector_v, vector_u, vector_w);
                            return true;
                        }
                    }
                    else {
//                        return false;
                    }
                }
                if (temp2 < t_max && temp2 > t_min) {
                    rec.t = temp2;
                    rec.p = vector_trans(r.point_at_parameter(rec.t), vector_v, vector_u, vector_w);
                    if (((rec.p.y()-center_vuw.y()) > -height_half_y) && ((rec.p.y()-center_vuw.y()) < height_half_y)) {
                        pc = rec.p - center_vuw;
                        rec.normal = unit_vector(vec3(2*bc_square*pc.x(), 2*ac_square*pc.y(), 2*ab_square*pc.z()));
                        if (dot(rec.normal, direction) > 0) {
                            rec.normal = -rec.normal;
                        }
                        rec.normal = vector_trans_back(rec.normal, vector_v, vector_u, vector_w);
                        rec.mat_ptr = ma;
                        if ((intercept_x == intercept_y) && (intercept_y == intercept_z)) {//cylinder
                            rec.v = (rec.p.y() - (center_vuw.y() - height_half_y)) / (2*height_half_y);
                            vec3 pc1 = vec3(rec.p.x()-center_vuw.x(), 0, rec.p.z()-center_vuw.z());
                            vec3 vx = vec3(0, 0, 1);
                            float u = acos(dot(pc1, vx) / (pc1.length()*vx.length())) / (2*M_PI);
//                            float u = acos((rec.p.x() - center_vuw.x()) / intercept_x) / (2*M_PI);
                            if ((rec.p.x() - center_vuw.x()) < 0) {
                                rec.u = 1-u;
                            }
                            else {
                                rec.u = u;
                            }
                            rec.p = vector_trans_back(rec.p, vector_v, vector_u, vector_w);
                            return true;
                        }
                        else {
                            rec.u = -1.0;
                            rec.v = -1.0;
                            rec.p = vector_trans_back(rec.p, vector_v, vector_u, vector_w);
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
