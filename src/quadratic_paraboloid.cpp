#include "quadratic_paraboloid.h"

#include <iostream>
using namespace std;

bool quadratic_paraboloid::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if QUADRATIC_PARABOLOID_LOG == 1
        std::cout << "-------------quadratic_paraboloid::hit----------------" << endl;
#endif // QUADRATIC_PARABOLOID_LOG
        float A = focus_directrix_q*r.direction().x()*r.direction().x()
                + focus_directrix_p*r.direction().z()*r.direction().z();
        float B = 2*focus_directrix_q*r.direction().x()*(r.origin().x()-center.x())
                + 2*focus_directrix_p*r.direction().z()*(r.origin().z()-center.z())
                - 2*focus_directrix_p*focus_directrix_q*r.direction().y();
        float C = focus_directrix_q*(r.origin().x()-center.x())*(r.origin().x()-center.x())
                + focus_directrix_p*(r.origin().z()-center.z())*(r.origin().z()-center.z())
                - 2*focus_directrix_p*focus_directrix_q*(r.origin().y()-center.y());
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
                        rec.normal = unit_vector(vec3(2*focus_directrix_q*pc.x(), -2*focus_directrix_p*focus_directrix_q, 2*focus_directrix_p*pc.z()));
                        if (dot(rec.normal, r.direction()) > 0) {
                            rec.normal = -rec.normal;
                        }
                        rec.mat_ptr = ma;
                        return true;
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
                        rec.normal = unit_vector(vec3(2*focus_directrix_q*pc.x(), -2*focus_directrix_p*focus_directrix_q, 2*focus_directrix_p*pc.z()));
                        if (dot(rec.normal, r.direction()) > 0) {
                            rec.normal = -rec.normal;
                        }
                        rec.mat_ptr = ma;
                        return true;
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
                        rec.normal = unit_vector(vec3(2*focus_directrix_q*pc.x(), -2*focus_directrix_p*focus_directrix_q, 2*focus_directrix_p*pc.z()));
                        if (dot(rec.normal, r.direction()) > 0) {
                            rec.normal = -rec.normal;
                        }
                        rec.mat_ptr = ma;
                        return true;
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
