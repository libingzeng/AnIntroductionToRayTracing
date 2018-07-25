#include "quadratic_ellipsoid.h"

#include <iostream>
using namespace std;

bool quadratic_ellipsoid::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if QUADRATIC_ELLIPSOID_LOG == 1
        std::cout << "-------------quadratic_ellipsoid::hit----------------" << endl;
#endif // QUADRATIC_ELLIPSOID_LOG
        float ab_square = intercept_x*intercept_x*intercept_y*intercept_y;
        float bc_square = intercept_y*intercept_y*intercept_z*intercept_z;
        float ac_square = intercept_x*intercept_x*intercept_z*intercept_z;
        float abc_square = intercept_x*intercept_x*intercept_y*intercept_y*intercept_z*intercept_z;

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

        float discriminant = B*B - 4*A*C;
        if (discriminant >= 0) {
            float temp;
            vec3 pc;
            temp = (-B - sqrt(discriminant)) / (2.0*A);
            if (temp < t_max && temp > t_min) {
                rec.t = temp;
                rec.p = r.point_at_parameter(rec.t);
                pc = rec.p - center;
                rec.normal = unit_vector(vec3(2*bc_square*pc.x(), 2*ac_square*pc.y(), 2*ab_square*pc.z()));
                rec.mat_ptr = ma;
                rec.u = -1.0;
                rec.v = -1.0;
                return true;
            }
            temp = (-B + sqrt(discriminant)) / (2.0*A);
            if (temp < t_max && temp > t_min) {
                rec.t = temp;
                rec.p = r.point_at_parameter(rec.t);
                pc = rec.p - center;
                rec.normal = unit_vector(vec3(2*bc_square*pc.x(), 2*ac_square*pc.y(), 2*ab_square*pc.z()));
                rec.mat_ptr = ma;
                rec.u = -1.0;
                rec.v = -1.0;
                return true;
            }
        }
        return false;
}
