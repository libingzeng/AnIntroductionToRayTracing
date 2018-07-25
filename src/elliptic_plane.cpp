#include "elliptic_plane.h"

#include <iostream>
using namespace std;

bool elliptic_plane::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if ELLIPTIC_PLANE_LOG == 1
        std::cout << "-------------elliptic_plane::hit----------------" << endl;
#endif // ELLIPTIC_PLANE_LOG
        if(dot(normal, r.direction()) == 0) {
            return false;
        }
        else {
            if (dot(normal, r.direction()) > 0) {
                rec.normal = - normal;
            }
            float temp = (dot(normal, center) - dot(normal, r.origin())) / dot(normal, r.direction());
            if (temp < t_max && temp > t_min) {
                rec.t = temp;
                rec.p = r.point_at_parameter(rec.t);
                if (((((rec.p.x()-center.x())*(rec.p.x()-center.x())) / (intercept_x*intercept_x))
                   + (((rec.p.y()-center.y())*(rec.p.y()-center.y())) / (intercept_y*intercept_y))
                   + (((rec.p.z()-center.z())*(rec.p.z()-center.z())) / (intercept_z*intercept_z))) <= 1) {
                    rec.mat_ptr = ma;
                    if ((intercept_x == intercept_y) && (intercept_y == intercept_z)) {
                        rec.v = sqrt(((rec.p.x()-center.x())*(rec.p.x()-center.x()) + (rec.p.y()-center.y())*(rec.p.y()-center.y())) / (intercept_x*intercept_x));
                        vec3 pc1 = vec3(rec.p.x()-center.x(), rec.p.y()-center.y(), 0);
                        vec3 vx = vec3(1, 0, 0);
                        float u = acos(dot(pc1, vx) / (pc1.length()*vx.length())) / (2*M_PI);
//                        float u = acos((rec.p.x()-center.x()) / sqrt((rec.p.x()-center.x())*(rec.p.x()-center.x()) + (rec.p.y()-center.y())*(rec.p.y()-center.y()))) / (2*M_PI);
                        if ((rec.p.y()-center.y()) < 0) {
                            rec.u = 1 - u;
                        }
                        else {
                            rec.u = u;
                        }
                    }
                    else {
                        rec.u = -1.0;
                        rec.v = -1.0;
                    }

                    return true;
                }
            }
            return false;
        }
}
