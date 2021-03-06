#include "tori_part.h"

#include <iostream>
using namespace std;

bool tori_part::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if TORI_PART_LOG == 1
        std::cout << "-------------tori_part::hit---1-------------" << endl;
#endif // TORI_PART_LOG
        vec3 oc = r.origin() - center;
        float A = dot(oc, oc);
        float B = 2*dot(oc, r.direction());
        float C = dot(r.direction(), r.direction());
        float Rr_square_p = radius_a*radius_a +radius_b*radius_b;
        float Rr_square_s_square = (radius_a*radius_a - radius_b*radius_b) * (radius_a*radius_a - radius_b*radius_b);
        float R_square = radius_a*radius_a;
        float a4 = C*C;
        float a3 = 2*B*C;
        float a2 = B*B + 2*A*C - 2*Rr_square_p*C + 4*R_square*r.direction().z()*r.direction().z();
        float a1 = 2*A*B - 2*Rr_square_p*B + 8*R_square*oc.z()*r.direction().z();
        float a0 = A*A - 2*Rr_square_p*A + 4*R_square*oc.z()*oc.z() + Rr_square_s_square;
/*
        float a4 = dot(r.direction(), r.direction()) * dot(r.direction(), r.direction());
        float a3 = 4*(dot(oc, r.direction())) * dot(r.direction(), r.direction());
        float a2 = 2*dot(r.direction(), r.direction()) * (dot(oc, oc) - (radius_a*radius_a+radius_b*radius_b))
                    + 4*(dot(oc, r.direction()))*(dot(oc, r.direction()))
                    + 4*radius_a*radius_a*r.direction().z()*r.direction().z();
        float a1 = 4*dot(oc, r.direction())*(dot(oc, oc) - (radius_a*radius_a+radius_b*radius_b))
                    + 8*radius_a*radius_a*oc.z()*r.direction().z();
        float a0 = (dot(oc, oc) - (radius_a*radius_a+radius_b*radius_b))*(dot(oc, oc) - (radius_a*radius_a+radius_b*radius_b))
                    - 4*radius_a*radius_a*(radius_b*radius_b-oc.z()*oc.z());
*/
//        float *roots = roots_quartic_equation(a4, a3, a2, a1, a0);

        float roots[5];
        roots_quartic_equation2(a4, a3, a2, a1, a0, roots);

        float temp;
        if (roots[0] > 0.0001) {
            for (int i=1; i<int(roots[0]); i++) {
                for (int j=i+1; j<int(roots[0])+1; j++) {
                    if (roots[i] > roots[j]) {
                        temp = roots[i];
                        roots[i] = roots[j];
                        roots[j] = temp;
                    }
                }
            }
            for (int k=1; k<int(roots[0])+1; k++) {
                if (roots[k] < t_max && roots[k] > t_min) {
                    rec.t = roots[k];
                    rec.p = r.point_at_parameter(rec.t);
                    vec3 pc = rec.p - center;
                    vec3 x = vec3(1.0, 0.0, 0.0);

                    float cos_theta = dot(pc, x) / (pc.length() * x.length());
                    float theta = acos(cos_theta);
                    if (rec.p.y() < center.y()) {
                        theta = 2*M_PI - theta;
                    }
                    if (((theta >= theta1) && (theta <= theta2)) || ((theta + 2*M_PI) <= theta2)) {
                        float nx = 4*pc.x()*pc.x()*pc.x() + 4*pc.x()*(pc.y()*pc.y()+pc.z()*pc.z()-(radius_a*radius_a+radius_b*radius_b));
                        float ny = 4*pc.y()*pc.y()*pc.y() + 4*pc.y()*(pc.x()*pc.x()+pc.z()*pc.z()-(radius_a*radius_a+radius_b*radius_b));
                        float nz = 4*pc.z()*pc.z()*pc.z() + 4*pc.z()*(pc.x()*pc.x()+pc.y()*pc.y()+radius_a*radius_a-radius_b*radius_b);
                        rec.normal = unit_vector(vec3(nx, ny, nz));
                        if(dot(r.direction(), rec.normal) > 0) {
                            rec.normal = - rec.normal;
                        }
                        rec.mat_ptr = ma;
                        rec.u = -1.0;
                        rec.v = -1.0;
//xcode                       delete [] roots;
                        return true;
                    }
                }
            }
        }
//xcode        delete [] roots;
        return false;
}
