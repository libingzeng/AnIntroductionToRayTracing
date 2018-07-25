#include "quartic_blend_cylinder.h"

#include <iostream>
using namespace std;

bool quartic_blend_cylinder::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if QUARTIC_BLEND_CYLINDER_LOG == 1
        std::cout << "-------------quartic_blend_cylinder::hit----------------" << endl;
#endif // QUARTIC_BLEND_CYLINDER_LOG
        float a1_square = intercept_x1*intercept_x1;
        float b1_square = intercept_z1*intercept_z1;
        float a1b1_square_a3 = a1_square*b1_square+intercept_s1;
        float a2_square = intercept_y2*intercept_y2;
        float b2_square = intercept_z2*intercept_z2;
        float a2b2_square_b3 = a2_square*b2_square+intercept_s2;
        float a3_square = intercept_s1*intercept_s1;
        float b3_square = intercept_s2*intercept_s2;
        float a3b3_sqare = a3_square*b3_square;
        float xd = r.direction().x();
        float yd = r.direction().y();
        float zd = r.direction().z();
        float xoc1 = r.origin().x() - center1.x();
        float zoc1 = r.origin().z() - center1.z();
        float yoc2 = r.origin().y() - center2.y();
        float zoc2 = r.origin().z() - center2.z();

        float A1 = b1_square*xd*xd + a1_square*zd*zd;
        float B1 = 2*b1_square*xd*xoc1 + 2*a1_square*zd*zoc1;
        float C1 = b1_square*xoc1*xoc1 + a1_square*zoc1*zoc1;

        float A2 = A1*A1*b3_square;
        float B2 = 2*A1*B1*b3_square;
        float C2 = (B1*B1 + 2*A1*C1 - 2*a1b1_square_a3*A1)*b3_square;
        float D2 = (2*B1*C1 - 2*a1b1_square_a3*B1)*b3_square;
        float E2 = (C1*C1 - 2*a1b1_square_a3*C1 + a1b1_square_a3*a1b1_square_a3)*b3_square;

        float A3 = b2_square*yd*yd + a2_square*zd*zd;
        float B3 = 2*b2_square*yd*yoc2 + 2*a2_square*zd*zoc2;
        float C3 = b2_square*yoc2*yoc2 + a2_square*zoc2*zoc2;

        float A4 = A3*A3*a3_square;
        float B4 = 2*A3*B3*a3_square;
        float C4 = (B3*B3 + 2*A3*C3 - 2*a2b2_square_b3*A3)*a3_square;
        float D4 = (2*B3*C3 - 2*a2b2_square_b3*B3)*a3_square;
        float E4 = (C3*C3 - 2*a2b2_square_b3*C3 + a2b2_square_b3*a2b2_square_b3)*a3_square - a3b3_sqare;


        float roots[5];
        roots_quartic_equation2(A2+A4, B2+B4, C2+C4, D2+D4, E2+E4, roots);

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
                    vec3 pc1 = rec.p - center1;
                    vec3 pc2 = rec.p - center2;
                    if (((pc1.y()>= -length_half_y1) && (pc1.y() <= length_half_y1)) && ((pc2.x() >= -length_half_x2) && (pc2.x() <= length_half_x2))) {
                        float nx = 2*b3_square*(b1_square*pc1.x()*pc1.x() + a1_square*pc1.z()*pc1.z() - a1b1_square_a3)*2*b1_square*pc1.x();
                        float ny = 2*a3_square*(b2_square*pc2.y()*pc2.y() + a2_square*pc2.z()*pc2.z() - a2b2_square_b3)*2*b2_square*pc2.y();
                        float nz = 2*b3_square*(b1_square*pc1.x()*pc1.x() + a1_square*pc1.z()*pc1.z() - a1b1_square_a3)*2*a1_square*pc1.z()
                                    + 2*a3_square*(b2_square*pc2.y()*pc2.y() + a2_square*pc2.z()*pc2.z() - a2b2_square_b3)*2*a2_square*pc2.z();
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
        }
        return false;
}
