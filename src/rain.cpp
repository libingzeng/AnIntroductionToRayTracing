#include "rain.h"
#include <iostream>
using namespace std;

bool rain::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if RAIN_LOG == 1
        std::cout << "-------------rain::hit---1-------------" << endl;
#endif // RAIN_LOG
        vec3 oc1 = r.origin() - center;
        vec3 dir1 = r.direction();
        vec3 oc = vec3(double(oc1.x()), double(oc1.y()), double(oc1.z()));
        vec3 dir = vec3(double(dir1.x()), double(dir1.y()), double(dir1.z()));
        double A = 4*z_a3*z_a3*y_a2*y_a2*y_a2*y_a2;
        double B = 4*x_a1*x_a1*y_a2*y_a2*y_a2*y_a2;
        double C = x_a1*x_a1*z_a3*z_a3;
        double D = -2*x_a1*x_a1*z_a3*z_a3*y_a2;
        double E = 2*x_a1*x_a1*z_a3*z_a3*y_a2*y_a2*y_a2;
        double F = -x_a1*x_a1*z_a3*z_a3*y_a2*y_a2*y_a2*y_a2;
        double a4 = dir.y()*dir.y()*dir.y()*dir.y()*C;
        double a3 = 4*oc.y()*dir.y()*dir.y()*dir.y()*C + dir.y()*dir.y()*dir.y()*D;
        double a2 = 6*oc.y()*oc.y()*dir.y()*dir.y()*C + 3*oc.y()*dir.y()*dir.y()*D + dir.z()*dir.z()*B + dir.x()*dir.x()*A;
        double a1 = 4*oc.y()*oc.y()*oc.y()*dir.y()*C + 3*oc.y()*oc.y()*dir.y()*D + 2*oc.z()*dir.z()*B + 2*oc.x()*dir.x()*A + dir.y()*E;
        double a0 = oc.y()*oc.y()*oc.y()*oc.y()*C + oc.y()*oc.y()*oc.y()*D + oc.y()*E + F + oc.z()*oc.z()*B + oc.x()*oc.x()*A;

        double roots[5];
        roots_quartic_equation2_rain(a4, a3, a2, a1, a0, roots);

        double temp;
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
            vec3 pc;
            double nx, ny, nz;
            for (int k=1; k<int(roots[0])+1; k++) {
                if ((float(roots[k]) < t_max) && (float(roots[k]) > t_min)) {
                    rec.t = float(roots[k]);
                    rec.p = r.point_at_parameter(rec.t);
                    pc = rec.p - center;
                    if (rec.p.y() < 1) {
                    //this is a bad trick, there are something wrong in the progam when r.direction().y() superimposes the drop picture.
                        continue;
                    }
                    nx = 2.0*A*double(pc.x());
                    ny = 4.0*C*double(pc.y())*double(pc.y())*double(pc.y()) + 3.0*D*double(pc.y())*double(pc.y()) + E;
                    nz = 2.0*B*double(pc.z());

                    double length = sqrt(nx*nx+ny*ny+nz*nz);
                    if (fabs(length) >= EPS) {
                    /* if the length is very small,
                    when we trasform the type from double to float,
                    we will get: nx=ny=nz=0.*/
                        nx = nx/length;
                        ny = ny/length;
                        nz = nz/length;
                    }

                    rec.normal = unit_vector(vec3(float(nx), float(ny), float(nz)));
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
