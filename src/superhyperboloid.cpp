#include "superhyperboloid.h"

#include <iostream>
#include <limits>
#include "float.h"
#include "log.h"

using namespace std;

bool ray_hit_box_h(const ray& r, const vec3& vertex_l, const vec3& vertex_h, float& t_near, float& t_far) {
        t_near = (numeric_limits<float>::min)();
        t_far = (numeric_limits<float>::max)();
        vec3 direction = r.direction();
        vec3 origin = r.origin();
        vec3 bl = vertex_l;
        vec3 bh = vertex_h;
        float array1[6];

        if(direction.x() == 0) {
            if((origin.x() < bl.x()) || (origin.x() > bh.x())) {
#if SUPERHYPERBOLOID_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin X0 is not between the slabs. return false" <<endl;
#endif // SUPERHYPERBOLOID_LOG
                return false;
            }
            array1[0] = (numeric_limits<float>::min)();
            array1[1] = (numeric_limits<float>::max)();
        }
        if(direction.y() == 0) {
            if((origin.y() < bl.y()) || (origin.y() > bh.y())) {
#if SUPERHYPERBOLOID_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin Y0 is not between the slabs. return false" <<endl;
#endif // SUPERHYPERBOLOID_LOG
                return false;
            }
            array1[2] = (numeric_limits<float>::min)();
            array1[3] = (numeric_limits<float>::max)();
        }
        if(direction.z() == 0) {
            if((origin.z() < bl.z()) || (origin.z() > bh.z())) {
#if SUPERHYPERBOLOID_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin Z0 is not between the slabs. return false" <<endl;
#endif // SUPERHYPERBOLOID_LOG
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
#if SUPERHYPERBOLOID_LOG == 1
            std::cout << "array1[" << i << "]:" << array1[i] <<endl;
            std::cout << "array1[" << i+1 << "]:" << array1[i+1] <<endl;
#endif // SUPERHYPERBOLOID_LOG
            if(array1[i] >= t_near) {t_near = array1[i];}
            if(array1[i+1] <= t_far) {t_far = array1[i+1];}
            if(t_near > t_far) {
#if SUPERHYPERBOLOID_LOG == 1
                std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_near > t_far. return false" <<endl;
#endif // SUPERHYPERBOLOID_LOG
                return false;
            }
            if(t_far < 0) {
#if SUPERHYPERBOLOID_LOG == 1
                std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_far < 0. return false" <<endl;
#endif // SUPERHYPERBOLOID_LOG
                return false;
            }
        }
        if (t_near != t_near) {
            t_near = t_near * 1;
        }
        return true;
}

bool get_superellipsoid_function_and_derivative_h(float a1, float a2, float a3, float e1, float e2, float s1, float s2, float xo, float yo, float zo, float xd, float yd, float zd, double t, double& f, double& fd) {
        double a1_r = double(a1);
        double a2_r = double(a2);
        double a3_r = double(a3);
        double e1_r = double(e1);
        double e2_r = double(e2);
        double xo_r = double(xo);
        double yo_r = double(yo);
        double zo_r = double(zo);
        double xd_r = double(xd);
        double yd_r = double(yd);
        double zd_r = double(zd);
        double pow_x, pow_y, pow_z, pow_x_d, pow_y_d, pow_z_d;
        double xd_a1, yd_a2, zd_a3;

        if ((xo_r+xd_r*t) < 0) {
            xd_a1 = -xd_r/a1_r;
        }
        else {
            xd_a1 = xd_r/a1_r;
        }
        if ((yo_r+yd_r*t) < 0) {
            yd_a2 = -yd_r/a2_r;
        }
        else {
            yd_a2 = yd_r/a2_r;
        }
        if ((zo_r+zd_r*t) < 0) {
            zd_a3 = -zd_r/a3_r;
        }
        else {
            zd_a3 = zd_r/a3_r;
        }

        if ((xo_r+xd_r*t) == 0) {
            pow_x = 0;
            pow_x_d = 0;
        }
        else {
            pow_x = pow(fabs(xo_r/a1_r + xd_r*t/a1_r), (2/e2_r));
            pow_x_d = pow(fabs(xo_r/a1_r + xd_r*t/a1_r), ((2/e2_r)-1));
        }
        if ((yo_r+yd_r*t) == 0) {
            pow_y = 0;
            pow_y_d = 0;
        }
        else {
            pow_y = pow(fabs(yo_r/a2_r + yd_r*t/a2_r), (2/e1_r));
            pow_y_d = pow(fabs(yo_r/a2_r + yd_r*t/a2_r), ((2/e1_r)-1));
        }
        if ((zo_r+zd_r*t) == 0) {
            pow_z = 0;
            pow_z_d = 0;
        }
        else {
            pow_z = pow(fabs(zo_r/a3_r + zd_r*t/a3_r), (2/e2_r));
            pow_z_d = pow(fabs(zo_r/a3_r + zd_r*t/a3_r), ((2/e2_r)-1));
        }

        if ((pow_x+s2*pow_z) == 0) {
            f = pow_y - 1;
            fd = (2/e1_r)*pow_y_d*yd_a2;
        }
        else if ((pow_x+s2*pow_z) < 0)  {
            f = -pow(-(pow_x+s2*pow_z), (e2_r/e1_r)) + s1*pow_y - 1;
            fd = (2/e1_r)*(-pow(-(pow_x+s2*pow_z), ((e2_r/e1_r)-1))*(pow_x_d*xd_a1+s2*pow_z_d*zd_a3) + s1*pow_y_d*yd_a2);
        }
        else {
            f = pow(pow_x+s2*pow_z, (e2_r/e1_r)) + s1*pow_y - 1;
            fd = (2/e1_r)*(pow(pow_x+s2*pow_z, ((e2_r/e1_r)-1))*(pow_x_d*xd_a1+s2*pow_z_d*zd_a3) + s1*pow_y_d*yd_a2);
        }
/*
        if (f != f) {
            f = f * 1;
        }
        if (fd != fd) {
            fd = fd * 1;
        }
*/
        return true;
}

bool get_roots_by_newton_iteration_h(const ray& r, vec3 c, float a1, float a2, float a3, float e1, float e2, float s1, float s2, int in, float tol, float *x0, float (&roots)[3]) {
        float xo = r.origin().x() - c.x();
        float yo = r.origin().y() - c.y();
        float zo = r.origin().z() - c.z();
        float xd = r.direction().x();
        float yd = r.direction().y();
        float zd = r.direction().z();
        double t_k, t_k1, ft_k, ft_d_k;
        int j=0, in_r;
        if (in > int(x0[0])) {
            in_r = int(x0[0]);
        }
        else {
            in_r = in;
        }

        for (int i=1; i<=in_r; i++) {
            t_k = double(x0[i]);
            for (int k=0; k<50; k++) {
                if (!(isnan(t_k))) {
                    get_superellipsoid_function_and_derivative_h(a1, a2, a3, e1, e2, s1, s2, xo, yo, zo, xd, yd, zd, t_k, ft_k, ft_d_k);
                    if ((ft_d_k != 0) && !(isnan(ft_k)) && !(isnan(ft_d_k))) {
                        t_k1 = t_k - ft_k/ft_d_k;
//                        if (fabs(t_k1) >= 1) {
                            if (fabs((t_k1 - t_k)/t_k1) < tol) {
                                if ((t_k1 >= x0[1]) && (t_k1 <= x0[in_r])) {
                                    roots[j+1] = float(t_k1);
                                    j++;
                                    break;
                                }
                                else {
                                    break;
                                }
                            }
                            else {
                                t_k = t_k1;
                            }
                    }
                    else {
                        break;
                    }
                }
                else {
                    break;
                }
            }

            if (j == 1) {
                break;
            }

        }
        roots[0] = float(j);
        return true;
}

bool superhyperboloid::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if SUPERHYPERBOLOID_LOG == 1
        std::cout << "-------------superhyperboloid::hit----------------" << endl;
#endif // SUPERHYPERBOLOID_LOG
        vec3 vertex_l[1], vertex_h[1];
        float intercept_hyper_x, intercept_hyper_z;
        float pow_y = pow((half_y/intercept_y), (2/p_e1));
        intercept_hyper_x = intercept_x*pow((pow_y+1), (p_e1/2));
        intercept_hyper_z = intercept_z*pow((pow_y+1), (p_e1/2));
        vertex_l[0] = vec3(center.x()-intercept_hyper_x, center.y()-half_y, center.z()-intercept_hyper_z);
        vertex_h[0] = vec3(center.x()+intercept_hyper_x, center.y()+half_y, center.z()+intercept_hyper_z);
        float roots[3] = {0.0, -1.0, -1.0};
        float x0[initial_number+1];
        float t_near = 0;
        float t_far = 0;
        if (ray_hit_box_h(r, vertex_l[0], vertex_h[0], t_near, t_far)) {
            for (int i=0; i<initial_number; i++) {
                x0[i+1] = t_near + i*(t_far - t_near)/(initial_number-1);
            }
            x0[0] = float(initial_number);
            get_roots_by_newton_iteration_h(r, center, intercept_x, intercept_y, intercept_z, p_e1, p_e2, sign1, sign2, initial_number, tolerance, x0, roots);
        }
        else {
            return false;
        }

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
            double nx, ny, nz, pow_x, pow_z, pow_x_d, pow_z_d, pow_y_d, pow_x_z_d;
            float d_a1 = intercept_x;
            float d_a2 = intercept_y;
            float d_a3 = intercept_z;
            vec3 pc;
            for (int k=1; k<int(roots[0])+1; k++) {
                if (roots[k] < t_max && roots[k] > t_min) {
                    rec.t = roots[k];
                    rec.p = r.point_at_parameter(rec.t);
                    if (((sign1 == -1)&&(sign2 == 1)) || ((sign1 == -1)&&(sign2 == -1))) {
                        if (!(((rec.p.y()-center.y()) >= -half_y) && ((rec.p.y()-center.y()) <= half_y))) {
                            continue;
                        }
                    }
                    pc = rec.p - center;
                    if (pc.x() < 0) {d_a1 = -intercept_x;}
                    if (pc.y() < 0) {d_a2 = -intercept_y;}
                    if (pc.z() < 0) {d_a3 = -intercept_z;}
                    if (pc.x() == 0) {
                        pow_x = 0;
                        pow_x_d = 0;
                    }
                    else {
                        pow_x = pow(double(fabs(pc.x()/d_a1)), double(2/p_e2));
                        pow_x_d = pow(double(fabs(pc.x()/d_a1)), double(2/p_e2-1));
                    }
                    if (pc.y() == 0) {
                        pow_y_d = 0;
                    }
                    else {
                        pow_y_d = pow(double(fabs(pc.y()/d_a2)), double(2/p_e1-1));
                    }
                    if (pc.z() == 0) {
                        pow_z = 0;
                        pow_z_d = 0;
                    }
                    else {
                        pow_z = pow(double(fabs(pc.z()/d_a3)), double(2/p_e2));
                        pow_z_d = pow(double(fabs(pc.z()/d_a3)), double(2/p_e2-1));
                    }
                    if ((pow_x+sign2*pow_z) == 0) {
                        pow_x_z_d = 0;
                    }
                    if ((pow_x+sign2*pow_z) < 0) {
                        pow_x_z_d = -pow(-(pow_x+sign2*pow_z), double(p_e2/p_e1-1));
                    }
                    else {
                        pow_x_z_d = pow(pow_x+sign2*pow_z, double(p_e2/p_e1-1));
                    }

                    nx = double(2/p_e1) * pow_x_z_d  * pow_x_d / d_a1;
                    ny = double(2/p_e1) *sign1* pow_y_d / d_a2;
                    nz = double(2/p_e1) * pow_x_z_d  * pow_z_d / d_a3;

                    if (isnan(nx)) {
                        nx = nx * 1;
                    }

                    nx = nx/sqrt(nx*nx+ny*ny+nz*nz);
                    ny = ny/sqrt(nx*nx+ny*ny+nz*nz);
                    nz = nz/sqrt(nx*nx+ny*ny+nz*nz);

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
            return false;
        }
        else {
            return false;
        }
        return false;
}
