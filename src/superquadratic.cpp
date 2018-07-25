#include "superquadratic.h"

#include <iostream>
#include <limits>
#include "float.h"
#include "log.h"

using namespace std;

bool ray_hit_box_q(const ray& r, const vec3& vertex_l, const vec3& vertex_h, float& t_near, float& t_far) {
        t_near = (numeric_limits<float>::min)();
        t_far = (numeric_limits<float>::max)();
        vec3 direction = r.direction();
        vec3 origin = r.origin();
        vec3 bl = vertex_l;
        vec3 bh = vertex_h;
        float array1[6];

        if(direction.x() == 0) {
            if((origin.x() < bl.x()) || (origin.x() > bh.x())) {
#if SUPERQUADRATIC_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin X0 is not between the slabs. return false" <<endl;
#endif // SUPERQUADRATIC_LOG
                return false;
            }
            array1[0] = (numeric_limits<float>::min)();
            array1[1] = (numeric_limits<float>::max)();
        }
        if(direction.y() == 0) {
            if((origin.y() < bl.y()) || (origin.y() > bh.y())) {
#if SUPERQUADRATIC_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin Y0 is not between the slabs. return false" <<endl;
#endif // SUPERQUADRATIC_LOG
                return false;
            }
            array1[2] = (numeric_limits<float>::min)();
            array1[3] = (numeric_limits<float>::max)();
        }
        if(direction.z() == 0) {
            if((origin.z() < bl.z()) || (origin.z() > bh.z())) {
#if SUPERQUADRATIC_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin Z0 is not between the slabs. return false" <<endl;
#endif // SUPERQUADRATIC_LOG
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
#if SUPERQUADRATIC_LOG == 1
            std::cout << "array1[" << i << "]:" << array1[i] <<endl;
            std::cout << "array1[" << i+1 << "]:" << array1[i+1] <<endl;
#endif // SUPERQUADRATIC_LOG
            if(array1[i] >= t_near) {t_near = array1[i];}
            if(array1[i+1] <= t_far) {t_far = array1[i+1];}
            if(t_near > t_far) {
#if SUPERQUADRATIC_LOG == 1
                std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_near > t_far. return false" <<endl;
#endif // SUPERQUADRATIC_LOG
                return false;
            }
            if(t_far < 0) {
#if SUPERQUADRATIC_LOG == 1
                std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_far < 0. return false" <<endl;
#endif // SUPERQUADRATIC_LOG
                return false;
            }
        }
        return true;
}

float get_superellipsoid_function(float a1, float a2, float a3, float er, float es, float et, float xo, float yo, float zo, float xd, float yd, float zd, float t) {
    return (  pow(fabs(xo/a1 + xd*t/a1), er)
            + pow(fabs(yo/a2 + yd*t/a2), es)
            + pow(fabs(zo/a3 + zd*t/a3), et) - 1);
}

bool get_superellipsoid_function_and_derivative_q(float a1, float a2, float a3, float er, float es, float et, float xo, float yo, float zo, float xd, float yd, float zd, double t, double& f, double& fd) {
        double a1_r = double(a1);
        double a2_r = double(a2);
        double a3_r = double(a3);
        double er_r = double(er);
        double es_r = double(es);
        double et_r = double(et);
        double xo_r = double(xo);
        double yo_r = double(yo);
        double zo_r = double(zo);
        double xd_r = double(xd);
        double yd_r = double(yd);
        double zd_r = double(zd);
        double pow_x, pow_y, pow_z, pow_x_d, pow_y_d, pow_z_d;
        double xd_a1, yd_a2, zd_a3;
        if ((xo_r+xd_r*t) < 0) { xd_a1 = -xd_r/a1_r; }
        else { xd_a1 = xd_r/a1_r; }
        if ((yo_r+yd_r*t) < 0) { yd_a2 = -yd_r/a2_r; }
        else { yd_a2 = yd_r/a2_r; }
        if ((zo_r+zd_r*t) < 0) { zd_a3 = -zd_r/a3_r; }
        else { zd_a3 = zd_r/a3_r; }
        if ((xo_r+xd_r*t) == 0) {
            pow_x = 0;
            pow_x_d = 0;
        }
        else {
            pow_x = pow(fabs(xo_r/a1_r + xd_r*t/a1_r), er_r);
            pow_x_d = pow(fabs(xo_r/a1_r + xd_r*t/a1_r), er_r-1);
        }
        if ((yo_r+yd_r*t) == 0) {
            pow_y = 0;
            pow_y_d = 0;
        }
        else {
            pow_y = pow(fabs(yo_r/a2_r + yd_r*t/a2_r), es_r);
            pow_y_d = pow(fabs(yo_r/a2_r + yd_r*t/a2_r), es_r-1);
        }
        if ((zo_r+zd_r*t) == 0) {
            pow_z = 0;
            pow_x_d = 0;
        }
        else {
            pow_z = pow(fabs(zo_r/a3_r + zd_r*t/a3_r), et_r);
            pow_z_d = pow(fabs(zo_r/a3_r + zd_r*t/a3_r), et_r-1);
        }
        f = pow_x + pow_y + pow_z - 1;
        fd = er_r*pow_x_d*xd_a1 + es_r*pow_y_d*yd_a2 + et_r*pow_z_d*zd_a3;

        if (isnan(f)) {
            f = f * 1;
        }
        if (isnan(fd)) {
            fd = fd * 1;
        }

        return true;
}

bool get_roots_by_newton_iteration_q(const ray& r, vec3 c, float a1, float a2, float a3, float er, float es, float et, float x0[8], float (&roots)[2]) {
        float xo = r.origin().x() - c.x();
        float yo = r.origin().y() - c.y();
        float zo = r.origin().z() - c.z();
        float xd = r.direction().x();
        float yd = r.direction().y();
        float zd = r.direction().z();
        double t_k, t_k1, ft_k, ft_d_k;
        int j=0;

        for (int i=0; i<4; i++) {
            t_k = double(x0[i]);
            for (int k=0; k<50; k++) {
                if (!(isnan(t_k))) {
                    get_superellipsoid_function_and_derivative_q(a1, a2, a3, er, es, et, xo, yo, zo, xd, yd, zd, t_k, ft_k, ft_d_k);
                    if ((ft_d_k != 0) && !(isnan(ft_k)) && !(isnan(ft_d_k))) {
                        t_k1 = t_k - ft_k/ft_d_k;
                        if (fabs(t_k1) >= 1) {
                            if (fabs((t_k1 - t_k)/t_k1) < 0.001) {
                                roots[j+1] = float(t_k1);
                                j++;
                                break;
                            }
                            else {
                                t_k = t_k1;
                            }
                        }
                        else {
                            if (fabs(t_k1 - t_k) < 0.001) {
                                roots[j+1] = float(t_k1);
                                j++;
                                break;
                            }
                            else {
                                t_k = t_k1;
                            }
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
        if (j == 0) {

        }
        return true;
}

bool superquadratic::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if SUPERQUADRATIC_LOG == 1
        std::cout << "-------------superquadratic::hit----------------" << endl;
#endif // SUPERQUADRATIC_LOG
        vec3 vertex_l[4], vertex_h[4];
        vertex_l[0] = vec3(center.x()-intercept_x, center.y()-intercept_y, center.z()-intercept_z);
        vertex_h[0] = vec3(center.x()+intercept_x, center.y()+intercept_y, center.z()+intercept_z);
        float roots[2], x0[8];
        float t_near = 0;
        float t_far = 0;
        for (int i=0; i<1; i++) {
            if (ray_hit_box_q(r, vertex_l[i], vertex_h[i], t_near, t_far)) {
                x0[0] = t_near;
                x0[1] = t_far;
                x0[2] = t_near+1*(t_far - t_near)/4;
                x0[3] = t_near+3*(t_far - t_near)/4;
                x0[4] = t_near+1*(t_far - t_near)/8;
                x0[5] = t_near+7*(t_far - t_near)/8;
                x0[6] = t_near+3*(t_far - t_near)/8;
                x0[7] = t_near+5*(t_far - t_near)/8;

                get_roots_by_newton_iteration_q(r, center, intercept_x, intercept_y, intercept_z, p_r, p_s, p_t, x0, roots);
                if (roots[0] > 0.0001) {
                    break;
                }

            }
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
            double nx, ny, nz;
            float d_a1 = intercept_x;
            float d_a2 = intercept_y;
            float d_a3 = intercept_z;
            vec3 pc;
            for (int k=1; k<int(roots[0])+1; k++) {
                if (roots[k] < t_max && roots[k] > t_min) {
                    rec.t = roots[k];
                    rec.p = r.point_at_parameter(rec.t);
                    pc = rec.p - center;
                    if (pc.x() < 0) {d_a1 = -intercept_x;}
                    if (pc.y() < 0) {d_a2 = -intercept_y;}
                    if (pc.z() < 0) {d_a3 = -intercept_z;}
                    if (pc.x() == 0) {
                        nx = 0;
                    }
                    else {
                        nx = double(p_r)*pow(double(fabs(pc.x()/d_a1)), double(p_r-1))/double(d_a1);
                    }
                    if (pc.y() == 0) {
                        ny = 0;
                    }
                    else {
                        ny = double(p_s)*pow(double(fabs(pc.y()/d_a2)), double(p_s-1))/double(d_a2);
                    }
                    if (pc.z() == 0) {
                        nz = 0;
                    }
                    else {
                        nz = double(p_t)*pow(double(fabs(pc.z()/d_a3)), double(p_t-1))/double(d_a3);
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
