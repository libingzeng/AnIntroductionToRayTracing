#include "blobs.h"

#include <iostream>
#include <limits>
#include "float.h"
#include "log.h"

using namespace std;

bool ray_hit_box_b(const ray& r, const vec3& vertex_l, const vec3& vertex_h, float& t_near, float& t_far) {
        t_near = (numeric_limits<float>::min)();
        t_far = (numeric_limits<float>::max)();
        vec3 direction = r.direction();
        vec3 origin = r.origin();
        vec3 bl = vertex_l;
        vec3 bh = vertex_h;
        float array1[6];

        if(direction.x() == 0) {
            if((origin.x() < bl.x()) || (origin.x() > bh.x())) {
#if BLOBS_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin X0 is not between the slabs. return false" <<endl;
#endif // BLOBS_LOG
                return false;
            }
            array1[0] = (numeric_limits<float>::min)();
            array1[1] = (numeric_limits<float>::max)();
        }
        if(direction.y() == 0) {
            if((origin.y() < bl.y()) || (origin.y() > bh.y())) {
#if BLOBS_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin Y0 is not between the slabs. return false" <<endl;
#endif // BLOBS_LOG
                return false;
            }
            array1[2] = (numeric_limits<float>::min)();
            array1[3] = (numeric_limits<float>::max)();
        }
        if(direction.z() == 0) {
            if((origin.z() < bl.z()) || (origin.z() > bh.z())) {
#if BLOBS_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin Z0 is not between the slabs. return false" <<endl;
#endif // BLOBS_LOG
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
#if BLOBS_LOG == 1
            std::cout << "array1[" << i << "]:" << array1[i] <<endl;
            std::cout << "array1[" << i+1 << "]:" << array1[i+1] <<endl;
#endif // BLOBS_LOG
            if(array1[i] >= t_near) {t_near = array1[i];}
            if(array1[i+1] <= t_far) {t_far = array1[i+1];}
            if(t_near > t_far) {
#if BLOBS_LOG == 1
                std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_near > t_far. return false" <<endl;
#endif // BLOBS_LOG
                return false;
            }
            if(t_far < 0) {
#if BLOBS_LOG == 1
                std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_far < 0. return false" <<endl;
#endif // BLOBS_LOG
                return false;
            }
        }
        if (t_near != t_near) {
            t_near = t_near * 1;
        }
        return true;
}

bool get_blobs_function_and_derivative_b(const ray& r, vec3 c1, vec3 c2, float b1, float b2, float r1, float r2, float s, double t, double& f, double& fd) {
        double xo1 = double(r.origin().x() - c1.x());
        double yo1 = double(r.origin().y() - c1.y());
        double zo1 = double(r.origin().z() - c1.z());
        double xo2 = double(r.origin().x() - c2.x());
        double yo2 = double(r.origin().y() - c2.y());
        double zo2 = double(r.origin().z() - c2.z());
        double xd = double(r.direction().x());
        double yd = double(r.direction().y());
        double zd = double(r.direction().z());
        double b1_r1 = double(b1/(r1*r1));
        double b2_r2 = double(b2/(r2*r2));
        double xo1_t = xo1+xd*t;
        double yo1_t = yo1+yd*t;
        double zo1_t = zo1+zd*t;
        double xo2_t = xo2+xd*t;
        double yo2_t = yo2+yd*t;
        double zo2_t = zo2+zd*t;
        double r1_2 = (xo1_t*xo1_t+yo1_t*yo1_t+zo1_t*zo1_t);
        double r2_2 = (xo2_t*xo2_t+yo2_t*yo2_t+zo2_t*zo2_t);
        double e1 = exp(b1_r1*(xo1_t*xo1_t+yo1_t*yo1_t+zo1_t*zo1_t)-double(b1));
        double e2 = exp(b2_r2*(xo2_t*xo2_t+yo2_t*yo2_t+zo2_t*zo2_t)-double(b2));

        f = e1 + e2 - double(s);
        fd = e1*b1_r1*2*(xo1_t*xd+yo1_t*yd+zo1_t*zd) + e2*b2_r2*2*(xo2_t*xd+yo2_t*yd+zo2_t*zd);

        if (e1 == 1.0) {
            f = f*1;
        }
        return true;
}

bool get_roots_by_newton_iteration_b(const ray& r, vec3 c1, vec3 c2, float b1, float b2, float r1, float r2, float s, int in, float tol, float *x0, float (&roots)[2]) {
        double t_k, t_k1, ft_k, ft_d_k;
        int j=0, in_r;
        if (in > int(x0[0])) {
            in_r = int(x0[0]);
        }
        else {
            in_r = in;
        }

        for (int i=1; i<in_r; i++) {
            t_k = double(x0[i]);
            for (int k=0; k<50; k++) {
                if (!(isnan(t_k))) {
                    get_blobs_function_and_derivative_b(r, c1, c2, b1, b2, r1, r2, s, t_k, ft_k, ft_d_k);
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
/*
                        }
                        else {
                            if (fabs(t_k1 - t_k) < tol) {
                                roots[j+1] = float(t_k1);
                                j++;
                                break;
                            }
                            else {
                                t_k = t_k1;
                            }
                        }
*/
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

bool blobs::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if BLOBS_LOG == 1
        std::cout << "-------------blobs::hit----------------" << endl;
#endif // BLOBS_LOG
        float box_blobs_x, box_blobs_y, box_blobs_z, center_y;
        box_blobs_x = ((radius1 > radius2)? radius1:radius2);
        box_blobs_y = (fabs(center1.y()-center2.y())+radius1+radius2)/2;
        box_blobs_z = ((radius1 > radius2)? radius1:radius2);
        if (center1.y() > center2.y()) {
            center_y = ((center1.y()+radius1)+(center2.y()-radius2))/2;
        }
        else {
            center_y = ((center2.y()+radius2)+(center1.y()-radius1))/2;
        }
        vec3 vertex_l[1], vertex_h[1];
        vertex_l[0] = vec3(center1.x()-box_blobs_x, center_y-box_blobs_y, center1.z()-box_blobs_z);
        vertex_h[0] = vec3(center1.x()+box_blobs_x, center_y+box_blobs_y, center1.z()+box_blobs_z);


        float roots[2] = {0.0, -1.0};
        float x0[initial_number+1];
        float t_near = 0;
        float t_far = 0;
        if (ray_hit_box_b(r, vertex_l[0], vertex_h[0], t_near, t_far)) {
            if (initial_number == 1) {
                x0[1] = t_near;
            }
            else {
                for (int i=0; i<initial_number; i++) {
                    x0[i+1] = t_near + i*(t_far - t_near)/(initial_number-1);
                }
            }
            x0[0] = float(initial_number);
            get_roots_by_newton_iteration_b(r, center1, center2, blob_p1, blob_p2, radius1, radius2, sum, initial_number, tolerance, x0, roots);
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
            vec3 pc1, pc2;
            double b1_r1 = double(blob_p1/(radius1*radius1));
            double b2_r2 = double(blob_p2/(radius2*radius2));
            double a1, a2, nx, ny, nz;
            for (int k=1; k<int(roots[0])+1; k++) {
                if (roots[k] < t_max && roots[k] > t_min) {
                    rec.t = roots[k];
                    rec.p = r.point_at_parameter(rec.t);
                    pc1 = rec.p - center1;
                    pc2 = rec.p - center2;
                    a1 = b1_r1*2*exp(b1_r1*double(dot(pc1, pc1))-double(blob_p1));
                    a2 = b2_r2*2*exp(b2_r2*double(dot(pc2, pc2))-double(blob_p2));

                    nx = a1*double(pc1.x())+a2*double(pc2.x());
                    ny = a1*double(pc1.y())+a2*double(pc2.y());
                    nz = a1*double(pc1.z())+a2*double(pc2.z());

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
