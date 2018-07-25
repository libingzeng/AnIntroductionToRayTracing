#include "parametric_surface.h"

#include <iostream>
#include <limits>
#include "float.h"
#include "log.h"

using namespace std;

bool ray_hit_box_para(const ray& r, const vec3& vertex_l, const vec3& vertex_h, float& t_near, float& t_far) {
        t_near = (numeric_limits<float>::min)();
        t_far = (numeric_limits<float>::max)();
        vec3 direction = r.direction();
        vec3 origin = r.origin();
        vec3 bl = vertex_l;
        vec3 bh = vertex_h;
        float array1[6];

        if(direction.x() == 0) {
            if((origin.x() < bl.x()) || (origin.x() > bh.x())) {
#if PARAMETRIC_SURFACE_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin X0 is not between the slabs. return false" <<endl;
#endif // PARAMETRIC_SURFACE_LOG
                return false;
            }
            array1[0] = (numeric_limits<float>::min)();
            array1[1] = (numeric_limits<float>::max)();
        }
        if(direction.y() == 0) {
            if((origin.y() < bl.y()) || (origin.y() > bh.y())) {
#if PARAMETRIC_SURFACE_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin Y0 is not between the slabs. return false" <<endl;
#endif // PARAMETRIC_SURFACE_LOG
                return false;
            }
            array1[2] = (numeric_limits<float>::min)();
            array1[3] = (numeric_limits<float>::max)();
        }
        if(direction.z() == 0) {
            if((origin.z() < bl.z()) || (origin.z() > bh.z())) {
#if PARAMETRIC_SURFACE_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin Z0 is not between the slabs. return false" <<endl;
#endif // PARAMETRIC_SURFACE_LOG
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
#if PARAMETRIC_SURFACE_LOG == 1
            std::cout << "array1[" << i << "]:" << array1[i] <<endl;
            std::cout << "array1[" << i+1 << "]:" << array1[i+1] <<endl;
#endif // PARAMETRIC_SURFACE_LOG
            if(array1[i] >= t_near) {t_near = array1[i];}
            if(array1[i+1] <= t_far) {t_far = array1[i+1];}
            if(t_near > t_far) {
#if PARAMETRIC_SURFACE_LOG == 1
            std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_near > t_far. return false" <<endl;
#endif // PARAMETRIC_SURFACE_LOG
                return false;
            }
            if(t_far < 0) {
#if PARAMETRIC_SURFACE_LOG == 1
            std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_far < 0. return false" <<endl;
#endif // VEC3_LOG
                return false;
            }
        }
        if (t_near != t_near) {
            t_near = t_near * 1;
        }
        return true;
}


#if TYPE == 3
bool bezier3_parametric(const float matrix_x[4][4], const float matrix_y[4][4], const float matrix_z[4][4], const float u, const float v, vec3& parametric) {
        float temp_x = 0.0;
        float temp_y = 0.0;
        float temp_z = 0.0;
        float pow_u_i, pow_v_j;
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                pow_u_i = pow(u, i);
                pow_v_j = pow(v, j);
                temp_x = temp_x + matrix_x[i][j]*pow_u_i*pow_v_j;
                temp_y = temp_y + matrix_y[i][j]*pow_u_i*pow_v_j;
                temp_z = temp_z + matrix_z[i][j]*pow_u_i*pow_v_j;
            }
        }
        parametric = vec3(temp_x, temp_y, temp_z);
        return true;
}

bool bezier3_parametric_f1(const ray& r, const float matrix_x[4][4], const float matrix_y[4][4], const float matrix_z[4][4], const float u, const float v, vec3& fu, vec3& fv, vec3& am) {
        float fu_temp_x = 0.0;
        float fu_temp_y = 0.0;
        float fu_temp_z = 0.0;
        float fv_temp_x = 0.0;
        float fv_temp_y = 0.0;
        float fv_temp_z = 0.0;
        float pow_u_i, pow_v_j;
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                if (i == 0) {
                    fu_temp_x = fu_temp_x + 0.0;
                    fu_temp_y = fu_temp_y + 0.0;
                    fu_temp_z = fu_temp_y + 0.0;
                }
                else {
                    pow_u_i = pow(u, i-1);
                    pow_v_j = pow(v, j);
                    fu_temp_x = fu_temp_x + i*matrix_x[i][j]*pow_u_i*pow_v_j;
                    fu_temp_y = fu_temp_y + i*matrix_y[i][j]*pow_u_i*pow_v_j;
                    fu_temp_z = fu_temp_z + i*matrix_z[i][j]*pow_u_i*pow_v_j;
                }
                if (j == 0) {
                    fv_temp_x = fv_temp_x + 0.0;
                    fv_temp_y = fv_temp_y + 0.0;
                    fv_temp_z = fv_temp_z + 0.0;
                }
                else {
                    pow_u_i = pow(u, i);
                    pow_v_j = pow(v, j-1);
                    fv_temp_x = fv_temp_x + j*matrix_x[i][j]*pow_u_i*pow_v_j;
                    fv_temp_y = fv_temp_y + j*matrix_y[i][j]*pow_u_i*pow_v_j;
                    fv_temp_z = fv_temp_z + j*matrix_z[i][j]*pow_u_i*pow_v_j;
                }
            }
        }
        fu = vec3(fu_temp_x, fu_temp_y, fu_temp_z);
        fv = vec3(fv_temp_x, fv_temp_y, fv_temp_z);
        am = vec3(-r.direction().x(), -r.direction().y(), -r.direction().z());
        return true;
}

bool bezier3_parametric_f1_f2(const float matrix_x[4][4], const float matrix_y[4][4], const float matrix_z[4][4], const float u, const float v, vec3& fu, vec3& fv, vec3& fuu, vec3& fuv, vec3& fvv) {
        float fu_temp_x = 0.0;
        float fu_temp_y = 0.0;
        float fu_temp_z = 0.0;
        float fv_temp_x = 0.0;
        float fv_temp_y = 0.0;
        float fv_temp_z = 0.0;
        float fuu_temp_x = 0.0;
        float fuu_temp_y = 0.0;
        float fuu_temp_z = 0.0;
        float fvv_temp_x = 0.0;
        float fvv_temp_y = 0.0;
        float fvv_temp_z = 0.0;
        float fuv_temp_x = 0.0;
        float fuv_temp_y = 0.0;
        float fuv_temp_z = 0.0;
        float pow_u_i, pow_v_j, pow_u_i2, pow_v_j2;
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                if (i == 0) {
                    fu_temp_x = fu_temp_x + 0.0;
                    fu_temp_y = fu_temp_y + 0.0;
                    fu_temp_z = fu_temp_z + 0.0;
                }
                else {
                    pow_u_i = pow(u, i-1);
                    pow_v_j = pow(v, j);
                    fu_temp_x = fu_temp_x + i*matrix_x[i][j]*pow_u_i*pow_v_j;
                    fu_temp_y = fu_temp_y + i*matrix_y[i][j]*pow_u_i*pow_v_j;
                    fu_temp_z = fu_temp_z + i*matrix_z[i][j]*pow_u_i*pow_v_j;
                    if (i == 1) {
                        fuu_temp_x = fuu_temp_x + 0.0;
                        fuu_temp_y = fuu_temp_y + 0.0;
                        fuu_temp_z = fuu_temp_z + 0.0;
                    }
                    else {
                        pow_u_i2 = pow(u, i-2);
                        fuu_temp_x = fuu_temp_x + i*(i-1)*matrix_x[i][j]*pow_u_i2*pow_v_j;
                        fuu_temp_y = fuu_temp_y + i*(i-1)*matrix_x[i][j]*pow_u_i2*pow_v_j;
                        fuu_temp_z = fuu_temp_z + i*(i-1)*matrix_x[i][j]*pow_u_i2*pow_v_j;
                    }
                }

                if (j == 0) {
                    fv_temp_x = fv_temp_x + 0.0;
                    fv_temp_y = fv_temp_y + 0.0;
                    fv_temp_z = fv_temp_z + 0.0;
                }
                else {
                    pow_u_i = pow(u, i);
                    pow_v_j = pow(v, j-1);
                    fv_temp_x = fv_temp_x + j*matrix_x[i][j]*pow_u_i*pow_v_j;
                    fv_temp_y = fv_temp_y + j*matrix_y[i][j]*pow_u_i*pow_v_j;
                    fv_temp_z = fv_temp_z + j*matrix_z[i][j]*pow_u_i*pow_v_j;
                    if (j == 1) {
                        fvv_temp_x = fvv_temp_x + 0.0;
                        fvv_temp_y = fvv_temp_y + 0.0;
                        fvv_temp_z = fvv_temp_z + 0.0;
                    }
                    else {
                        pow_v_j2 = pow(v, j-2);
                        fvv_temp_x = fvv_temp_x + j*(j-1)*matrix_x[i][j]*pow_u_i*pow_v_j2;
                        fvv_temp_y = fvv_temp_y + j*(j-1)*matrix_x[i][j]*pow_u_i*pow_v_j2;
                        fvv_temp_z = fvv_temp_z + j*(j-1)*matrix_x[i][j]*pow_u_i*pow_v_j2;
                    }
                }

                if ((i == 0) || (j == 0)) {
                    fuv_temp_x = fuv_temp_x + 0.0;
                    fuv_temp_y = fuv_temp_y + 0.0;
                    fuv_temp_z = fuv_temp_z + 0.0;
                }
                else {
                    pow_u_i = pow(u, i-1);
                    pow_v_j = pow(v, j-1);
                    fuv_temp_x = fuv_temp_x + i*j*matrix_x[i][j]*pow_u_i*pow_v_j;
                    fuv_temp_y = fuv_temp_y + i*j*matrix_y[i][j]*pow_u_i*pow_v_j;
                    fuv_temp_z = fuv_temp_z + i*j*matrix_z[i][j]*pow_u_i*pow_v_j;
                }
            }
        }
        fu = vec3(fu_temp_x, fu_temp_y, fu_temp_z);
        fv = vec3(fv_temp_x, fv_temp_y, fv_temp_z);
        fuu = vec3(fuu_temp_x, fuu_temp_y, fuu_temp_z);
        fvv = vec3(fvv_temp_x, fvv_temp_y, fvv_temp_z);
        fuv = vec3(fuv_temp_x, fuv_temp_y, fuv_temp_z);
        return true;
}

bool bezier3_parametric_bl_bh(const float matrix_x[4][4], const float matrix_y[4][4], const float matrix_z[4][4], const float u1, const float u2, const float v1, const float v2, vec3& bl, vec3& bh) {
        vec3 p[4];
        float x[4], y[4], z[4];
        float temp;
        bezier3_parametric(matrix_x, matrix_y, matrix_z, u1, v1, p[0]);
        bezier3_parametric(matrix_x, matrix_y, matrix_z, u1, v2, p[1]);
        bezier3_parametric(matrix_x, matrix_y, matrix_z, u2, v1, p[2]);
        bezier3_parametric(matrix_x, matrix_y, matrix_z, u2, v2, p[3]);
        for (int i=0; i<4; i++) {
            x[i] = p[i].x();
            y[i] = p[i].y();
            z[i] = p[i].z();
        }
        for (int i=0; i<3; i++) {
            for (int j=i+1; j<4; j++) {
                if (x[i] > x[j]) {
                    temp = x[i];
                    x[i] = x[j];
                    x[j] = temp;
                }
                if (y[i] > y[j]) {
                    temp = y[i];
                    y[i] = y[j];
                    y[j] = temp;
                }
                if (z[i] > z[j]) {
                    temp = z[i];
                    z[i] = z[j];
                    z[j] = temp;
                }
            }
        }

        bl = vec3(x[0], y[0], z[0]);
        bh = vec3(x[3], y[3], z[3]);

        return true;
}
#endif // TYPE


#if TYPE == 2
bool horn_parametric(const vec3 cen, const float a, const float b, const float c, const float u, const float v, vec3& parametric) {
        float sin_u = sin(2*M_PI*u);
        float sin_v = sin(v);
        float cos_u = cos(2*M_PI*u);
        float cos_v = cos(v);
        sin_u = (fabs(sin_u)<1e-6)? 0.0:sin_u;
        sin_v = (fabs(sin_v)<1e-6)? 0.0:sin_v;
        cos_u = (fabs(cos_u)<1e-6)? 0.0:cos_u;
        cos_v = (fabs(cos_v)<1e-6)? 0.0:cos_v;
        parametric = vec3(a*(2.0+u*cos_v)*sin_u+cen.x(),
                          c*((2+u*cos_v)*cos_u+2*u)+cen.y(),
                          b*u*sin_v+cen.z());
        return true;
}

bool horn_parametric_f1(const ray& r, const float a, const float b, const float c, const float u, const float v, vec3& fu, vec3& fv, vec3& am) {
        float sin_u = sin(2*M_PI*u);
        float sin_v = sin(v);
        float cos_u = cos(2*M_PI*u);
        float cos_v = cos(v);
        sin_u = (fabs(sin_u)<1e-6)? 0.0:sin_u;
        sin_v = (fabs(sin_v)<1e-6)? 0.0:sin_v;
        cos_u = (fabs(cos_u)<1e-6)? 0.0:cos_u;
        cos_v = (fabs(cos_v)<1e-6)? 0.0:cos_v;
        fu = vec3(a*(2.0*M_PI*(2.0+u*cos_v)*cos_u+cos_v*sin_u),
                  c*(-2.0*M_PI*(2.0+u*cos_v)*sin_u+cos_v*cos_u+2.0),
                  b*(sin_v));
        fv = vec3(a*(-u*sin_u*sin_v),
                  c*(-u*cos_u*sin_v),
                  b*(u*cos_v));
        am = vec3(-r.direction().x(), -r.direction().y(), -r.direction().z());
        return true;
}

bool horn_parametric_f1_f2(const float a, const float b, const float c, const float u, const float v, vec3& fu, vec3& fv, vec3& fuu, vec3& fuv, vec3& fvv) {
        float sin_u = sin(2*M_PI*u);
        float sin_v = sin(v);
        float cos_u = cos(2*M_PI*u);
        float cos_v = cos(v);
        sin_u = (fabs(sin_u)<1e-6)? 0.0:sin_u;
        sin_v = (fabs(sin_v)<1e-6)? 0.0:sin_v;
        cos_u = (fabs(cos_u)<1e-6)? 0.0:cos_u;
        cos_v = (fabs(cos_v)<1e-6)? 0.0:cos_v;
        fu = vec3(a*(2.0*M_PI*(2.0+u*cos_v)*cos_u+cos_v*sin_u),
                  c*(-2.0*M_PI*(2.0+u*cos_v)*sin_u+cos_v*cos_u+2.0),
                  b*(sin_v));
        fv = vec3(a*(-u*sin_u*sin_v),
                  c*(-u*cos_u*sin_v),
                  b*(u*cos_v));
        fuu = vec3(a*(-4.0*M_PI*M_PI*(2.0+u*cos_v)*sin_u+4.0*M_PI*cos_v*cos_u),
                   c*(-4.0*M_PI*(cos_v*sin_u+M_PI*(2.0+u*cos_v)*cos_u)),
                   0.0);
        fuv = vec3(-a*sin_v*(sin_u+2*M_PI*u*cos_u),
                   -c*sin_v*(cos_u-2*M_PI*u*sin_u),
                   b*cos_v);
        fvv = vec3(-a*u*sin_u*cos_v,
                   -c*u*cos_u*cos_v,
                   -b*u*sin_v);
        return true;
}

bool horn_parametric_bl_bh(const vec3 cen, const float a, const float b, const float c, const float u1, const float u2, const float v1, const float v2, vec3& bl, vec3& bh) {
        float x[4], z[4], y[2];
        float sin_u1 = sin(2*M_PI*u1);
        float sin_u2 = sin(2*M_PI*u2);
        float sin_v1 = sin(v1);
        float sin_v2 = sin(v2);
        float cos_u1 = cos(2*M_PI*u1);
        float cos_u2 = cos(2*M_PI*u2);
        float cos_v1 = cos(v1);
        float cos_v2 = cos(v2);
        sin_u1 = (fabs(sin_u1)<1e-6)? 0.0:sin_u1;
        sin_u2 = (fabs(sin_u2)<1e-6)? 0.0:sin_u2;
        sin_v1 = (fabs(sin_v1)<1e-6)? 0.0:sin_v1;
        sin_v2 = (fabs(sin_v2)<1e-6)? 0.0:sin_v2;
        cos_u1 = (fabs(cos_u1)<1e-6)? 0.0:cos_u1;
        cos_u2 = (fabs(cos_u2)<1e-6)? 0.0:cos_u2;
        cos_v1 = (fabs(cos_v1)<1e-6)? 0.0:cos_v1;
        cos_v2 = (fabs(cos_v2)<1e-6)? 0.0:cos_v2;

        float temp;
        x[0] = a*((2+u1*cos_v1)*sin_u1)+cen.x();
        x[1] = a*((2+u1*cos_v2)*sin_u1)+cen.x();
        x[2] = a*((2+u2*cos_v1)*sin_u2)+cen.x();
        x[3] = a*((2+u2*cos_v2)*sin_u2)+cen.x();
        for (int i=0; i<3; i++) {
            for (int j=i+1; j<4; j++) {
                if (x[i] > x[j]) {
                    temp = x[i];
                    x[i] = x[j];
                    x[j] = temp;
                }
            }
        }
        y[0] = c*((2+u1*cos_v1)*cos_u1+2*u1)+cen.y();
        y[1] = c*((2+u1*cos_v2)*cos_u1+2*u1)+cen.y();
        y[2] = c*((2+u2*cos_v1)*cos_u2+2*u2)+cen.y();
        y[3] = c*((2+u2*cos_v2)*cos_u2+2*u2)+cen.y();
        for (int i=0; i<3; i++) {
            for (int j=i+1; j<4; j++) {
                if (y[i] > y[j]) {
                    temp = y[i];
                    y[i] = y[j];
                    y[j] = temp;
                }
            }
        }
        z[0] = b*u1*sin_v1+cen.z();
        z[1] = b*u1*sin_v2+cen.z();
        z[2] = b*u2*sin_v1+cen.z();
        z[3] = b*u2*sin_v2+cen.z();
        for (int i=0; i<3; i++) {
            for (int j=i+1; j<4; j++) {
                if (z[i] > z[j]) {
                    temp = z[i];
                    z[i] = z[j];
                    z[j] = temp;
                }
            }
        }

        bl = vec3(x[0], y[0], z[0]);
        bh = vec3(x[3], y[3], z[3]);

        return true;
}
#endif // TYPE

#if TYPE == 1
bool sphere_parametric(const vec3 cen, const float a, const float b, const float c, const float u, const float v, vec3& parametric) {
        float sin_u = sin(u);
        float sin_v = sin(v);
        float cos_u = cos(u);
        float cos_v = cos(v);
        sin_u = (fabs(sin_u)<1e-6)? 0.0:sin_u;
        sin_v = (fabs(sin_v)<1e-6)? 0.0:sin_v;
        cos_u = (fabs(cos_u)<1e-6)? 0.0:cos_u;
        cos_v = (fabs(cos_v)<1e-6)? 0.0:cos_v;
        parametric = vec3(a*sin_u*sin_v+cen.x(), c*cos_u+cen.y(), b*sin_u*cos_v+cen.z());
        return true;
}

bool sphere_parametric_f1(const ray& r, const float a, const float b, const float c, const float u, const float v, vec3& fu, vec3& fv, vec3& am) {
        float sin_u = sin(u);
        float sin_v = sin(v);
        float cos_u = cos(u);
        float cos_v = cos(v);
        sin_u = (fabs(sin_u)<1e-6)? 0.0:sin_u;
        sin_v = (fabs(sin_v)<1e-6)? 0.0:sin_v;
        cos_u = (fabs(cos_u)<1e-6)? 0.0:cos_u;
        cos_v = (fabs(cos_v)<1e-6)? 0.0:cos_v;
        fu = vec3(a*cos_u*sin_v, -c*sin_u, b*cos_u*cos_v);
        fv = vec3(a*sin_u*cos_v, 0, -b*sin_u*sin_v);
        am = vec3(-r.direction().x(), -r.direction().y(), -r.direction().z());
        return true;
}

bool sphere_parametric_f1_f2(const float a, const float b, const float c, const float u, const float v, vec3& fu, vec3& fv, vec3& fuu, vec3& fuv, vec3& fvv) {
        float sin_u = sin(u);
        float sin_v = sin(v);
        float cos_u = cos(u);
        float cos_v = cos(v);
        sin_u = (fabs(sin_u)<1e-6)? 0.0:sin_u;
        sin_v = (fabs(sin_v)<1e-6)? 0.0:sin_v;
        cos_u = (fabs(cos_u)<1e-6)? 0.0:cos_u;
        cos_v = (fabs(cos_v)<1e-6)? 0.0:cos_v;
        fu = vec3(a*cos_u*sin_v, -c*sin_u, b*cos_u*cos_v);
        fv = vec3(a*sin_u*cos_v, 0, -b*sin_u*sin_v);
        fuu = vec3(-a*sin_u*sin_v, -c*cos_u, -b*sin_u*cos_v);
        fuv = vec3(a*cos_u*cos_v, 0, -b*cos_u*sin_v);
        fvv = vec3(-a*sin_u*sin_v, 0, -b*sin_u*cos_v);
        return true;
}

bool sphere_parametric_bl_bh(const vec3 cen, const float a, const float b, const float c, const float u1, const float u2, const float v1, const float v2, vec3& bl, vec3& bh) {
        float x[4], z[4], y[2];
        float sin_u1 = sin(u1);
        float sin_u2 = sin(u2);
        float sin_v1 = sin(v1);
        float sin_v2 = sin(v2);
        float cos_u1 = cos(u1);
        float cos_u2 = cos(u2);
        float cos_v1 = cos(v1);
        float cos_v2 = cos(v2);
        sin_u1 = (fabs(sin_u1)<1e-6)? 0.0:sin_u1;
        sin_u2 = (fabs(sin_u2)<1e-6)? 0.0:sin_u2;
        sin_v1 = (fabs(sin_v1)<1e-6)? 0.0:sin_v1;
        sin_v2 = (fabs(sin_v2)<1e-6)? 0.0:sin_v2;
        cos_u1 = (fabs(cos_u1)<1e-6)? 0.0:cos_u1;
        cos_u2 = (fabs(cos_u2)<1e-6)? 0.0:cos_u2;
        cos_v1 = (fabs(cos_v1)<1e-6)? 0.0:cos_v1;
        cos_v2 = (fabs(cos_v2)<1e-6)? 0.0:cos_v2;

        float temp;
        x[0] = a*sin_u1*sin_v1+cen.x();
        x[1] = a*sin_u1*sin_v2+cen.x();
        x[2] = a*sin_u2*sin_v1+cen.x();
        x[3] = a*sin_u2*sin_v2+cen.x();
        for (int i=0; i<3; i++) {
            for (int j=i+1; j<4; j++) {
                if (x[i] > x[j]) {
                    temp = x[i];
                    x[i] = x[j];
                    x[j] = temp;
                }
            }
        }
        z[0] = b*sin_u1*cos_v1+cen.z();
        z[1] = b*sin_u1*cos_v2+cen.z();
        z[2] = b*sin_u2*cos_v1+cen.z();
        z[3] = b*sin_u2*cos_v2+cen.z();
        for (int i=0; i<3; i++) {
            for (int j=i+1; j<4; j++) {
                if (z[i] > z[j]) {
                    temp = z[i];
                    z[i] = z[j];
                    z[j] = temp;
                }
            }
        }
        y[0] = ((c*cos_u1+cen.y())<(c*cos_u2+cen.y()))? (c*cos_u1+cen.y()):(c*cos_u2+cen.y());
        y[1] = ((c*cos_u1+cen.y())<(c*cos_u2+cen.y()))? (c*cos_u2+cen.y()):(c*cos_u1+cen.y());

        bl = vec3(x[0], y[0], z[0]);
        bh = vec3(x[3], y[1], z[3]);

        return true;
}
#endif // TYPE


#if ((TYPE == 1) || (TYPE ==2))
int get_root_by_newton_iteration(const ray& r, vec3 cen, float a, float b, float c, vec3 x0, float tao, float uid, float vid, vec3& root) {
#endif // TYPE
#if TYPE == 3
int get_root_by_newton_iteration(const ray& r, const float matrix_x[4][4], const float matrix_y[4][4], const float matrix_z[4][4], vec3 x0, float tao, float uid, float vid, vec3& root) {
#endif // TYPE
        vec3 fun0, dis0, fu, fv, am, fun1, dis1, fuu, fuv, fvv, normal;
        float m_d0[3][1], m_j[3][3], m_j_i[3][3], m_j_i_d0[3][1], m_x_k0[3][1], m_x_k1[3][1];
        float d0, d1, d2, h;
        get_matrix_3_1(x0, m_x_k0);
#if TYPE == 1
        sphere_parametric(cen, a, b, c, m_x_k0[0][0], m_x_k0[1][0], fun0);
#endif // TYPE
#if TYPE == 2
        horn_parametric(cen, a, b, c, m_x_k0[0][0], m_x_k0[1][0], fun0);
#endif // TYPE
#if TYPE == 3
        bezier3_parametric(matrix_x, matrix_y, matrix_z, m_x_k0[0][0], m_x_k0[1][0], fun0);
#endif // TYPE
        dis0 = fun0 - r.point_at_parameter(m_x_k0[2][0]);

        for (int i=0; i<20; i++) {
#if TYPE == 1
            sphere_parametric_f1(r, a, b, c, m_x_k0[0][0], m_x_k0[1][0], fu, fv, am);
#endif // TYPE
#if TYPE == 2
            horn_parametric_f1(r, a, b, c, m_x_k0[0][0], m_x_k0[1][0], fu, fv, am);
#endif // TYPE
#if TYPE == 3
            bezier3_parametric_f1(r, matrix_x, matrix_y, matrix_z, m_x_k0[0][0], m_x_k0[1][0], fu, fv, am);
#endif // TYPE
            get_matrix_3_1(dis0, m_d0);
            get_matrix_3_3(fu, fv, am, m_j);
            if (get_matrix_inverse_3_3(m_j, m_j_i)) {
                matrix_3_3_multiply_3_1(m_j_i, m_d0, m_j_i_d0);
                matrix_3_1_minus_3_1(m_x_k0, m_j_i_d0, m_x_k1);
//                if (m_x_k1[2][0] > 1e-16)
                {
#if TYPE == 1
                    sphere_parametric(cen, a, b, c, m_x_k1[0][0], m_x_k1[1][0], fun1);
#endif // TYPE
#if TYPE == 2
                    horn_parametric(cen, a, b, c, m_x_k1[0][0], m_x_k1[1][0], fun1);
#endif // TYPE
#if TYPE == 3
                    bezier3_parametric(matrix_x, matrix_y, matrix_z, m_x_k1[0][0], m_x_k1[1][0], fun1);
#endif // TYPE
                    dis1 = fun1 - r.point_at_parameter(m_x_k1[2][0]);
                    if (dis1.length() < tao) {
#if TYPE == 1
                        sphere_parametric_f1_f2(a, b, c, m_x_k1[0][0], m_x_k1[1][0], fu, fv, fuu, fuv, fvv);
#endif // TYPE
#if TYPE == 2
                        horn_parametric_f1_f2(a, b, c, m_x_k1[0][0], m_x_k1[1][0], fu, fv, fuu, fuv, fvv);
#endif // TYPE
#if TYPE == 3
                        bezier3_parametric_f1_f2(matrix_x, matrix_y, matrix_z, m_x_k1[0][0], m_x_k1[1][0], fu, fv, fuu, fuv, fvv);
#endif // TYPE
                        normal = unit_vector(cross(fu, fv));
                        d0 = dot(normal, fuu);
                        d1 = dot(normal, fuv);
                        d2 = dot(normal, fvv);
                        h = 0.5*(d0*uid*uid+2*d1*uid*vid+d2*vid*vid);
                        if (fabs(h) < tao) {
                            if (root.z() == 0.0) {
                                root = vec3(m_x_k1[0][0], m_x_k1[1][0], m_x_k1[2][0]);
                            }
                            else {
                                if (m_x_k1[2][0] < root.z()) {
                                    root = vec3(m_x_k1[0][0], m_x_k1[1][0], m_x_k1[2][0]);
                                }
                            }
                            return 1;
                        }
                        else {
                            root = vec3(m_x_k1[0][0], m_x_k1[1][0], m_x_k1[2][0]);
                            return 2;//0;
                        }
                    }
                    else {
                        m_x_k0[0][0] = m_x_k1[0][0];
                        m_x_k0[1][0] = m_x_k1[1][0];
                        m_x_k0[2][0] = m_x_k1[2][0];
                        dis0 = dis1;
                    }
                }
            }
            else {
                return 0;
            }
        }
        return 0;
}


bool parametric_surface::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if PARAMETRIC_SURFACE_LOG == 1
        std::cout << "-------------parametric_surface::hit----------------" << endl;
#endif // PARAMETRIC_SURFACE_LOG
#if TYPE == 1
        #define UN 4
        #define VN 8
#endif // TYPE
#if TYPE == 2
        #define UN 200
        #define VN 200
#endif // TYPE
#if TYPE == 3
        #define UN 16
        #define VN 16
#endif // TYPE
        float ui1, ui2, vi1, vi2, t_near, t_far, tao, u1, v1, uid, vid;
        vec3 bl, bh, fu, fv, am, x0;
        vec3 root = vec3(0.0, 0.0, 0.0);
        int rn = 2;
        int root_num = 0;
        bool hitted1 = 0;
        float t_near_a, t_far_a;
        t_near_a = -1.0;
        t_far_a = t_near_a;
        bool hitted2 = 0;
        float ui1_a[4], ui2_a[4], vi1_a[4], vi2_a[4];
#if TYPE == 1
        bl = vec3(center.x()-intercept_x, center.y()-intercept_y, center.z()-intercept_z);
        bh = vec3(center.x()+intercept_x, center.y()+intercept_y, center.z()+intercept_z);
#endif // TYPE
#if TYPE == 2
        bl = vec3(center.x()-intercept_x*2, center.y()-intercept_y*4, center.z()-intercept_z*2);
        bh = vec3(center.x()+intercept_x*2, center.y()+intercept_y*4, center.z()+intercept_z*2);
#endif // TYPE
#if TYPE == 3
        bl = vec3(bezier3_bl.x(), bezier3_bl.y(), bezier3_bl.z());
        bh = vec3(bezier3_bh.x(), bezier3_bh.y(), bezier3_bh.z());
#endif // TYPE
        if(ray_hit_box_para(r, bl, bh, t_near, t_far)) {
            for (int i=0; i<UN; i++) {
#if TYPE == 1
                ui1 = float(i)*M_PI/float(UN);
                ui2 = float(i+1)*M_PI/float(UN);
#endif // TYPE
#if TYPE == 2
                ui1 = float(i)/float(UN);
                ui2 = float(i+1)/float(UN);
#endif // TYPE
#if TYPE == 3
                ui1 = float(i)/float(UN);
                ui2 = float(i+1)/float(UN);
#endif // TYPE
                for (int j=0; j<VN; j++) {
#if (TYPE == 1) || (TYPE == 2)
                    vi1 = float(j)*2*M_PI/float(VN);
                    vi2 = float(j+1)*2*M_PI/float(VN);
#endif // TYPE
#if TYPE == 3
                    vi1 = float(j)/float(VN);
                    vi2 = float(j+1)/float(VN);
#endif // TYPE
                    u1 = (ui1+ui2)/2;
                    v1 = (vi1+vi2)/2;
#if TYPE == 1
                    sphere_parametric_bl_bh(center, intercept_x, intercept_y, intercept_z, ui1, ui2, vi1, vi2, bl, bh);
#endif // TYPE
#if TYPE == 2
                    horn_parametric_bl_bh(center, intercept_x, intercept_y, intercept_z, ui1, ui2, vi1, vi2, bl, bh);
#endif // TYPE
#if TYPE == 3
                    bezier3_parametric_bl_bh(matrix_c_x, matrix_c_y, matrix_c_z, ui1, ui2, vi1, vi2, bl, bh);
#endif // TYPE
                    hitted1 = ray_hit_box_para(r, bl, bh, t_near, t_far);
                    if(hitted1 && (t_near > 1e-7)) {
                        while (rn == 2) {
                            tao = fabs(2*rho*t_near*sin(theta*(M_PI/180)/2)/num);
                            x0 = vec3(u1, v1, t_near);
                            uid = ui2-ui1;
                            vid = vi2-vi1;
#if (TYPE == 1) || (TYPE == 2)
                            rn = get_root_by_newton_iteration(r, center, intercept_x, intercept_y, intercept_z, x0, tao, uid, vid, root);
#endif // TYPE
#if TYPE == 3
                            rn = get_root_by_newton_iteration(r, matrix_c_x, matrix_c_y, matrix_c_z, x0, tao, uid, vid, root);
#endif // TYPE
                            if (rn == 1) {
                                root_num ++ ;
                            }
                            else if (rn == 2) {
                                ui1_a[0] = ui1;
                                ui2_a[0] = ui1+uid/2;
                                vi1_a[0] = vi1;
                                vi2_a[0] = vi1+vid/2;

                                ui1_a[1] = ui1;
                                ui2_a[1] = ui1+uid/2;
                                vi1_a[1] = vi1+vid/2;
                                vi2_a[1] = vi2;

                                ui1_a[2] = ui1+uid/2;
                                ui2_a[2] = ui2;
                                vi1_a[2] = vi1;
                                vi2_a[2] = vi1+vid/2;

                                ui1_a[3] = ui1+uid/2;
                                ui2_a[3] = ui2;
                                vi1_a[3] = vi1+vid/2;
                                vi2_a[3] = vi2;

                                for (int j=0; j<4; j++) {
                                    ui1 = ui1_a[j];
                                    ui2 = ui2_a[j];
                                    vi1 = vi1_a[j];
                                    vi2 = vi2_a[j];
#if TYPE == 1
                                    sphere_parametric_bl_bh(center, intercept_x, intercept_y, intercept_z, ui1, ui2, vi1, vi2, bl, bh);
#endif // TYPE
#if TYPE == 2
                                    horn_parametric_bl_bh(center, intercept_x, intercept_y, intercept_z, ui1, ui2, vi1, vi2, bl, bh);
#endif // TYPE
#if TYPE == 3
                                    bezier3_parametric_bl_bh(matrix_c_x, matrix_c_y, matrix_c_z, ui1, ui2, vi1, vi2, bl, bh);
#endif // TYPE
                                    hitted2 = ray_hit_box_para(r, bl, bh, t_near_a, t_far_a);
                                    if(hitted2 && (t_near_a > 1e-7)) {
    //                                    if (t_near_a == t_near) {
                                            break;
    //                                    }
                                    }
                                }
                                if (hitted2) {
                                    t_near = root.z(); //t_near_a;
                                    u1 = root.x();
                                    v1 = root.y();
                                }
                                else {
                                    break;
                                }
                            }
                            else if (rn == 0) {
                                return false;
                            }
                        }
                        if (root_num == 2) {
                            break;
                        }
                   }
                }
                if (root_num == 2) {
                    break;
                }
            }
        }
        else {
            return false;
        }


        if (root.z() < t_max && root.z() > t_min) {
            rec.t = root.z();
            rec.p = r.point_at_parameter(rec.t);
//            vec3 pc = rec.p - center;
#if TYPE == 1
            sphere_parametric_f1(r, intercept_x, intercept_y, intercept_z, root.x(), root.y(), fu, fv, am);
#endif // TYPE
#if TYPE == 2
            horn_parametric_f1(r, intercept_x, intercept_y, intercept_z, root.x(), root.y(), fu, fv, am);
#endif // TYPE
#if TYPE == 3
            bezier3_parametric_f1(r, matrix_c_x, matrix_c_y, matrix_c_z, root.x(), root.y(), fu, fv, am);
#endif // TYPE
            rec.normal = unit_vector(cross(fu, fv));
            if(dot(r.direction(), rec.normal) > 0) {
                rec.normal = - rec.normal;
            }
            rec.mat_ptr = ma;
            rec.u = -1.0;
            rec.v = -1.0;
            return true;
        }
        return false;
}
