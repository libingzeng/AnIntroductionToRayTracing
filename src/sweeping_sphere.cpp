#include "sweeping_sphere.h"


bool sweeping_sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if SWEEPING_SPHERE_LOG == 1
        std::cout << "-------------sweeping_sphere::hit---1-------------" << endl;
#endif // SWEEPING_SPHERE_LOG
        float t_near, t_far;
        if (ray_hit_box_general(r, sweeping_bl, sweeping_bh, t_near, t_far)) {
            float xx0 = r.origin().x();
            float yy0 = r.origin().y();
            float zz0 = r.origin().z();
            float xxd = r.direction().x();
            float yyd = r.direction().y();
            float zzd = r.direction().z();
            if (fabs(xxd) < 1e-6) {
                xxd = 1e-6;
            }
            if (fabs(yyd) < 1e-6) {
                yyd = 1e-6;
            }
            if (fabs(zzd) < 1e-6) {
                zzd = 1e-6;
            }
            float f20 = xxd*xxd + yyd*yyd + zzd*zzd;
            float f13, f12, f11, f10, aa2, aa1, aa0, cc5, cc4, cc3, cc2, cc1, cc0, dd4, dd3, dd2, dd1, dd0;
            float f06, f05, f04, f03, f02, f01, f00, bb5, bb4, bb3, bb2, bb1, bb0;
            float cx3, cx2, cx1, cx0;
            float cy3, cy2, cy1, cy0;
            float cz3, cz2, cz1, cz0;
            float cr3, cr2, cr1, cr0;
            float gg10[11], hh10[11], ff10[11], ee10[11], roots[11], roots_t[31][3], uuu, f0d, f1d, ttt, temp0, temp1, temp2;
            float tol = 1e-6;
            int num_roots_t = 0;
            float px, py, pz, cx, cy, cz;
            int num;

            for (int i=0; i<5; i++) {//5
                cx3 = matrix_c_cx[3][i];
                cx2 = matrix_c_cx[2][i];
                cx1 = matrix_c_cx[1][i];
                cx0 = matrix_c_cx[0][i];
                cy3 = matrix_c_cy[3][i];
                cy2 = matrix_c_cy[2][i];
                cy1 = matrix_c_cy[1][i];
                cy0 = matrix_c_cy[0][i];
                cz3 = matrix_c_cz[3][i];
                cz2 = matrix_c_cz[2][i];
                cz1 = matrix_c_cz[1][i];
                cz0 = matrix_c_cz[0][i];
                cr3 = matrix_c_r[3][i];
                cr2 = matrix_c_r[2][i];
                cr1 = matrix_c_r[1][i];
                cr0 = matrix_c_r[0][i];
                f13 = -2*(xxd*cx3 + yyd*cy3 + zzd*cz3);
                f12 = -2*(xxd*cx2 + yyd*cy2 + zzd*cz2);
                f11 = -2*(xxd*cx1 + yyd*cy1 + zzd*cz1);
                f10 = -2*(xxd*cx0 + yyd*cy0 + zzd*cz0) + 2*(xx0*xxd+yy0*yyd+zz0*zzd);

                aa2 = 3*f13;
                aa1 = 2*f12;
                aa0 = f11;
                f06 = cx3*cx3 + cy3*cy3 + cz3*cz3 - cr3*cr3;
                f05 = 2*(cx2*cx3 + cy2*cy3 + cz2*cz3 - cr2*cr3);
                f04 = 2*(cx1*cx3 + cy1*cy3 + cz1*cz3 - cr1*cr3) + (cx2*cx2 + cy2*cy2 + cz2*cz2 - cr2*cr2);
                f03 = 2*(cx1*cx2 + cy1*cy2 + cz1*cz2 - cr1*cr2) + 2*((cx0-xx0)*cx3 + (cy0-yy0)*cy3 + (cz0-zz0)*cz3 - cr0*cr3);
                f02 = 2*((cx0-xx0)*cx2 + (cy0-yy0)*cy2 + (cz0-zz0)*cz2 - cr0*cr2) + (cx1*cx1 + cy1*cy1 + cz1*cz1 - cr1*cr1);
                f01 = 2*((cx0-xx0)*cx1 + (cy0-yy0)*cy1 + (cz0-zz0)*cz1 - cr0*cr1);
                f00 = (cx0-xx0)*(cx0-xx0) + (cy0-yy0)*(cy0-yy0) + (cz0-zz0)*(cz0-zz0) - cr0*cr0;
                bb5 = 6*f06;
                bb4 = 5*f05;
                bb3 = 4*f04;
                bb2 = 3*f03;
                bb1 = 2*f02;
                bb0 = f01;
                hh10[0]  = f20*(bb5*bb5);
                hh10[1]  = f20*(2*bb4*bb5);
                hh10[2]  = f20*(2*bb3*bb5 + bb4*bb4);
                hh10[3]  = f20*(2*bb2*bb5 + 2*bb3*bb4);
                hh10[4]  = f20*(2*bb1*bb5 + 2*bb2*bb4 + bb3*bb3);
                hh10[5]  = f20*(2*bb0*bb5 + 2*bb1*bb4 + 2*bb2*bb3);
                hh10[6]  = f20*(2*bb0*bb4 + 2*bb1*bb3 + bb2*bb2);
                hh10[7]  = f20*(2*bb0*bb3 + 2*bb1*bb2);
                hh10[8]  = f20*(2*bb0*bb2 + bb1*bb1);
                hh10[9]  = f20*(2*bb0*bb1);
                hh10[10] = f20*(bb0*bb0);
                cc5 = f13*aa2;
                cc4 = f13*aa1 + f12*aa2;
                cc3 = f13*aa0 + f12*aa1 + f11*aa2;
                cc2 = f12*aa0 + f11*aa1 + f10*aa2;
                cc1 = f11*aa0 + f10*aa1;
                cc0 = f10*aa0;
                ff10[0]  = bb5*cc5;
                ff10[1]  = bb5*cc4 + bb4*cc5;
                ff10[2]  = bb5*cc3 + bb4*cc4 + bb3*cc5;
                ff10[3]  = bb5*cc2 + bb4*cc3 + bb3*cc4 + bb2*cc5;
                ff10[4]  = bb5*cc1 + bb4*cc2 + bb3*cc3 + bb2*cc4 + bb1*cc5;
                ff10[5]  = bb5*cc0 + bb4*cc1 + bb3*cc2 + bb2*cc3 + bb1*cc4+ bb0*cc5;
                ff10[6]  = bb4*cc0 + bb3*cc1 + bb2*cc2 + bb1*cc3 + bb0*cc4;
                ff10[7]  = bb3*cc0 + bb2*cc1 + bb1*cc2 + bb0*cc3;
                ff10[8]  = bb2*cc0 + bb1*cc1+ bb0*cc2;
                ff10[9]  = bb1*cc0 + bb0*cc1;
                ff10[10] = bb0*cc0;
                dd4 = aa2*aa2;
                dd3 = 2*aa1*aa2;
                dd2 = 2*aa0*aa2 + aa1*aa1;
                dd1 = 2*aa0*aa1;
                dd0 = aa0*aa0;
                ee10[0]  = dd4*f06;
                ee10[1]  = dd4*f05 + dd3*f06;
                ee10[2]  = dd4*f04 + dd3*f05 + dd2*f06;
                ee10[3]  = dd4*f03 + dd3*f04 + dd2*f05 + dd1*f06;
                ee10[4]  = dd4*f02 + dd3*f03 + dd2*f04 + dd1*f05 + dd0*f06;
                ee10[5]  = dd4*f01 + dd3*f02 + dd2*f03 + dd1*f04 + dd0*f05;
                ee10[6]  = dd4*f00 + dd3*f01 + dd2*f02 + dd1*f03 + dd0*f04;
                ee10[7]  = dd3*f00 + dd2*f01 + dd1*f02 + dd0*f03;
                ee10[8]  = dd2*f00 + dd1*f01 + dd0*f02;
                ee10[9]  = dd1*f00 + dd0*f01;
                ee10[10] = dd0*f00;
                for (int k=0; k<11; k++) {
                    gg10[k] = hh10[k] - ff10[k] + ee10[k];
                    if (fabs(gg10[k])<1e-6) {
                        gg10[k] = 1e-6;
                    }
                }
                roots_equation_10th(gg10, 0, 1, tol, roots);
                if (int(roots[0]) >= 1) {
                    for (int j=1; j<(int(roots[0])+1); j++) {
                        uuu = roots[j];//the real roots u which is the parameter of cubic b-spline.
                        if (fabs(uuu)>1e-6) {
                            f0d = bb5*pow(uuu,5) + bb4*pow(uuu,4) + bb3*pow(uuu,3) + bb2*pow(uuu,2) + bb1*uuu + bb0;
                            f1d = aa2*pow(uuu,2) + aa1*uuu + aa0;
                            if (fabs(f1d)>1e-6) {
                                ttt = -f0d/f1d;//note that b-spline is the trajectory of sphere center and the distance from hit-point to origin of ray is determined by this fomula
                                roots_t[num_roots_t+1][0] = ttt;//store the distance from hit-point to origin.
                                if (roots_t[num_roots_t+1][0] > 0) {
                                    roots_t[num_roots_t+1][1] = i;//store the b-spline curve number that the corresponding sphere center lies on
                                    roots_t[num_roots_t+1][2] = roots[j];//store the b-spline parameter of the center of the hitted sphere.
                                    num_roots_t ++;
                                }
                            }
                        }
                    }
                }
            }
            roots_t[0][0] = float(num_roots_t);
            if (roots_t[0][0] != 0.0) {
                for (int i=1; i<int(roots_t[0][0]); i++) {
                    for (int j=i+1; j<int(roots_t[0][0])+1; j++) {
                        if (roots_t[i][0] > roots_t[j][0]) {
                            temp0 = roots_t[i][0];
                            roots_t[i][0] = roots_t[j][0];
                            roots_t[j][0] = temp0;
                            temp1 = roots_t[i][1];
                            roots_t[i][1] = roots_t[j][1];
                            roots_t[j][1] = temp1;
                            temp2 = roots_t[i][2];
                            roots_t[i][2] = roots_t[j][2];
                            roots_t[j][2] = temp2;
                        }
                    }
                }
                if ((roots_t[1][0] < t_max) && (roots_t[1][0] > t_min)) {
                    rec.t = roots_t[1][0];
                    rec.p = r.point_at_parameter(rec.t);
                    px = rec.p.x();
                    py = rec.p.y();
                    pz = rec.p.z();
                    num = int(roots_t[1][1]);
                    uuu = roots_t[1][2];
                    cx = matrix_c_cx[0][num] + matrix_c_cx[1][num]*uuu + matrix_c_cx[2][num]*uuu*uuu + matrix_c_cx[3][num]*uuu*uuu*uuu;
                    cy = matrix_c_cy[0][num] + matrix_c_cy[1][num]*uuu + matrix_c_cy[2][num]*uuu*uuu + matrix_c_cy[3][num]*uuu*uuu*uuu;
                    cz = matrix_c_cz[0][num] + matrix_c_cz[1][num]*uuu + matrix_c_cz[2][num]*uuu*uuu + matrix_c_cz[3][num]*uuu*uuu*uuu;
                    rec.normal = unit_vector(vec3(px-cx, py-cy, pz-cz));
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
