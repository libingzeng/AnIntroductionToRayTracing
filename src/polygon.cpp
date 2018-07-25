#include "polygon.h"
#include "vec2.h"
#include "log.h"

#include <iostream>
#include <fstream>

using namespace std;

bool in_polygon_test(vec2 *vertexes_uv, int number) {
        int sh, nsh;
        int nc = 0;
        if(vertexes_uv[0].v() < 0) { sh = -1;}
        else { sh = 1;}
        for(int j=0; j<number; j++) {
            if(vertexes_uv[j+1].v() < 0) { nsh = -1;}
            else { nsh = 1;}
            if(sh != nsh) {
                if((vertexes_uv[j].u() > 0) && (vertexes_uv[j+1].u() >0)) { nc++;}
                else
                    if((vertexes_uv[j].u() > 0) || (vertexes_uv[j+1].u() >0)) {
                        if(vertexes_uv[j].u() - (vertexes_uv[j].v())*(vertexes_uv[j+1].u()-vertexes_uv[j].u())/(vertexes_uv[j+1].v()-vertexes_uv[j].v()) > 0) { nc++;}
                    }
            }
            sh = nsh;
        }
        if((nc)%(2)) {return true;}
        else {return false;}
}

bool in_polygon_test2(vec2 *vertexes_uv, int number) {
        int sh, nsh;
        int nc = 0;
        if(vertexes_uv[0].v() < 0) { sh = -1;}
        else { sh = 1;}
        for(int j=0; j<number; j++) {
            if(vertexes_uv[j+1].v() < 0) { nsh = -1;}
            else { nsh = 1;}
            if(sh != nsh) {
                if((vertexes_uv[j].u() > 0) && (vertexes_uv[j+1].u() >0)) {
                    if(vertexes_uv[j].v() > vertexes_uv[j+1].v()) { nc = nc + 1;}
                    else { nc = nc - 1;}
                }
                else
                    if((vertexes_uv[j].u() > 0) || (vertexes_uv[j+1].u() >0)) {
                        if(vertexes_uv[j].u() - (vertexes_uv[j].v())*(vertexes_uv[j+1].u()-vertexes_uv[j].u())/(vertexes_uv[j+1].v()-vertexes_uv[j].v()) > 0) {
                            if(vertexes_uv[j].v() > vertexes_uv[j+1].v()) { nc = nc + 1;}
                            else { nc = nc - 1;}
                        }
                    }
            }
            sh = nsh;
        }
        if(nc != 0) {return true;}
        else {return false;}
}

bool polygon::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
        vec3 poly_n;
        for(int i=0; i<number-2; i++) {
            poly_n = unit_vector(cross((vertexes[i]-vertexes[i+1]), (vertexes[i+1]-vertexes[i+2])));//determine the normal of the plane
            if (dot(poly_n, r.direction()) > 0) {
                poly_n = - poly_n;
            }
            if(!vector_equ(poly_n, vec3(0,0,0))) {
                break;
            }
        }
        float poly_d = -(dot(poly_n, vertexes[0]));//determine the distance from the origin to the plane
        float vd = dot(poly_n, r.direction());
        float v0 = -(dot(poly_n, r.origin()) + poly_d);

        if(vd == 0) {//the ray is parallel to the polygon plane
            return false;
        }
        else {
            float t = v0/vd;//determine t and intersection pi
            vec3 pi = r.point_at_parameter(t);

            /*find the dominant coordinate, X, Y, or Z?
            i=1: means that X is the dominant coordinate;
            i=2: means that Y is the dominant coordinate;
            i=3: means that Z is the dominant coordinate;
            */
            float temp = fabs(poly_n.x());
            int i = 1;
            if(temp <= fabs(poly_n.y())) {
                temp = fabs(poly_n.y());
                i++;
            }
            if(temp <= fabs(poly_n.z())) {
                i++;
            }

            /*throw the dorminant coordinate of 3-d vector, then we get 2-d vector in uv-plane*/
            vec2 vertexes_uv[number+1];
            switch (i) {
            case 1:
                for(int i=0; i<number; i++) {
                    vertexes_uv[i] = vec2(vertexes[i].y(),vertexes[i].z());
                }
                vertexes_uv[number] = vec2(pi.y(),pi.z());
                break;
            case 2:
                for(int i=0; i<number; i++) {
                    vertexes_uv[i] = vec2(vertexes[i].x(),vertexes[i].z());
                }
                vertexes_uv[number] = vec2(pi.x(),pi.z());
                break;
           case 3:
                for(int i=0; i<number; i++) {
                    vertexes_uv[i] = vec2(vertexes[i].x(),vertexes[i].y());
                }
                vertexes_uv[number] = vec2(pi.x(),pi.y());
                break;
            }

            /*move intersection uv-coordinate to origin.
            so all the vertexes substract intersection uv-coordinate.*/
            for(int i=0; i<number; i++) {
                vertexes_uv[i] = vertexes_uv[i] - vertexes_uv[number];
            }
            vertexes_uv[number] = vertexes_uv[0];
            //set the first vertex to the last position of the array, so that we get the whole vertexes loop
            if(in_polygon_test2(vertexes_uv,number)) {//check if the intersection locates inside the polygon or not
                if (t < t_max && t > t_min) {
                    rec.t = t;
                    rec.p = r.point_at_parameter(rec.t);
                    rec.normal = poly_n;
                    rec.mat_ptr = ma;

                    if (number <= 4) {
                        vec3 p00 = vertexes[0];
                        vec3 p10 = vertexes[1];
                        vec3 p11 = vertexes[2];
                        vec3 p01 = vertexes[0];
                        if (number == 4) {
                            p01 = vertexes[3];
                        }
                        vec3 pa = p00-p01+p11-p10;
                        vec3 pb = p10-p00;
                        vec3 pc = p01-p00;
                        vec3 pd = p00;
                        vec3 pn = rec.normal;
                        vec3 na = cross(pa, pn);
                        vec3 nb = cross(pb, pn);
                        vec3 nc = cross(pc, pn);
                        float du0 = dot(nc, pd);
                        float du1 = dot(na, pd) + dot(nc, pb);
                        float du2 = dot(na, pb);
                        float dv0 = dot(nb, pd);
                        float dv1 = dot(na, pd) + dot(nb, pc);
                        float dv2 = dot(na, pc);
                        vec3 pi = rec.p;
                        float Au = du2;
                        float Bu = du1 - dot(na, pi);
                        float Cu = du0 - dot(nc, pi);
                        float Av = dv2;
                        float Bv = dv1 - dot(na, pi);
                        float Cv = dv0 - dot(nb, pi);

                        if (Au == 0) {
                            rec.u = -Cu/Bu;
                        }
                        else {
                            float u_temp = (-Bu + sqrt(Bu*Bu-4*Au*Cu)) / (2*Au);
                            if ((u_temp >= 0) && (u_temp <= 1)) {
                                rec.u = u_temp;
                            }
                            else {
                                rec.u = (-Bu - sqrt(Bu*Bu-4*Au*Cu)) / (2*Au);
                            }
                        }


                        if (Av == 0) {
                            rec.v = -Cv/Bv;
                        }
                        else {
                            float v_temp = (-Bv + sqrt(Bv*Bv-4*Av*Cv)) / (2*Av);
                            if ((v_temp >= 0) && (v_temp <= 1)) {
                                rec.v = v_temp;
                            }
                            else {
                                rec.v = (-Bv - sqrt(Bv*Bv-4*Av*Cv)) / (2*Av);
                            }
                        }
                    }
                    return true;
                }
                 return false;
            }
            return false;
        }
}
