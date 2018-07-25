#include "triangle.h"
#include "vec2.h"
#include "log.h"

#include <iostream>
#include <fstream>

using namespace std;

bool triangle::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if TRIANGLE_LOG == 1
        ofstream outfile( ".\\results\\log_triangle.txt", ios_base::app);
        outfile << "-------------triangle::hit----------------" << endl;
#endif // TRIANGLE_LOG
        vec3 tri_n = unit_vector(cross((vertex1-vertex2), (vertex2-vertex3)));
        float tri_d = -(dot(tri_n, vertex1));
        float vd = dot(tri_n, r.direction());
        float v0 = -(dot(tri_n, r.origin()) + tri_d);

        if(vd == 0) {
#if TRIANGLE_LOG == 1
            outfile << "-------------triangle::hit-------1---------" << endl;
#endif // TRIANGLE_LOG
            return false;
        }
        else {
#if TRIANGLE_LOG == 1
            outfile << "-------------triangle::hit-------2---------" << endl;
#endif // TRIANGLE_LOG
            float t = v0/vd;
            vec3 pi = r.point_at_parameter(t);

            float temp = abs(tri_n.x());
            int i = 1;
            if(temp <= abs(tri_n.y())) {
                temp = abs(tri_n.y());
                i++;
            }
            if(temp <= abs(tri_n.z())) {
                i++;
            }
#if TRIANGLE_LOG == 1
            outfile << "-------------triangle::hit-------2.2---------" << endl;
            outfile << "triangle normal: " << tri_n << endl;
            outfile << "vertex1: " << vertex1 << endl;
            outfile << "vertex2: " << vertex2 << endl;
            outfile << "vertex3: " << vertex3 << endl;
            outfile << "pi: " << pi << endl;

#endif // TRIANGLE_LOG

            vec2 vertexs[4];
            switch (i) {
            case 1:
                vertexs[0] = vec2(vertex1.y(),vertex1.z());
                vertexs[1] = vec2(vertex2.y(),vertex2.z());
                vertexs[2] = vec2(vertex3.y(),vertex3.z());
                vertexs[3] = vec2(pi.y(),pi.z());
                break;
            case 2:
                vertexs[0] = vec2(vertex1.x(),vertex1.z());
                vertexs[1] = vec2(vertex2.x(),vertex2.z());
                vertexs[2] = vec2(vertex3.x(),vertex3.z());
                vertexs[3] = vec2(pi.x(),pi.z());
                break;
           case 3:
                vertexs[0] = vec2(vertex1.x(),vertex1.y());
                vertexs[1] = vec2(vertex2.x(),vertex2.y());
                vertexs[2] = vec2(vertex3.x(),vertex3.y());
                vertexs[3] = vec2(pi.x(),pi.y());
                break;
            }
#if TRIANGLE_LOG == 1
            outfile << "-------------triangle::hit-------2.2.1---------" << endl;
            outfile << "vertexs[0]: " << vertexs[0] << endl;
            outfile << "vertexs[1]: " << vertexs[1] << endl;
            outfile << "vertexs[2]: " << vertexs[2] << endl;
            outfile << "vertexs[3]: " << vertexs[3] << endl;

#endif // TRIANGLE_LOG
            vertexs[0] = vertexs[0] - vertexs[3];
            vertexs[1] = vertexs[1] - vertexs[3];
            vertexs[2] = vertexs[2] - vertexs[3];
            vertexs[3] = vertexs[0];
#if TRIANGLE_LOG == 1
            outfile << "-------------triangle::hit-------2.3---------" << endl;
            outfile << "vertexs[0]:" << vertexs[0].u() << " " << vertexs[0].v() << endl;
            outfile << "vertexs[1]:" << vertexs[1].u() << " " << vertexs[1].v() << endl;
            outfile << "vertexs[2]:" << vertexs[2].u() << " " << vertexs[2].v() << endl;
            outfile << "vertexs[3]:" << vertexs[3].u() << " " << vertexs[3].v() << endl;
#endif // TRIANGLE_LOG

            int sh, nsh;
            int nc = 0;
            if(vertexs[0].v() < 0) { sh = -1;}
            else { sh = 1;}
            for(int j=0; j<3; j++) {
                if(vertexs[j+1].v() < 0) { nsh = -1;}
                else { nsh = 1;}
                if(sh != nsh) {
                    if((vertexs[j].u() > 0) && (vertexs[j+1].u() >0)) { nc++;}
                    else
                        if((vertexs[j].u() > 0) || (vertexs[j+1].u() >0)) {
                            if(vertexs[j].u() - (vertexs[j].v())*(vertexs[j+1].u()-vertexs[j].u())/(vertexs[j+1].v()-vertexs[j].v()) > 0) { nc++;}
                        }
                }
                sh = nsh;
            }
#if TRIANGLE_LOG == 1
            outfile << "-------------triangle::hit-------2.5---------nc:" << nc << endl;
#endif // TRIANGLE_LOG
            if((nc)%(2)) {
                if (t < t_max && t > t_min) {
                    rec.t = t;
                    rec.p = r.point_at_parameter(rec.t);
                    rec.normal = tri_n;
                    rec.mat_ptr = ma;
#if TRIANGLE_LOG == 1
                    outfile << "-------------triangle::hit-------3---------" << endl;
#endif // TRIANGLE_LOG
                    return true;
                }
#if TRIANGLE_LOG == 1
                    outfile << "-------------triangle::hit-------4---------" << endl;
#endif // TRIANGLE_LOG
                 return false;
            }
#if TRIANGLE_LOG == 1
            outfile << "-------------triangle::hit-------5---------" << endl;
#endif // TRIANGLE_LOG
            return false;
        }
}
