#include <iostream>
#include <limits>
#include "float.h"

#include "box.h"
#include "log.h"

using namespace std;

bool box::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
        float t_near = (numeric_limits<float>::min)();
        float t_far = (numeric_limits<float>::max)();
        int near_flag, far_flag;
        vec3 direction = r.direction();
        vec3 origin = r.origin();
        vec3 bl = vertex_l;
        vec3 bh = vertex_h;
        float array1[6];

        if(direction.x() == 0) {
            if((origin.x() < bl.x()) || (origin.x() > bh.x())) {
#if BOX_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin X0 is not between the slabs. return false" <<endl;
#endif // BOX_LOG
                return false;
            }
            array1[0] = (numeric_limits<float>::min)();
            array1[1] = (numeric_limits<float>::max)();
        }
        if(direction.y() == 0) {
            if((origin.y() < bl.y()) || (origin.y() > bh.y())) {
#if BOX_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin Y0 is not between the slabs. return false" <<endl;
#endif // BOX_LOG
                return false;
            }
            array1[2] = (numeric_limits<float>::min)();
            array1[3] = (numeric_limits<float>::max)();
        }
        if(direction.z() == 0) {
            if((origin.z() < bl.z()) || (origin.z() > bh.z())) {
#if BOX_LOG == 1
                std::cout << "the ray is parallel to the planes and the origin Z0 is not between the slabs. return false" <<endl;
#endif // BOX_LOG
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
#if BOX_LOG == 1
            std::cout << "array1[" << i << "]:" << array1[i] <<endl;
            std::cout << "array1[" << i+1 << "]:" << array1[i+1] <<endl;
#endif // BOX_LOG
            if(array1[i] >= t_near) {t_near = array1[i]; near_flag = i;}
            if(array1[i+1] <= t_far) {t_far = array1[i+1]; far_flag = i+1;}
            if(t_near > t_far) {
#if BOX_LOG == 1
                std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_near > t_far. return false" <<endl;
#endif // BOX_LOG
                return false;
            }
            if(t_far < 0) {
#if BOX_LOG == 1
                std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_far < 0. return false" <<endl;
#endif // BOX_LOG
                return false;
            }
        }

#if BOX_LOG == 1
        std::cout << "t_near: " << t_near << "   near_flag: " << near_flag <<endl;
        std::cout << "t_far: " << t_far << "   far_flag: " << far_flag <<endl;
        std::cout << "t_near,parameters: " << origin+direction*t_near << "   t_far,parameters: " << origin+direction*t_far <<endl;
        std::cout << "pass all of the tests. return ture" <<endl;
#endif // BOX_LOG
        vec3 normals_choose[6];
        if (csg) {
            rec.t = t_near;
            rec.p = r.point_at_parameter(rec.t);
            rec.mat_ptr = ma;
            for(int j=0; j<6; j++) {
                normals_choose[j] = vec3(0,0,0);
            }
            for(int i=0; i<6; i++) {
                if(dot(normals[i], r.direction()) < 0) {
                    normals_choose[i] = normals[i];
                }
            }
            for(int k=near_flag; k<6; k++) {
                if(!vector_equ(normals_choose[k], vec3(0,0,0))) {
                    rec.normal = normals_choose[k];
                    break;
                }
            }

            rec.t2 = t_far;
            rec.p2 = r.point_at_parameter(rec.t2);
            for(int j=0; j<6; j++) {
                normals_choose[j] = vec3(0,0,0);
            }
            for(int i=0; i<6; i++) {
                if(dot(normals[i], r.direction()) > 0) {
                    normals_choose[i] = normals[i];
                }
            }
            for(int k=far_flag; k<6; k++) {
                if(!vector_equ(normals_choose[k], vec3(0,0,0))) {
                    rec.normal2 = normals_choose[k];
                    break;
                }
            }
            rec.u = -1.0;
            rec.v = -1.0;
            return true;
        }
        else {
            if (t_near < t_max && t_near > t_min) {
                rec.t = t_near;
                rec.p = r.point_at_parameter(rec.t);
                rec.mat_ptr = ma;

                for(int j=0; j<6; j++) {
                    normals_choose[j] = vec3(0,0,0);
                }
                for(int i=0; i<6; i++) {
                    if(dot(normals[i], r.direction()) < 0) {
                        normals_choose[i] = normals[i];
                    }
                }
                for(int k=near_flag; k<6; k++) {
                    if(!vector_equ(normals_choose[k], vec3(0,0,0))) {
                        rec.normal = normals_choose[k];
                        break;
                    }
                }
                return true;
            }
        }

        return false;
}
