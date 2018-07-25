#include "ray.h"

#include <iostream>
#include <limits>
#include "float.h"

bool ray_hit_box_general(const ray& r, const vec3& vertex_l, const vec3& vertex_h, float& t_near, float& t_far) {
        t_near = (numeric_limits<float>::min)();
        t_far = (numeric_limits<float>::max)();
        vec3 direction = r.direction();
        vec3 origin = r.origin();
        vec3 bl = vertex_l;
        vec3 bh = vertex_h;
        float array1[6];

        if(direction.x() == 0) {
            if((origin.x() < bl.x()) || (origin.x() > bh.x())) {
                return false;
            }
            array1[0] = (numeric_limits<float>::min)();
            array1[1] = (numeric_limits<float>::max)();
        }
        if(direction.y() == 0) {
            if((origin.y() < bl.y()) || (origin.y() > bh.y())) {
                return false;
            }
            array1[2] = (numeric_limits<float>::min)();
            array1[3] = (numeric_limits<float>::max)();
        }
        if(direction.z() == 0) {
            if((origin.z() < bl.z()) || (origin.z() > bh.z())) {
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
            if(array1[i] >= t_near) {t_near = array1[i];}
            if(array1[i+1] <= t_far) {t_far = array1[i+1];}
            if(t_near > t_far) {
                return false;
            }
            if(t_far < 0) {
                return false;
            }
        }
        if (t_near != t_near) {
            t_near = t_near * 1;
        }
        return true;
}
