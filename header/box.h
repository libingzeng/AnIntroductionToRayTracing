#ifndef BOX_H
#define BOX_H

#include "hitable.h"

class box : public hitable
{
    public:
        box() {}
        box(vec3 vl, vec3 vh, material *m, bool csg) : vertex_l(vl), vertex_h(vh), ma(m) ,csg(csg) {
            normals[0] = vec3(-1, 0, 0);//left
            normals[1] = vec3(1, 0, 0);//right
            normals[2] = vec3(0, 1, 0);//up
            normals[3] = vec3(0, -1, 0);//down
            normals[4] = vec3(0, 0, 1);//front
            normals[5] = vec3(0, 0, -1);//back
        }
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 vertex_l;
        vec3 vertex_h;
        vec3 normals[6];
        material *ma;
        bool csg;
};

#endif // BOX_H
