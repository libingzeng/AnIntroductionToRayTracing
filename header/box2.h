#ifndef BOX2_H
#define BOX2_H

#include "hitable.h"

class box2 : public hitable
{
    public:
        box2() {}
        box2(vec3 u, float an, float a, float b, float c, vec3 p, material *m) {

            vector_u = unit_vector(u);
            vector_v = unit_vector(get_vector_v(vector_u, an));
            vector_w = unit_vector(cross(vector_u, vector_v));

            vertex_l = vector_trans(p, vector_u, vector_v, vector_w);
            vertex_h = vector_trans((p + a*vector_u + c*vector_v - b*vector_w), vector_u, vector_v, vector_w);

            normals[0] = vector_trans_back(vec3(-1, 0, 0), vector_u, vector_v, vector_w);//left
            normals[1] = vector_trans_back(vec3(1, 0, 0), vector_u, vector_v, vector_w);//right
            normals[2] = vector_trans_back(vec3(0, 1, 0), vector_u, vector_v, vector_w);;//up
            normals[3] = vector_trans_back(vec3(0, -1, 0), vector_u, vector_v, vector_w);;//down
            normals[4] = vector_trans_back(vec3(0, 0, 1), vector_u, vector_v, vector_w);;//front
            normals[5] = vector_trans_back(vec3(0, 0, -1), vector_u, vector_v, vector_w);;//back

            ma = m;
        }
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 vector_u;
        vec3 vector_v;
        vec3 vector_w;
        vec3 vertex_l;
        vec3 vertex_h;
        vec3 normals[6];
        material *ma;
};

#endif // BOX2_H
