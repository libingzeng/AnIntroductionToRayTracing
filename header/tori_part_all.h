#ifndef TORI_PART_ALL_H
#define TORI_PART_ALL_H


#include "hitable.h"
#include "material.h"
#include "log.h"

class tori_part_all : public hitable
{
    public:
        tori_part_all() {}
        tori_part_all(vec3 cen, float ra, float rb, material *m, float t1, float t2, vec3 u, float an) {
             center = cen;
             radius_a = ra;
             radius_b = rb;
             ma = m;
             theta1 = t1*M_PI/180;
             theta2 = t2*M_PI/180;
             vector_u = u;
             get_vector_vw(vector_u, an, vector_v, vector_w);

        }
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float radius_a;
        float radius_b;
        material *ma;
        float theta1, theta2;
        vec3 vector_u;
        vec3 vector_v;
        vec3 vector_w;
};

#endif // TORI_PART_ALL_H
