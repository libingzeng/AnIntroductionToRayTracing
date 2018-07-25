#ifndef TORI_PART_H
#define TORI_PART_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class tori_part : public hitable
{
    public:
        tori_part() {}
        tori_part(vec3 cen, float ra, float rb, material *m, float t1, float t2) {
             center = cen;
             radius_a = ra;
             radius_b = rb;
             ma = m;
             theta1 = t1*M_PI/180;
             theta2 = t2*M_PI/180;
        }
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float radius_a;
        float radius_b;
        material *ma;
        float theta1, theta2;
};

#endif // TORI_PART_H
