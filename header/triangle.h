#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "hitable.h"


class triangle : public hitable
{
    public:
        triangle() {}
        triangle(vec3 v1, vec3 v2, vec3 v3, material *m) : vertex1(v1), vertex2(v2), vertex3(v3), ma(m) {}
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 vertex1, vertex2, vertex3;
        material *ma;
};

#endif // TRIANGLE_H
