#ifndef TORI_H
#define TORI_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class tori : public hitable
{
    public:
        tori() {}
        tori(vec3 cen, float ra, float rb, material *m) : center(cen), radius_a(ra), radius_b(rb), ma(m) {}
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float radius_a;
        float radius_b;
        material *ma;
};

#endif // TORI_H
