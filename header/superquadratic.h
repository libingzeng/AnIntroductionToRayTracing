#ifndef SUPERQUADRATIC_H
#define SUPERQUADRATIC_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class superquadratic : public hitable
{
    public:
        superquadratic() {}
        superquadratic(vec3 cen, float a1, float a2, float a3, float r, float s, float t, material *m) : center(cen), intercept_x(a1), intercept_y(a2), intercept_z(a3), p_r(r), p_s(s), p_t(t), ma(m) {}
/*
f(x,y,z)=(|x/a1|)^r + (|y/a2|)^s + (|z/a3|)^t -1 = 0
*/
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float intercept_x, intercept_y, intercept_z;
        float p_r, p_s, p_t;
        material *ma;
};

#endif // SUPERQUADRATIC_H
