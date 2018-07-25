#ifndef QUADRATIC_ELLIPSOID_H
#define QUADRATIC_ELLIPSOID_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class quadratic_ellipsoid : public hitable
{
    public:
        quadratic_ellipsoid() {}
        quadratic_ellipsoid(vec3 cen, float a, float b, float c, material *m) : center(cen), intercept_x(a), intercept_y(b), intercept_z(c), ma(m) {}
/*
(x-xc)^2/a^2 + (y-yc)^2/b^2 + (z-zc)^2/c^2 = 1
*/
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float intercept_x;
        float intercept_y;
        float intercept_z;
        material *ma;
};

#endif // QUADRATIC_ELLIPSOID_H
