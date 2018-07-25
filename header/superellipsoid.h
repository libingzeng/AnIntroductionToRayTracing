#ifndef SUPERELLIPSOID_H
#define SUPERELLIPSOID_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class superellipsoid : public hitable
{
    public:
        superellipsoid() {}
        superellipsoid(vec3 cen, float a1, float a2, float a3, float e1, float e2, int in, float tol,  material *m) :
            center(cen), intercept_x(a1), intercept_y(a2), intercept_z(a3), p_e1(e1), p_e2(e2),
            initial_number(in), tolerance(tol), ma(m) {}
/*
f(x,y,z)=( (x/a1)^(2/e2) + (z/a3)^(2/e2) )^(e2/e1) + (y/a2)^(2/e1) -1 = 0
in: initial number
tol: tolerance
*/
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float intercept_x, intercept_y, intercept_z;
        float p_e1, p_e2;
        int initial_number;
        float tolerance;
        material *ma;
};

#endif // SUPERELLIPSOID_H
