#ifndef SUPERHYPERBOLOID_H
#define SUPERHYPERBOLOID_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class superhyperboloid : public hitable
{
    public:
        superhyperboloid() {}
        superhyperboloid(vec3 cen, float a1, float a2, float a3, float e1, float e2, float s1, float s2, float hy, int in, float tol,  material *m) :
            center(cen), intercept_x(a1), intercept_y(a2), intercept_z(a3), p_e1(e1), p_e2(e2),
            sign1(s1), sign2(s2), half_y(hy), initial_number(in), tolerance(tol), ma(m) {}
/*
f(x,y,z)=( (x/a1)^(2/e2) + s2*(z/a3)^(2/e2) )^(e2/e1) + s1*(y/a2)^(2/e1) -1 = 0
in: initial number
tol: tolerance
s1,s2: 1, 1: superellipsoid
s1,s2:-1, 1: superhyperboloids of one sheet
s1,s2:-1,-1: superhyperboloids of two sheets
hy: half height of the surface in y-direction
NOTE: there are something wrong with two sheets situation.
*/
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float intercept_x, intercept_y, intercept_z;
        float p_e1, p_e2;
        float sign1, sign2;
        float half_y;
        int initial_number;
        float tolerance;
        material *ma;
};

#endif // SUPERHYPERBOLOID_H
