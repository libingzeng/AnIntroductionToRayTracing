#ifndef QUADRATIC_H
#define QUADRATIC_H

#include "hitable.h"
#include "material.h"
#include "log.h"


class quadratic : public hitable
{
    public:
        quadratic() {}
        quadratic(vec3 cen, float a, float b, float c, float s1, float s2, float hy, material *m) : center(cen), intercept_x(a), intercept_y(b), intercept_z(c), sign1(s1), sign2(s2), height_half_y(hy), ma(m) {}
/*
(x-xc)^2/a^2 + s1*(y-yc)^2/b^2 + (z-zc)^2/c^2 = s2
s1=  -1,s2=  1:  hyperboloid of one sheet
s1=  -1,s2= -1:  hyperboloid of two sheets
s1=  -1,s2=  0:  elliptic cone (when a=c, it's a cone)
s1=   0,s2=  1:  elliptic cylinder (when a=c, it's a cylinder)
*/
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float intercept_x;
        float intercept_y;
        float intercept_z;
        float sign1, sign2, height_half_y;
        material *ma;
};

#endif // QUADRATIC_H
