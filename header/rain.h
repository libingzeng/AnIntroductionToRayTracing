#ifndef RAIN_H
#define RAIN_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class rain : public hitable
{
    public:
        rain() {}
        rain(vec3 cen, float a1, float a2, float a3, material *m) : center(cen), x_a1(a1), y_a2(a2), z_a3(a3), ma(m) {}
/*
parametric equation:
    x=a1*cos(u)*sin(v)*(1-cos(v))/2
    y=a2*cos(v)
    z=a3*sin(u)*sin(v)*(1-cos(v))/2
imlipcit equation:
    4*(x/a1)^2+4*(z/a3)^2+(y/a2)^4-2*(y/a2)^3+2*(y/a2)-1=0
*/
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float x_a1, y_a2, z_a3;
        material *ma;
};
#endif // RAIN_H
