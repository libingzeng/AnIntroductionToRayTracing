#ifndef QUADRATIC_PARABOLOID_H
#define QUADRATIC_PARABOLOID_H

#include "hitable.h"


class quadratic_paraboloid : public hitable
{
    public:
        quadratic_paraboloid() {}
        quadratic_paraboloid(vec3 cen, float p, float q, float hy, material *m) : center(cen), focus_directrix_p(p), focus_directrix_q(q), height_half_y(hy), ma(m) {}
/*
(x-xc)^2/2*p + (z-zc)^2/2*q = (y-yc)
(p>0, q>0)
*/
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float focus_directrix_p;
        float focus_directrix_q;
        float height_half_y;
        material *ma;
};

#endif // QUADRATIC_PARABOLOID_H
