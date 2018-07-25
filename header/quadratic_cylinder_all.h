#ifndef QUADRATIC_CYLINDER_ALL_H
#define QUADRATIC_CYLINDER_ALL_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class quadratic_cylinder_all : public hitable
{
    public:
        quadratic_cylinder_all() {}
        quadratic_cylinder_all(vec3 cen, float a, float b, float c, float hy, material *m, vec3 u, float an) {
            center = cen;
            intercept_x = a;
            intercept_y = b;
            intercept_z = c;
            height_half_y = hy;
            ma = m;
            vector_u = u;
            get_vector_vw(vector_u, an, vector_v, vector_w);
        }
/*
(x-xc)^2/a^2 + 0*(y-yc)^2/b^2 + (z-zc)^2/c^2 = 1
*/
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float intercept_x;
        float intercept_y;
        float intercept_z;
        float height_half_y;
        material *ma;
        vec3 vector_u;
        vec3 vector_v;
        vec3 vector_w;
};

#endif // QUADRATIC_CYLINDER_ALL_H
