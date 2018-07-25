#ifndef ELLIPTIC_PLANE_H
#define ELLIPTIC_PLANE_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class elliptic_plane : public hitable
{
    public:
        elliptic_plane() {}
        elliptic_plane(vec3 cen, vec3 n, float a, float b, float c, material *m) : center(cen), normal(n), intercept_x(a), intercept_y(b), intercept_z(c), ma(m) {}
/*
normal: the normal of the plane;
center: the center of the ellipse or circle;
intercept_x/y/z: the intercept of the ellipse or circle on x/y/z axis
(for circle, a=b=c=radius of the circle.)
*/
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center, normal;
        float intercept_x, intercept_y, intercept_z;
        material *ma;
};

#endif // ELLIPTIC_PLANE_H
