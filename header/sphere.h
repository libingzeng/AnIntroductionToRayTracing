#ifndef SPHERE_H
#define SPHERE_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class sphere: public hitable{
    public:
        sphere() {}
        sphere(vec3 cen, float r, material *m, bool in, bool csg) : center(cen), radius(r), ma(m), inverse(in), csg(csg) {}
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center;
        float radius;
        material *ma;
        bool inverse;
        bool csg;
};
#endif // SPHERE_H
