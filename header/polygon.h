#ifndef POLYGON_H
#define POLYGON_H

#include "hitable.h"


class polygon : public hitable
{
    public:
        polygon() {}
        polygon(vec3 *v, int n, material *m) : vertexes(v), number(n), ma(m) {}
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 *vertexes;
        int number;
        material *ma;
};

#endif // POLYGON_H
