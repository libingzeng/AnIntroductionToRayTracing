#ifndef MATERIAL_H
#define MATERIAL_H

#include "hitable.h"

class material
{
    public:
        virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const = 0;
};

#endif // MATERIAL_H
