#ifndef DIELECTRIC_H
#define DIELECTRIC_H

#include "material.h"
#include "metal.h"
#include "log.h"


class dielectric : public material
{
    public:
        dielectric(float ri) : ref_idx(ri) {}
        virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const;
        float ref_idx;
};

#endif // DIELECTRIC_H
