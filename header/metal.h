#ifndef METAL_H
#define METAL_H

#include "material.h"

vec3 reflect(const vec3& v, const vec3& n);

class metal : public material
{
    public:
        metal(const vec3& a, float f) : albedo(a) { if(f < 1) fuzz = f; else fuzz = 1;}
        virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const;
        vec3 albedo;
        float fuzz;
};

#endif // METAL_H
