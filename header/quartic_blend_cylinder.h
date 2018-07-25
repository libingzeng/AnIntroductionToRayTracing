#ifndef QUARTIC_BLEND_CYLINDER_H
#define QUARTIC_BLEND_CYLINDER_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class quartic_blend_cylinder : public hitable
{
    public:
        quartic_blend_cylinder() {}
        quartic_blend_cylinder(vec3 cen1, float a1, float b1, float lh1, vec3 cen2, float a2, float b2, float lh2, float a3, float b3, material *m) {
            center1 = cen1;
            intercept_x1 = a1;
            intercept_z1 = b1;
            length_half_y1 = lh1;
            center2 = cen2;
            intercept_y2 = a2;
            intercept_z2 = b2;
            length_half_x2 = lh2;
            intercept_s1 = a3;
            intercept_s2 = b3;
            ma = m;
        }

        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;

        vec3 center1;
        float intercept_x1, intercept_z1, length_half_y1;
        vec3 center2;
        float intercept_y2, intercept_z2, length_half_x2, intercept_s1, intercept_s2;
        material *ma;
};

#endif // QUARTIC_BLEND_CYLINDER_H
