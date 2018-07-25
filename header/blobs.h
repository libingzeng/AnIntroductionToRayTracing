#ifndef BLOBS_H
#define BLOBS_H

#include "hitable.h"
#include "material.h"
#include "log.h"

class blobs : public hitable
{
    public:
        blobs() {}
        blobs(vec3 cen1, vec3 cen2, float b1, float r1, float b2, float r2, float s, int in, float tol, material *m) :
            center1(cen1), center2(cen2), blob_p1(b1), radius1(r1), blob_p2(b2), radius2(r2), sum(s),
            initial_number(in), tolerance(tol), ma(m) {}
/*
f(x,y,z)=   exp((B1/(R1^2) * ((x-x1)^2+(y-y1)^2+(z-z1)^2) - B1)
          + exp((B2/(R2^2) * ((x-x2)^2+(y-y2)^2+(z-z2)^2) - B2) - s = 0
NOTE: in our program, x1=x2, z1=z2
s: should be bigger than 1
in: initial number
tol: tolerance
*/
        virtual bool hit(const ray& r, float tmin, float tmax, hit_record& rec) const;
        vec3 center1, center2;
        float blob_p1, radius1, blob_p2, radius2, sum;
        int initial_number;
        float tolerance;
        material *ma;
};

#endif // BLOBS_H
