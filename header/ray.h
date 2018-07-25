#ifndef RAY_H
#define RAY_H
#include "vec3.h"

class ray
{
    public:
        ray() {}
        ray(const vec3& a, const vec3& b) { A = a; B = b; }
        vec3 origin() const     { return A; }
        vec3 direction() const { return B; }
        vec3 point_at_parameter(float t) const { return A + t*B; }

        vec3 A;
        vec3 B;

/*
        virtual ~ray();

    protected:

    private:
*/
};

bool ray_hit_box_general(const ray& r, const vec3& vertex_l, const vec3& vertex_h, float& t_near, float& t_far);

#endif // RAY_H
