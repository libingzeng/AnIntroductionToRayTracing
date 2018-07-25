#ifndef HITABLE_H
#define HITABLE_H

#include "ray.h"

class material;

struct hit_record{
    float t, u, v;
    vec3 p;
    vec3 normal;
    material *mat_ptr;
    float t2;// for csg
    vec3 p2;// for csg
    vec3 normal2;// for csg
};

class hitable
{
    public:
        virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec) const = 0;
};
/*
"hitable" is an "abstract class" for anything a ray might hit.
Ҳ���ǿ��ܱ�����ײ�ϵ��κζ��������ණ���Ĺ��Ծ��ǡ�������hit����
*/
#endif // HITABLE_H
