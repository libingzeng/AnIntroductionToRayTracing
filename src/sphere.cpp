#include "sphere.h"

#include <iostream>
using namespace std;

bool sphere::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
#if SPHERE_LOG == 1
        std::cout << "-------------sphere::hit----------------" << endl;
#endif // SPHERE_LOG
        vec3 oc = r.origin() - center;
        float a = dot(r.direction(), r.direction());
        float b = 2.0 * dot(oc, r.direction());
        float c = dot(oc, oc) - radius*radius;
        float discriminant = b*b - 4*a*c;

        if (discriminant > 0) {
            if (csg) {
                rec.t = (-b - sqrt(discriminant)) / (2.0*a);
                rec.p = r.point_at_parameter(rec.t);
                rec.normal = unit_vector((rec.p - center) / radius);
                rec.t2 = (-b + sqrt(discriminant)) / (2.0*a);
                rec.p2 = r.point_at_parameter(rec.t2);
                rec.normal2 = unit_vector((rec.p2 - center) / radius);
                rec.mat_ptr = ma;
                rec.u = -1.0;
                rec.v = -1.0;
                return true;
            }
            else {
                float temp = (-b - sqrt(discriminant)) / (2.0*a);
                if (temp < t_max && temp > t_min) {
                    rec.t = temp;
                    rec.p = r.point_at_parameter(rec.t);
                    rec.normal = unit_vector((rec.p - center) / radius);
                    rec.mat_ptr = ma;
                    if (inverse) {
                        vec3 pole = vec3(0, 1, 0);
                        vec3 equator = vec3(0, 0, 1);
                        float u, v;
                        float phi = acos(-dot(rec.normal, pole));
                        v = phi / M_PI;
                        float theta = acos((dot(equator, rec.normal)) / sin(phi)) / (2*M_PI);
                        if (dot(cross(pole, equator), rec.normal) > 0) {
                            u = theta;
                        }
                        else {
                            u = 1 - theta;
                        }
                        rec.u = u;
                        rec.v = v;
                    }
                    else {
                        rec.u = -1.0;
                        rec.v = -1.0;
                    }

    //                rec.c = center;
    //                rec.r = radius;

    //                std::cout << "-------------sphere::hit---1-------------" << endl;
                    return true;
                }
                temp = (-b + sqrt(discriminant)) / (2.0*a);
                if (temp < t_max && temp > t_min) {
                    rec.t = temp;
                    rec.p = r.point_at_parameter(rec.t);
                    rec.normal = unit_vector((rec.p - center) / radius);
                    rec.mat_ptr = ma;

                    vec3 pole = vec3(0, 1, 0);
                    vec3 equator = vec3(0, 0, 1);
                    float u, v;
                    float phi = acos(-dot(rec.normal, pole));
                    v = phi / M_PI;
                    float theta = acos((dot(equator, rec.normal)) / sin(phi)) / (2*M_PI);
                    if (dot(cross(pole, equator), rec.normal) > 0) {
                        u = theta;
                    }
                    else {
                        u = 1 - theta;
                    }
                    rec.u = u;
                    rec.v = v;

    //                rec.c = center;
    //                rec.r = radius;
    //                std::cout << "-------------sphere::hit---2-------------" << endl;
                    return true;
                }
            }
        }
//        std::cout << "-------------sphere::hit---3-------------" << endl;
        return false;
}/*
sphere::sphere()
{
    //ctor
}

sphere::~sphere()
{
    //dtor
}
*/
