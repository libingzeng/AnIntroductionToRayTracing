#include "metal.h"
#include "lambertian.h"

#include <iostream>
using namespace std;

vec3 reflect(const vec3& v, const vec3& n) {
    return v - 2*dot(v,n)*n;
}

bool metal::scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
//    std::cout << "-------------metal::scatter----------------" << endl;
    vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
    scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere());
    attenuation = albedo;
    return true;
//    return (dot(scattered.direction(), rec.normal) > 0);
}

/*
metal::metal()
{
    //ctor
}

metal::~metal()
{
    //dtor
}
*/
