#include "lambertian.h"

#include <iostream>
using namespace std;

vec3 random_in_unit_sphere() {
    vec3 p;
    do {
        p = 2.0*vec3((rand()%(100)/(float)(100)),
                    (rand()%(100)/(float)(100)),
                    (rand()%(100)/(float)(100)))
            - vec3(1,1,1);
    } while (p.squared_length() >= 1.0);
    return p;
}

vec3 get_albedo(float u, float v) {

    if ((v >= 0) && (v < 0.5)) {
        if ((u >= 0) && (u < 0.25)) {
            return vec3(0, 0, 0);
        }
        else if ((u >= 0.25) && (u < 0.5)) {
            return vec3(0, 0, 1);
        }
        else if ((u >= 0.5) && (u < 0.75)) {
            return vec3(0, 1, 0);
        }
        else {
            return vec3(0, 1, 1);
        }
    }
    else if ((v >= 0.5) && (v <= 1.0)) {
        if ((u >= 0) && (u < 0.25)) {
            return vec3(1, 0, 0);
        }
        else if ((u >= 0.25) && (u < 0.5)) {
            return vec3(1, 0, 1);
        }
        else if ((u >= 0.5) && (u < 0.75)) {
            return vec3(1, 1, 0);
        }
        else if ((u >= 0.75) && (u <= 1.0)) {
            return vec3(1, 1, 1);
        }
    }
    return vec3(-1, -1, -1);//xcode
}

bool lambertian::scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
//    std::cout << "-------------lambertian::scatter----------------" << endl;
    vec3 target = rec.p + rec.normal + random_in_unit_sphere();
    scattered = ray(rec.p, target-rec.p);
    if ((rec.u >= 0.0) && (rec.u <= 1.0) && (rec.v >= 0.0) && (rec.v <= 1.0)) {
        attenuation = get_albedo(rec.u, rec.v);
    }
    else {
        attenuation = albedo;
//        attenuation = get_albedo(rec.u, rec.v);
    }
    return true;
}

/*
lambertian::lambertian()
{
    //ctor
}

lambertian::~lambertian()
{
    //dtor
}
*/
