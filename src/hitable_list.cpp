#include "hitable_list.h"

#include <iostream>
using namespace std;

bool hitable_list::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
//        std::cout << "-------------hitable_list::hit----------------" << endl;
        hit_record temp_rec;
        bool hit_anything = false;
        double closest_so_far = t_max;
        for (int i = 0; i < list_size; i++) {
            if (list[i]->hit(r, t_min, closest_so_far, temp_rec)){
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }
//        std::cout << "-------------hitable_list::hit---return----------------" << endl;
        return hit_anything;
}
/*
hitable_list::hitable_list()
{
    //ctor
}

hitable_list::~hitable_list()
{
    //dtor
}
*/
