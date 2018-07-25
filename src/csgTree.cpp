#include "csgTree.h"
extern vec3 lookfrom;

bool csgDifference(const ray& r, bool hit_sphere, bool hit_box, hit_record rec_sphere, hit_record rec_box, float t_min, hit_record& rec) {
        float t[4], temp_t;
        int num = 0;
        vec3 normal[4], temp_normal;
        material *mat_ptr[4], *temp_mat_ptr;
        if (hit_sphere) {
            if (hit_box) {
                if (rec_sphere.t > t_min) {
                    t[num] = rec_sphere.t;
                    normal[num] = rec_sphere.normal;
                    mat_ptr[num] = rec_sphere.mat_ptr;
                    num ++;
                }
                if (rec_sphere.t2 > t_min) {
                    t[num] = rec_sphere.t2;
                    normal[num] = rec_sphere.normal2;
                    mat_ptr[num] = rec_sphere.mat_ptr;
                    num ++;
                }
                if (rec_box.t > t_min) {
                    t[num] = rec_box.t;
                    normal[num] = rec_box.normal;
                    mat_ptr[num] = rec_box.mat_ptr;
                    num ++;
                }
                if (rec_box.t2 > t_min) {
                    t[num] = rec_box.t2;
                    normal[num] = rec_box.normal2;
                    mat_ptr[num] = rec_box.mat_ptr;
                    num ++;
                }
                for (int i=0; i<(num-1); i++) {
                    for (int j=i+1; j<num; j++) {
                        if (t[i] > t[j]) {
                            temp_t = t[i];
                            t[i] = t[j];
                            t[j] = temp_t;
                            temp_normal = normal[i];
                            normal[i] = normal[j];
                            normal[j] = temp_normal;
                            temp_mat_ptr = mat_ptr[i];
                            mat_ptr[i] = mat_ptr[j];
                            mat_ptr[j] = temp_mat_ptr;
                        }
                    }
                }
                if (fabs(t[0]-rec_box.t)<1e-6) {
                    if (fabs(t[3]-rec_box.t2)<1e-6) {
                        return false;
                    }
                    else {
                        rec.t = rec_box.t2;
                        rec.p = r.point_at_parameter(rec.t);
                        rec.normal = rec_box.normal2;
                        if(dot(r.direction(), rec.normal) > 0) {
                            rec.normal = - rec.normal;
                        }
                        rec.mat_ptr = mat_ptr[0];
                        rec.u = -1.0;
                        rec.v = -1.0;

                        rec.t2 = t[3];
                        rec.p2 = r.point_at_parameter(rec.t2);
                        rec.normal2 = normal[3];
                        return true;
                    }
                }
                else {
                    rec.t = t[0];
                    rec.p = r.point_at_parameter(rec.t);
                    rec.normal = normal[0];
                    if(dot(r.direction(), rec.normal) > 0) {
                        rec.normal = - rec.normal;
                    }
                    rec.mat_ptr = mat_ptr[0];
                    rec.u = -1.0;
                    rec.v = -1.0;

// the interval between rec_box.t and rec_box.t2 is out of the set.
// the valid interval should be [t[0], rec_box.t] U [rec_box.t2, t[3]].
// but, here, we store [t[0], t[3]] as the interval.
                    rec.t2 = t[3];
                    rec.p2 = r.point_at_parameter(rec.t2);
                    rec.normal2 = normal[3];
                    return true;
                }
            }
            else if (!hit_box) {
                if (rec_sphere.t > t_min) {
                    t[num] = rec_sphere.t;
                    normal[num] = rec_sphere.normal;
                    mat_ptr[num] = rec_sphere.mat_ptr;
                    num ++;
                }
                if (rec_sphere.t2 > t_min) {
                    t[num] = rec_sphere.t2;
                    normal[num] = rec_sphere.normal2;
                    mat_ptr[num] = rec_sphere.mat_ptr;
                    num ++;
                }
                for (int i=0; i<(num-1); i++) {
                    for (int j=i+1; j<num; j++) {
                        if (t[i] > t[j]) {
                            temp_t = t[i];
                            t[i] = t[j];
                            t[j] = temp_t;
                            temp_normal = normal[i];
                            normal[i] = normal[j];
                            normal[j] = temp_normal;
                            temp_mat_ptr = mat_ptr[i];
                            mat_ptr[i] = mat_ptr[j];
                            mat_ptr[j] = temp_mat_ptr;
                        }
                    }
                }
                if (t[0] > t_min) {
                    rec.t = t[0];
                    rec.p = r.point_at_parameter(rec.t);
                    rec.normal = normal[0];
                    if(dot(r.direction(), rec.normal) > 0) {
                        rec.normal = - rec.normal;
                    }
                    rec.mat_ptr = mat_ptr[0];
                    rec.u = -1.0;
                    rec.v = -1.0;

                    rec.t2 = t[1];
                    rec.p2 = r.point_at_parameter(rec.t2);
                    rec.normal2 = normal[1];
                    return true;
                }
            }
        }
        return false;
}

bool csgUnion(const ray& r, bool hit_sphere, bool hit_box, hit_record rec_sphere, hit_record rec_box, float t_min, hit_record& rec) {
        float t[4], temp_t;
        int num = 0;
        vec3 normal[4], temp_normal;
        material *mat_ptr[4], *temp_mat_ptr;
        if (hit_sphere || hit_box) {
            if (hit_sphere && hit_box) {
                if (rec_sphere.t > t_min) {
                    t[num] = rec_sphere.t;
                    normal[num] = rec_sphere.normal;
                    mat_ptr[num] = rec_sphere.mat_ptr;
                    num ++;
                }
                if (rec_sphere.t2 > t_min) {
                    t[num] = rec_sphere.t2;
                    normal[num] = rec_sphere.normal2;
                    mat_ptr[num] = rec_sphere.mat_ptr;
                    num ++;
                }
                if (rec_box.t > t_min) {
                    t[num] = rec_box.t;
                    normal[num] = rec_box.normal;
                    mat_ptr[num] = rec_box.mat_ptr;
                    num ++;
                }
                if (rec_box.t2 > t_min) {
                    t[num] = rec_box.t2;
                    normal[num] = rec_box.normal2;
                    mat_ptr[num] = rec_box.mat_ptr;
                    num ++;
                }
            }
            else if (hit_sphere && !hit_box) {
                if (rec_sphere.t > t_min) {
                    t[num] = rec_sphere.t;
                    normal[num] = rec_sphere.normal;
                    mat_ptr[num] = rec_sphere.mat_ptr;
                    num ++;
                }
                if (rec_sphere.t2 > t_min) {
                    t[num] = rec_sphere.t2;
                    normal[num] = rec_sphere.normal2;
                    mat_ptr[num] = rec_sphere.mat_ptr;
                    num ++;
                }
            }
            else if (!hit_sphere && hit_box) {
                if (rec_box.t > t_min) {
                    t[num] = rec_box.t;
                    normal[num] = rec_box.normal;
                    mat_ptr[num] = rec_box.mat_ptr;
                    num ++;
                }
                if (rec_box.t2 > t_min) {
                    t[num] = rec_box.t2;
                    normal[num] = rec_box.normal2;
                    mat_ptr[num] = rec_box.mat_ptr;
                    num ++;
                }
            }
            for (int i=0; i<(num-1); i++) {
                for (int j=i+1; j<num; j++) {
                    if (t[i] > t[j]) {
                        temp_t = t[i];
                        t[i] = t[j];
                        t[j] = temp_t;
                        temp_normal = normal[i];
                        normal[i] = normal[j];
                        normal[j] = temp_normal;
                        temp_mat_ptr = mat_ptr[i];
                        mat_ptr[i] = mat_ptr[j];
                        mat_ptr[j] = temp_mat_ptr;
                    }
                }
            }
            if (t[0] > t_min) {
                rec.t = t[0];
                rec.p = r.point_at_parameter(rec.t);
                rec.normal = normal[0];
                if(dot(r.direction(), rec.normal) > 0) {
                    rec.normal = - rec.normal;
                }
                rec.mat_ptr = mat_ptr[0];
                rec.u = -1.0;
                rec.v = -1.0;

                rec.t2 = t[num-1];
                rec.p2 = r.point_at_parameter(rec.t2);
                rec.normal2 = normal[num-1];
                return true;
            }
        }
        return false;
}

bool csgIntersection(const ray& r, bool hit_sphere, bool hit_box, hit_record rec_sphere, hit_record rec_box, float t_min, hit_record& rec) {
        float t1[2], t2[2], t[4], temp_t;
        vec3 normal[4], temp_normal;
        material *mat_ptr[4], *temp_mat_ptr;
        if (hit_sphere && hit_box) {
            if (rec_sphere.t < rec_sphere.t2) {
                t1[0] = rec_sphere.t;
                t1[1] = rec_sphere.t2;
            }
            else {
                t1[0] = rec_sphere.t2;
                t1[1] = rec_sphere.t;
            }
            if (rec_box.t < rec_box.t2) {
                t2[0] = rec_box.t;
                t2[1] = rec_box.t2;
            }
            else {
                t2[0] = rec_box.t2;
                t2[1] = rec_box.t;
            }
            if ((t1[1]>t2[0]) && (t1[0]<t2[1])) {
                t[0] = rec_sphere.t;
                normal[0] = rec_sphere.normal;
                mat_ptr[0] = rec_sphere.mat_ptr;
                t[1] = rec_sphere.t2;
                normal[1] = rec_sphere.normal2;
                mat_ptr[1] = rec_sphere.mat_ptr;
                t[2] = rec_box.t;
                normal[2] = rec_box.normal;
                mat_ptr[2] = rec_box.mat_ptr;
                t[3] = rec_box.t2;
                normal[3] = rec_box.normal2;
                mat_ptr[3] = rec_box.mat_ptr;
                for (int i=0; i<3; i++) {
                    for (int j=i+1; j<4; j++) {
                        if (t[i] > t[j]) {
                            temp_t = t[i];
                            t[i] = t[j];
                            t[j] = temp_t;
                            temp_normal = normal[i];
                            normal[i] = normal[j];
                            normal[j] = temp_normal;
                            temp_mat_ptr = mat_ptr[i];
                            mat_ptr[i] = mat_ptr[j];
                            mat_ptr[j] = temp_mat_ptr;
                        }
                    }
                }
                if (t[1] > t_min) {
                    rec.t = t[1];
                    rec.p = r.point_at_parameter(rec.t);
                    rec.normal = normal[1];
                    if(dot(r.direction(), rec.normal) > 0) {
                        rec.normal = - rec.normal;
                    }
                    rec.mat_ptr = mat_ptr[1];
                    rec.u = -1.0;
                    rec.v = -1.0;

                    rec.t2 = t[2];
                    rec.p2 = r.point_at_parameter(rec.t2);
                    rec.normal2 = normal[2];
                    return true;
                }
            }
        }
        return false;
}

bool csgTree::hit(const ray& r, float t_min, float t_max, hit_record& rec) const {
        if (!vector_equ(r.origin(), lookfrom)) {
        // this is a bad trick for avoiding reflect or refract rays from csg hit itself.
            return false;
        }
        hit_record rec_sphere, rec_box;
        bool hit_sphere = false;
        bool hit_box = false;
        hit_sphere = root->left->info.solid->hit(r, t_min, t_max, rec_sphere);
        hit_box = root->right->info.solid->hit(r, t_min, t_max, rec_box);
        if (root->info.operation == 1) {//
            return (csgUnion(r, hit_sphere, hit_box, rec_sphere, rec_box, t_min, rec));
        }
        if (root->info.operation == 2) {
            return (csgIntersection(r, hit_sphere, hit_box, rec_sphere, rec_box, t_min, rec));
        }
        if (root->info.operation == 3) {
            return (csgDifference(r, hit_sphere, hit_box, rec_sphere, rec_box, t_min, rec));
        }
        if (root->info.operation == 4) {
            return (csgDifference(r, hit_box, hit_sphere, rec_box, rec_sphere, t_min, rec));
        }
        return false;
}
