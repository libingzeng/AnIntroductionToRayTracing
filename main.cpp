#define testNumber 46
/*
1: output the first image
2: test "int &ri£¬int& ri£¬int *&pri"
3: output the first image by using vector.
4: test the class "vec3"
5: test the function "inline vec3& vec3::operator+=(const vec3 &v)"
6: test rays
7: test rays and add a sphere
8: visualize the normals of sphere.
9: several spheres
10: abstract class
11: antialiasing
12: diffuse materials
13: metal
14: dielectric
15: positionable camera and defocus blur
16: random scene
17: ray/box intersection
18: triangle
19: polygon
20: box
21: linear equation
22: get_vector_v
23: box2()
24: any box
25: ellipsoid
26: inverse mapping
27: find the roots of quartic equation
28: tori
29: cylinder all
30: get_vector_vw()
31: tori_part_all
32: acos(z_n, u_xoz)
33: quartic_blend_cylinder
34: superellipsoid
35: rain
36: parametric surface
37: matrix
38: read file
39: translational sweeping
40: roots number of equation of 6th degree or 10th degree
41: rotational sweeping
42: sphere sweeping
43: binary tree
44: binary tree 2
45: intervals
46: CSG
*/

#if testNumber == 1 /*1: output the first image*/

    #include <iostream>
    #include <fstream>

    using namespace std;

    int main()
    {
        int nx = 200;
        int ny = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/FirstImage.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        for (int j = ny-1; j >= 0; j--)
        {
            for (int i = 0; i < nx; i++)
            {
                float r = float(i) / float(nx);
                float g = float(j) / float(ny);
                float b = 0.2;
                int ir = int (255.99*r);
                int ig = int (255.99*g);
                int ib = int (255.99*b);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
    }

#elif testNumber == 2 /*2: test "int &ri£¬int& ri£¬int *&pri"*/

    #include <iostream>

    using namespace std;

    int main()
    {
        int ival = 111;

        int *pi1 = &ival;
        cout << "int *pi1 = &ival--pi1:" << pi1 << "--*pi1:" << *pi1 << "--&pi1:" << &pi1 << endl;
        int* pi2 = &ival;
        cout << "int* pi2 = &ival--pi2:" << pi2 << "--*pi2:" << *pi2 << "--&pi2:" << &pi2  << endl;
        int * pi3 = &ival;
        cout << "int * pi3 = &ival--pi3:" << pi3 << "--*pi3:" << *pi3 << "--&pi3:" << &pi3 << endl << endl;

        int &ri1 = ival;
        cout << "int &ri1 = ival--ri1:" << ri1 << "--*ri1:(error)" << "--&ri1:" << &ri1 << endl;
        int& ri2 = ival;
        cout << "int& ri2 = ival--ri2:" << ri2 << "--*ri2:(error)" << "--&ri2:" << &ri2 << endl;
        int & ri3 = ival;
        cout << "int & ri3 = ival--ri3:" << ri3 << "--*ri3:(error)" << "--&ri3:" << &ri3 << endl << endl;

        int * const &pr1 = &ival;
        cout << "int * const &pr1 = &ival--pr1:" << pr1 << "--*pr1:" << *pr1 << "--&pr1:" << &pr1 << endl << endl;

    //    int &*rp1 = &ival;
        cout << "int &*rp1 = &ival--error)" << endl << endl;


        return 0;
    }

#elif testNumber == 3 /*3: output the first image by using vector.*/

    #include <iostream>
    #include <fstream>
    #include "vec3.h"

    using namespace std;

    int main()
    {
        int nx = 200;
        int ny = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/FirstImage2.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";
        for (int j = ny-1; j >= 0; j--)
        {
            for (int i = 0; i < nx; i++)
            {
                vec3 col(float(i) / float(nx), float(j) / float(ny), 0.2);

                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
    }

#elif testNumber == 4 /*4: test the class "vec3"*/

    #include <iostream>
    #include "vec3.h"
    #include "log.h"

    using namespace std;

    int main()
    {
        vec3 lookfrom(-2,2,1);
        vec3 lookat(0,0,-1);
        vec3 vup(0,1,0);
        float vfov = 90;
        float aspect = 2;

        vec3 u, v, w;
        vec3 u1, v1, w1;
        vec3 u2, v2, w2;
        vec3 lower_left_corner;
        vec3 horizontal;
        vec3 vertical;
        vec3 origin;

        float theta = vfov*M_PI/180;
        float half_height = tan(theta/2);
        float half_width = aspect * half_height;
        origin = lookfrom;
        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);
        lower_left_corner = origin - half_width*u - half_height*v -w;
        horizontal = 2*half_width*u;
        vertical = 2*half_height*v;

        std::cout << "lookfrom: " << lookfrom << endl;
        std::cout << "lookat: " << lookat << endl;
        std::cout << "vup: " << vup << endl;
        std::cout << "vfov: " << vfov << endl;
        std::cout << "aspect: " << aspect << endl;
        std::cout << "-------------------------------------"<< endl;
        std::cout << "theta = vfov*M_PI/180: " << theta << endl;
        std::cout << "half_height = tan(theta/2): " << half_height << endl;
        std::cout << "half_width = aspect * half_height: " << half_width << endl;
        std::cout << "origin = lookfrom: " << origin << endl;
        std::cout << "w = unit_vector(lookfrom - lookat): " << w << endl;
        std::cout << "u = unit_vector(cross(vup, w)): " << u << endl;
        std::cout << "v = cross(w, u): " << v << endl;
        std::cout << "dot(w, u): " << dot(w, u) << endl;
        std::cout << "dot(u, v): " << dot(u, v) << endl;
        std::cout << "dot(v, w): " << dot(v, w) << endl;
        std::cout << "dot(vup, v): " << dot(vup, v) << endl;
        std::cout << "lower_left_corner = origin - half_width*u - half_height*v -w: " << lower_left_corner << endl;
        std::cout << "horizontal = 2*half_width*u: " << horizontal << endl;
        std::cout << "vertical = 2*half_height*v: " << vertical << endl;
/*
#if MAIN_LOG2 ==1
        std::cout << endl;
        std::cout << "v2: " << v2.e[0] << " " << v2.e[1] << " " << v2.e[2] << endl;
        std::cout << "v2-length: " << v2.length() << endl;
        std::cout << "v2-squared_length: " << v2.squared_length() << endl;

        std::cout << endl;
        std::cout << "v1*v2: " << (v1*v2).e[0] << " " << (v1*v2).e[1] << " " << (v1*v2).e[2] << endl;
        std::cout << "v1*v2-length: " << (v1*v2).length() << endl;
        std::cout << "v1*v2-squared_length: " << (v1*v2).squared_length() << endl;


        std::cout << endl;
        std::cout << "v2*v1: " << (v2*v1).e[0] << " " << (v2*v1).e[1] << " " << (v2*v1).e[2] << endl;
        std::cout << "v2*v1-length: " << (v2*v1).length() << endl;
        std::cout << "v2*v1-squared_length: " << (v2*v1).squared_length() << endl;


        std::cout << endl;
        std::cout << "v2/v1: " << (v2/v1).e[0] << " " << (v2/v1).e[1] << " " << (v2/v1).e[2] << endl;
        std::cout << "v2/v1-length: " << (v2/v1).length() << endl;
        std::cout << "v2/v1-squared_length: " << (v2/v1).squared_length() << endl;

        std::cout << endl;
        std::cout << "---------------------------------------" << endl;

        v1.operator+=(v2);
        std::cout << endl;
        std::cout << "v1=v1+v2: " << v1.e[0] << " " << v1.e[1] << " " << v1.e[2] << endl;
        std::cout << "v1=v1+v2-length: " << v1.length() << endl;
        std::cout << "v1=v1+v2-squared_length: " << v1.squared_length() << endl;

        v1.operator-=(v2);
        std::cout << endl;
        std::cout << "v1=v1-v2: " << v1.e[0] << " " << v1.e[1] << " " << v1.e[2] << endl;
        std::cout << "v1=v1-v2-length: " << v1.length() << endl;
        std::cout << "v1=v1-v2-squared_length: " << v1.squared_length() << endl;

        v1.operator*=(v2);
        std::cout << endl;
        std::cout << "v1=v1*v2: " << v1.e[0] << " " << v1.e[1] << " " << v1.e[2] << endl;
        std::cout << "v1=v1*v2-length: " << v1.length() << endl;
        std::cout << "v1=v1*v2-squared_length: " << v1.squared_length() << endl;

        v1.operator/=(v2);
        std::cout << endl;
        std::cout << "v1=v1/v2: " << v1.e[0] << " " << v1.e[1] << " " << v1.e[2] << endl;
        std::cout << "v1=v1/v2-length: " << v1.length() << endl;
        std::cout << "v1=v1/v2-squared_length: " << v1.squared_length() << endl;

        v1.operator*=(4);
        std::cout << endl;
        std::cout << "v1=v1*4: " << v1.e[0] << " " << v1.e[1] << " " << v1.e[2] << endl;
        std::cout << "v1=v1*4-length: " << v1.length() << endl;
        std::cout << "v1=v1*4-squared_length: " << v1.squared_length() << endl;

        v1.operator/=(2);
        std::cout << endl;
        std::cout << "v1=v1/2: " << v1.e[0] << " " << v1.e[1] << " " << v1.e[2] << endl;
        std::cout << "v1=v1/2-length: " << v1.length() << endl;
        std::cout << "v1=v1/2-squared_length: " << v1.squared_length() << endl;

        std::cout << endl;
        std::cout << "---------------------------------------" << endl;

        std::cout << endl;
        std::cout << "v1: " << v1.e[0] << " " << v1.e[1] << " " << v1.e[2] << endl;
        std::cout << "v1-length: " << v1.length() << endl;
        std::cout << "v1-squared_length: " << v1.squared_length() << endl;

        std::cout << endl;
        std::cout << "v2: " << v2.e[0] << " " << v2.e[1] << " " << v2.e[2] << endl;
        std::cout << "v2-length: " << v2.length() << endl;
        std::cout << "v2-squared_length: " << v2.squared_length() << endl;

        vec3 vt;

        float f = vt.dot(v1, v2);
        std::cout << endl;
        std::cout << "f=vt.dot(v1, v2): " << f << endl;

        vt = vt.cross(v1, v2);
        std::cout << endl;
        std::cout << "vt=vt.cross(v1, v2): " << vt.e[0] << " " << vt.e[1] << " " << vt.e[2] << endl;
        std::cout << "vt=vt.cross(v1, v2)-length: " << vt.length() << endl;
        std::cout << "vt=vt.cross(v1, v2)-squared_length: " << vt.squared_length() << endl;

        std::cout << endl;
        std::cout << "---------------------------------------" << endl;

        v1.make_unit_vector();
        v2.make_unit_vector();
        std::cout << endl;
        std::cout << "v1.make_unit_vector: " << v1.e[0] << " " << v1.e[1] << " " << v1.e[2] << endl;
        std::cout << "v1.make_unit_vector-length: " << v1.length() << endl;
        std::cout << "v1.make_unit_vector-squared_length: " << v1.squared_length() << endl;

        std::cout << endl;
        std::cout << "v2.make_unit_vector: " << v2.e[0] << " " << v2.e[1] << " " << v2.e[2] << endl;
        std::cout << "v2.make_unit_vector-length: " << v2.length() << endl;
        std::cout << "v2.make_unit_vector-squared_length: " << v2.squared_length() << endl;
#endif // MAIN_LOG2
*/
        return 0;
    }

#elif testNumber == 5 /*5: test the function "inline vec3& vec3::operator+=(const vec3 &v)"*/

    #include <iostream>
    #include "vec3.h"

    using namespace std;

    int main()
    {

        vec3 v1(1.0, 2.0, 3.0);
        vec3 v2(4.0, 5.0, 6.0);
        std::cout << "v1: " << v1.e[0] << " " << v1.e[1] << " " << v1.e[2] << endl;
        std::cout << "v2: " << v2.e[0] << " " << v2.e[1] << " " << v2.e[2] << endl;

        (v1.operator+=(v2)).e[0] += 100;

        std::cout << "v1==============================: " << v1.e[0] << " " << v1.e[1] << " " << v1.e[2] << endl;
        std::cout << "(v1.operator+=(v2)).e[0] += 1000: " << (v1.operator+=(v2)).e[0] << " " << (v1.operator+=(v2)).e[1] << " " << (v1.operator+=(v2)).e[2] << endl;

        return 0;
    }


#elif testNumber == 6 /*6: test rays*/

    #include <iostream>
    #include <fstream>
    #include "ray.h"

    using namespace std;

    vec3 color(const ray& r)
    {
        vec3 unit_direction = unit_vector(r.direction());
        float t = 0.5*(unit_direction.y() + 1.0);
        return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
//        return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(1.0, 0, 0.7);//white, pink
    }

    int main()
    {
        int nx = 200;
        int ny = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/RaysBackgroundY.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        vec3 lower_left_corner(-2.0, -1.0, -1.0);
        vec3 horizontal(4.0, 0.0, 0.0);
        vec3 vertical(0.0, 2.0, 0.0);
        vec3 origin(0.0, 0.0, 0.0);

        for (int j = ny-1; j >= 0; j--)
        {
            for (int i = 0; i < nx; i++)
            {
                float u = float(i) / float(nx);
                float v = float(j) / float(ny);
                ray r(origin, lower_left_corner + u*horizontal + v*vertical);
                vec3 col = color(r);
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }

    }


#elif testNumber == 7 /*7: test rays and add a sphere*/

    #include <iostream>
    #include <fstream>
    #include "ray.h"

    using namespace std;

    bool hit_sphere(const vec3& center, float radius, const ray& r)
    {
        vec3 oc = r.orgin() - center;
        float a = oc.dot(r.direction(), r.direction());
        float b = 2.0 * oc.dot(oc, r.direction());
        float c = oc.dot(oc, oc) - radius*radius;
        float discriminant = b*b - 4*a*c;
        return (discriminant > 0);
    }

    vec3 color(const ray& r)
    {
        if(hit_sphere(vec3(0,0,-1), 0.5, r))
            return vec3(1, 0, 0);

        vec3 unit_direction = unit_vector(r.direction());
        float t = 0.5*(unit_direction.y() + 1.0);
        return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
//        return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(1.0, 0, 0.7);//white, pink
    }

    int main()
    {
        int nx = 200;
        int ny = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/AddSphere.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        vec3 lower_left_corner(-2.0, -1.0, -1.0);
        vec3 horizontal(4.0, 0.0, 0.0);
        vec3 vertical(0.0, 2.0, 0.0);
        vec3 origin(0.0, 0.0, 0.0);

        for (int j = ny-1; j >= 0; j--)
        {
            for (int i = 0; i < nx; i++)
            {
                float u = float(i) / float(nx);
                float v = float(j) / float(ny);
                ray r(origin, lower_left_corner + u*horizontal + v*vertical);
                vec3 col = color(r);
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }

    }

#elif testNumber == 8 /*8: visualize the normals of sphere.*/

    #include <iostream>
    #include <fstream>
    #include "ray.h"

    using namespace std;

    float hit_sphere(const vec3& center, float radius, const ray& r){
        vec3 oc = r.orgin() - center;
        float a = oc.dot(r.direction(), r.direction());
        float b = 2.0 * oc.dot(oc, r.direction());
        float c = oc.dot(oc, oc) - radius*radius;
        float discriminant = b*b - 4*a*c;

        if (discriminant < 0){
            return -1.0;
        }
        else{
            return (-b - sqrt(discriminant)) / (2.0*a);
        }
    }

    vec3 color(const ray& r){
        float t = hit_sphere(vec3(0,0,-1), 0.5, r);
        if (t > 0.0){
            vec3 N = unit_vector(r.point_at_parameter(t) - vec3(0,0,-1));
            return 0.5*vec3(N.x()+1, N.y()+1, N.z()+1);
        }

        vec3 unit_direction = unit_vector(r.direction());
        t = 0.5*(unit_direction.y() + 1.0);
        return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
//        return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(1.0, 0, 0.7);//white, pink
    }

    int main(){
        int nx = 200;
        int ny = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/VisualizeNormals.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        vec3 lower_left_corner(-2.0, -1.0, -1.0);
        vec3 horizontal(4.0, 0.0, 0.0);
        vec3 vertical(0.0, 2.0, 0.0);
        vec3 origin(0.0, 0.0, 0.0);

        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                float u = float(i) / float(nx);
                float v = float(j) / float(ny);
                ray r(origin, lower_left_corner + u*horizontal + v*vertical);
                vec3 col = color(r);
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }

    }

#elif testNumber == 9 /*9: several spheres*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"

    using namespace std;

    vec3 color(const ray& r, hitable *world) {

        hit_record rec;
        if (world->hit(r, 0.0, (numeric_limits<float>::max)(), rec)) {
            return 0.5*vec3(rec.normal.x()+1, rec.normal.y()+1, rec.normal.z()+1);
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/SeveralSpheres.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        vec3 lower_left_corner(-2.0, -1.0, -1.0);
        vec3 horizontal(4.0, 0.0, 0.0);
        vec3 vertical(0.0, 2.0, 0.0);
        vec3 origin(0.0, 0.0, 0.0);
        hitable *list[2];
        list[0] = new sphere(vec3(0,0,-1), 0.5);
        list[1] = new sphere(vec3(0,-100.5,-1), 100);
        hitable *world = new hitable_list(list,2);

        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                float u = float(i) / float(nx);
                float v = float(j) / float(ny);
                ray r(origin, lower_left_corner + u*horizontal + v*vertical);

//                vec3 p = r.point_at_parameter(2.0);
                vec3 col = color(r, world);
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }

    }


#elif testNumber == 10 /*10: abstract class*/

    #include "test_extend1.h"
    #include "test_extend2.h"

    int main(){
        test_base *te1 = new test_extend1();
        test_base *te2 = new test_extend2();

        te1->vfb();
        te2->vfb();

        return 0;
    }
#elif testNumber == 11 /*11: antialiasing*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"

    using namespace std;

    vec3 color(const ray& r, hitable *world) {

        hit_record rec;
        if (world->hit(r, 0.0, (numeric_limits<float>::max)(), rec)) {
            return 0.5*vec3(rec.normal.x()+1, rec.normal.y()+1, rec.normal.z()+1);
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/Antialiasing.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        hitable *list[2];
        list[0] = new sphere(vec3(0,0,-1), 0.5);
        list[1] = new sphere(vec3(0,-100.5,-1), 100);
        hitable *world = new hitable_list(list,2);

        camera cam;

        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
                    ray r = cam.get_ray(u, v);

    //                vec3 p = r.point_at_parameter(2.0);
                    col += color(r, world);
                }
                col /= float(ns);
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }

    }

#elif testNumber == 12 /*12: diffuse materials*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>
//    #include <iomanip>

    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"

    using namespace std;
/*
    extern int counter1 = 0;
    extern int counter2 = 0;
    extern int counter3 = 0;
    extern int counter4 = 0;
*/
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

    vec3 color(const ray& r, hitable *world) {

//        counter1++;//5883059/59002/32920
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
/*            counter2++;//3883059/39002/12920
            if ((r.origin().x() == (float)(0)) &&
                (r.origin().y() == (float)(0)) &&
                (r.origin().z() == (float)(0)))
                counter3++;//1019651/10194/2615

            ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/log/ns100_ball2_t0+[log].txt", ios_base::app);
            outfile << "counter2: " << counter2
                    << "===ray:origin: ( " << setprecision(18) <<r.origin().x()<< ", " << r.origin().y() << ", " << r.origin().z() << "), " <<setprecision(6)
                    << "===ray:direction: ( " << r.direction().x() << ", " << r.direction().y() << ", " << r.direction().z() << "), "
                    << "===center,radius: ( " << rec.c.x() << ", " << rec.c.y() << ", " << rec.c.z() << " ), " << rec.r
                    << "===hitpoint(t): " << rec.t << endl;
*/
            vec3 target = rec.p + rec.normal + random_in_unit_sphere();
            return 0.5*color( ray(rec.p, target-rec.p), world);
        }
        else {
//            counter4++;//2000000/20000/20000
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;
        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/ns100_ball2_t0+[image].txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        hitable *list[2];
        list[0] = new sphere(vec3(0,0,-1), 0.5);
        list[1] = new sphere(vec3(0,-100.5,-1), 100);
        hitable *world = new hitable_list(list,2);

        camera cam;

        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
                    ray r = cam.get_ray(u, v);

    //                vec3 p = r.point_at_parameter(2.0);
                    col += color(r, world);
                }
                col /= float(ns);
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
/*
        std::cout << "the number of all the rays(original rays and reflected rays):    counter1: " << counter1 << endl;
        std::cout << "the number of the rays that hit the sphere:                      counter2: " << counter2 << endl;
        std::cout << "the number of the rays that hit the sphere and come from origin: counter3: " << counter3 << endl;
        std::cout << "the number of the rays that are set color:                       counter4: " << counter4 << endl;
*/
    }
#elif testNumber == 13 /*13: metal*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"

    using namespace std;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
//        std::cout << "-------------color()--1----------------" << endl;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
//            std::cout << "-------------color()--2----------------" << endl;
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()--3----------------" << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
//                std::cout << "-------------color()--4----------------" << endl;
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;
        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/metal.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        hitable *list[4];
        list[0] = new sphere(vec3(0,0,-1), 0.5, new lambertian(vec3(0.8, 0.3, 0.3)));
        list[1] = new sphere(vec3(0,-100.5,-1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[2] = new sphere(vec3(1,0,-1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 0.3));
        list[3] = new sphere(vec3(-1,0,-1), 0.5, new metal(vec3(0.8, 0.8, 0.8), 0.3));
        hitable *world = new hitable_list(list,4);
        camera cam;
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
    }
#elif testNumber == 14 /*14: dielectric*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/dielectric_reflection.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        hitable *list[5];
        list[0] = new sphere(vec3(0,0,-1), 0.5, new lambertian(vec3(0.1, 0.2, 0.5)));
        list[1] = new sphere(vec3(0,-100.5,-1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[2] = new sphere(vec3(1,0,-1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 0.0));
        list[3] = new sphere(vec3(-1,0,-1), 0.5, new dielectric(1.5));
        list[4] = new sphere(vec3(-1,0,-1), -0.45, new dielectric(1.5));
        hitable *world = new hitable_list(list,5);
        camera cam;
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }
#elif testNumber == 15 /*15: positionable camera and defocus blur*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/positionable_camera.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";
/*
        float R = cos(M_PI/4);
        hitable *list[2];
        list[0] = new sphere(vec3(-R,0,-1), R, new lambertian(vec3(0, 0, 1)));
        list[1] = new sphere(vec3( R,0,-1), R, new lambertian(vec3(1, 0, 0)));
//        list[2] = new sphere(vec3(1,0,-1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 0.0));
//        list[3] = new sphere(vec3(-1,0,-1), 0.5, new dielectric(1.5));
//        list[4] = new sphere(vec3(-1,0,-1), -0.3, new dielectric(1.5));
        hitable *world = new hitable_list(list,2);
*/
        hitable *list[5];
        list[0] = new sphere(vec3(0,0,-1), 0.5, new lambertian(vec3(0.1, 0.2, 0.5)));
        list[1] = new sphere(vec3(0,-100.5,-1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[2] = new sphere(vec3(1,0,-1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 0.0));
        list[3] = new sphere(vec3(-1,0,-1), 0.5, new dielectric(1.5));
        list[4] = new sphere(vec3(-1,0,-1), -0.45, new dielectric(1.5));
        hitable *world = new hitable_list(list,5);

        vec3 lookfrom(-2,2,1);
        vec3 lookat(0,0,-1);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 2.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }
#elif testNumber == 16 /*16: random scene*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    int sphere_counter = 0;

    hitable *random_scene_my() {
        int n = 500;
        hitable **list = new hitable *[n+1];
        list[0] = new sphere(vec3(0,-1000,0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));
        int i = 1;
        for (int a = -11; a < 11; a++) {
            for (int b = -11; b < 11; b++) {
                float choose_mat = (rand()%(100)/(float)(100));
                vec3 center(a+0.9*(rand()%(100)/(float)(100)),0.2,b+0.9*(rand()%(100)/(float)(100)));
                if ((center-vec3(4,0.2,0)).length() > 0.9) {
                    if (choose_mat < 0.8) {     //diffuse
                        list[i++] = new sphere(center, 0.2,
                                               new lambertian(vec3((rand()%(100)/(float)(100))*(rand()%(100)/(float)(100)),
                                                                   (rand()%(100)/(float)(100))*(rand()%(100)/(float)(100)),
                                                                   (rand()%(100)/(float)(100))*(rand()%(100)/(float)(100)))));
                    }
                    else if (choose_mat < 0.95) {
                        list[i++] = new sphere(center, 0.2,
                                               new metal(vec3(0.5*(1+(rand()%(100)/(float)(100))),
                                                              0.5*(1+(rand()%(100)/(float)(100))),
                                                              0.5*(1+(rand()%(100)/(float)(100)))),
                                                         0.5*(1+(rand()%(100)/(float)(100)))));
                    }
                    else {
                        list[i++] = new sphere(center, 0.2,
                                               new dielectric(1.5));
                    }
                }
            }
        }

        list[i++] = new sphere(vec3(-6, 2, -6), 2.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));
        list[i++] = new sphere(vec3(6, 2, -6), 2.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

        list[i++] = new sphere(vec3(0, 2, -7), 2.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

        list[i++] = new sphere(vec3(-2, 1, -4), 1.0, new dielectric(1.5));
        list[i++] = new sphere(vec3(2, 1, -4), 1.0, new dielectric(1.5));
        list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));

        list[i++] = new sphere(vec3(-4, 1, -2), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
        list[i++] = new sphere(vec3(4, 1, -2), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));

        list[i++] = new sphere(vec3(-6, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));
        list[i++] = new sphere(vec3(6, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));
/* */
        return new hitable_list(list, i);
    }

    hitable *random_scene() {
        int n = 500;
        hitable **list = new hitable *[n+1];
        list[0] = new sphere(vec3(0,-1000,0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));
        int i = 1;
        for (int a = -11; a < 11; a++) {
            for (int b = -11; b < 11; b++) {
                float choose_mat = (rand()%(100)/(float)(100));
                vec3 center(a+0.9*(rand()%(100)/(float)(100)),0.2,b+0.9*(rand()%(100)/(float)(100)));
                if ((center-vec3(4,0.2,0)).length() > 0.9) {
                    sphere_counter++;
                    if (choose_mat < 0.8) {     //diffuse
                        list[i++] = new sphere(center, 0.2,
                                               new lambertian(vec3((rand()%(100)/(float)(100))*(rand()%(100)/(float)(100)),
                                                                   (rand()%(100)/(float)(100))*(rand()%(100)/(float)(100)),
                                                                   (rand()%(100)/(float)(100))*(rand()%(100)/(float)(100)))));
                    }
                    else if (choose_mat < 0.95) {
                        list[i++] = new sphere(center, 0.2,
                                               new metal(vec3(0.5*(1+(rand()%(100)/(float)(100))),
                                                              0.5*(1+(rand()%(100)/(float)(100))),
                                                              0.5*(1+(rand()%(100)/(float)(100)))),
                                                         0.5*(1+(rand()%(100)/(float)(100)))));
                    }
                    else {
                        list[i++] = new sphere(center, 0.2,
                                               new dielectric(1.5));
                    }
                }
            }
        }

        list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
        list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
        list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

        return new hitable_list(list, i);
    }


    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/random_scene.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";
/*
        hitable *list[5];
        list[0] = new sphere(vec3(0,0,-1), 0.5, new lambertian(vec3(0.1, 0.2, 0.5)));
        list[1] = new sphere(vec3(0,-100.5,-1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[2] = new sphere(vec3(1,0,-1), 0.5, new metal(vec3(0.8, 0.6, 0.2), 0.0));
        list[3] = new sphere(vec3(-1,0,-1), 0.5, new dielectric(1.5));
        list[4] = new sphere(vec3(-1,0,-1), -0.45, new dielectric(1.5));
        hitable *world = new hitable_list(list,5);
*/
        hitable *world = random_scene();

        vec3 lookfrom(13,2,3);
        vec3 lookat(0,0,0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float u = float(i + (rand()%(100)/(float)(100))) / float(nx);
                    float v = float(j + (rand()%(100)/(float)(100))) / float(ny);
                    ray r = cam.get_ray(u, v);
                    col += color(r, world, 0);
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;
     std::cout << "sphere_counter: " << sphere_counter << endl;

    }
#elif testNumber == 17 /*17: ray/box intersection*/
    #include <iostream>
    #include <fstream>
    #include "ray.h"

    using namespace std;

    bool ray_box_intersection(const vec3& origin, const vec3& direction, const vec3& bl, const vec3& bh) {
        /*use for loop and array to write the procedure*/
        float t_near = -10000;
        float t_far = 10000;
        int near_flag, far_flag;

        if(direction.x() == 0) {
            if((origin.x() < bl.x()) || (origin.x() > bh.x())) {
                std::cout << "the ray is parallel to the planes and the origin X0 is not between the slabs. return false" <<endl;
                return false;
            }
        }
        if(direction.y() == 0) {
            if((origin.y() < bl.y()) || (origin.y() > bh.y())) {
                std::cout << "the ray is parallel to the planes and the origin Y0 is not between the slabs. return false" <<endl;
                return false;
            }
        }
        if(direction.z() == 0) {
            if((origin.z() < bl.z()) || (origin.z() > bh.z())) {
                std::cout << "the ray is parallel to the planes and the origin Z0 is not between the slabs. return false" <<endl;
                return false;
            }
        }

        float array1[6];
        array1[0] = (bl.x()-origin.x())/direction.x();
        array1[1] = (bh.x()-origin.x())/direction.x();
        array1[2] = (bl.y()-origin.y())/direction.y();
        array1[3] = (bh.y()-origin.y())/direction.y();
        array1[4] = (bl.z()-origin.z())/direction.z();
        array1[5] = (bh.z()-origin.z())/direction.z();

        for (int i=0; i<6; i++){
            if(array1[i] > array1[i+1]) {
                float t = array1[i];
                array1[i] = array1[i+1];
                array1[i+1] = t;
            }
            std::cout << "array1[" << i << "]:" << array1[i] <<endl;
            std::cout << "array1[" << i+1 << "]:" << array1[i+1] <<endl;
            if(array1[i] >= t_near) {t_near = array1[i]; near_flag = i;}
            if(array1[i+1] <= t_far) {t_far = array1[i+1]; far_flag = i+1;}
            if(t_near > t_far) {
                std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_near > t_far. return false" <<endl;
                return false;
            }
            if(t_far < 0) {
                std::cout << "No.(0=X;2=Y;4=Z):" << i << "  :t_far < 0. return false" <<endl;
                return false;
            }
            i++;
        }

        std::cout << "t_near: " << t_near << "   near_flag: " << near_flag <<endl;
        std::cout << "t_far: " << t_far << "   far_flag: " << far_flag <<endl;
        std::cout << "t_near,parameters: " << origin+direction*t_near << "   t_far,parameters: " << origin+direction*t_far <<endl;
        std::cout << "pass all of the tests. return ture" <<endl;
        return true;
    }

    int main() {

#if 0
        vec3 origin(0,4,2);
        vec3 direction(0.218,-0.436,0.873);
        vec3 bl(-1,2,1);
        vec3 bh(3,3,3);
#endif // 0
#if 0
        vec3 origin(0,0,0);
        vec3 direction(0.3,0.6,-1);
        vec3 bl(-2,0,-4);
        vec3 bh(2,4,-8);
#endif // 0
#if 1
        vec3 origin(6,-1,1.6);
        vec3 direction(-1,1.4,-0.2);
        vec3 bl(0,0,2);
        vec3 bh(4,2,0);

#endif // 1

        bool sec = ray_box_intersection(origin, direction, bl, bh);
        std::cout << "ray_box_intersection(0=No;1=Yes): " << sec <<endl;
    }
#elif testNumber == 18 /*18: triangle*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(1,1,1);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/triangle.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        hitable *list[5];
        list[0] = new triangle(vec3(-1,1.5,-1), vec3(-3,1.5,-1), vec3(-2,3.5,-1), new lambertian(vec3(0.3, 0.8, 0.0)));
        list[1] = new triangle(vec3(3,-0.5,-2), vec3(-1,-0.5,-4), vec3(0,5,-3), new metal(vec3(0.8, 0.6, 0.5), 0.0));
//        list[2] = new triangle(vec3(4,-0.5,4), vec3(-4,-0.5,4), vec3(0,4,4), new dielectric(4.5));
        list[2] = new triangle(vec3(1.5,-0.5,4), vec3(-1.5,-0.5,4), vec3(0,2,4), new dielectric(4.5));
        list[3] = new sphere(vec3(0,-100.5,-1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[4] = new sphere(vec3(-2,0.5,-1), 1, new lambertian(vec3(0.5, 0.7, 0.6)));
        hitable *world = new hitable_list(list,5);

        vec3 lookfrom(0,0,10);
        vec3 lookat(0,1,-1);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 30, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }
#elif testNumber == 19 /*19: polygon*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(1,1,1);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/polygon.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

//triangle1, the light pink metal one
        vec3 vertexes3_1[3];
        vertexes3_1[0] = vec3(3,-0.5,-2);
        vertexes3_1[1] = vec3(-1,-0.5,-4);
        vertexes3_1[2] = vec3(0,5,-3);
//triangle2, the green lambertian one
        vec3 vertexes3_2[3];
        vertexes3_2[0] = vec3(-1,1.5,-1);
        vertexes3_2[1] = vec3(-3,1.5,-1);
        vertexes3_2[2] = vec3(-2,3.5,-1);

//triangle3, the smaller dielectric one
        vec3 vertexes3_3[3];
        vertexes3_3[0] = vec3(1.5,-0.5,4);
        vertexes3_3[1] = vec3(-1.5,-0.5,4);
        vertexes3_3[2] = vec3(0,2,4);

/*
//triangle3, the bigger dielectric one
        vec3 vertexes3_4[3];
        vertexes3_4[0] = vec3(4,-0.5,4);
        vertexes3_4[1] = vec3(-4,-0.5,4);
        vertexes3_4[2] = vec3(0,4,4);
*/
/*//pentagon
        vec3 vertexes5[5];
        vertexes5[0] = vec3(3.3000,1.3000,0.0000);
        vertexes5[1] = vec3(2.3489,0.6090,0.0000);
        vertexes5[2] = vec3(2.7122,-0.5090,0.0000);
        vertexes2[3] = vec3(3.8878,-0.5090,0.0000);
        vertexes5[4] = vec3(4.2511,0.6090,0.0000);
*/
//pentacle
        vec3 vertexes5[5];
        vertexes5[0] = vec3(3.3000,1.3000,0.0000);
        vertexes5[1] = vec3(2.7122,-0.5090,0.0000);
        vertexes5[2] = vec3(4.2511,0.6090,0.0000);
        vertexes5[3] = vec3(2.3489,0.6090,0.0000);
        vertexes5[4] = vec3(3.8878,-0.5090,0.0000);

        hitable *list[6];
        list[0] = new polygon(vertexes3_2, 3, new lambertian(vec3(0.3, 0.8, 0.0)));
        list[1] = new sphere(vec3(0,-100.5,-1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[2] = new sphere(vec3(-2,0.5,-1), 1, new lambertian(vec3(0.5, 0.7, 0.6)));
        list[3] = new polygon(vertexes3_1, 3, new metal(vec3(0.8, 0.6, 0.5), 0.0));
        list[4] = new polygon(vertexes5, 5, new lambertian(vec3(1.0, 0.0, 0.0)));
        list[5] = new polygon(vertexes3_3, 3, new dielectric(9.0));
//        list[5] = new polygon(vertexes3_4, 3, new dielectric(9.0));
        hitable *world = new hitable_list(list,6);

        vec3 lookfrom(0,0,10);
        vec3 lookat(0,1,-1);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 30, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }
#elif testNumber == 20 /*20: box*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "box.h"
    #include "box2.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/box.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

//triangle2, the green lambertian one
        vec3 vertexes3_2[3];
        vertexes3_2[0] = vec3(1.5,0.5,1.0);
        vertexes3_2[1] = vec3(2.5,0.5,1.0);
        vertexes3_2[2] = vec3(2.0,2.0,1.0);

        hitable *list[7];
        list[0] = new sphere(vec3(0.0,-100.5,-1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
//        list[1] = new box(vec3(-2.0,-0.5,4.0), vec3(-1.0,1.0,2.0), new lambertian(vec3(0.0, 1.0, 0.5)));
        list[1] = new box2(vec3(1, 0.5, -0.5), 0, 1, 2, 1.5, vec3(-2.0,-0.5,4.0), new lambertian(vec3(0.0, 1.0, 0.5)));
        list[2] = new box(vec3(-0.25,-0.5,0.0), vec3(0.75,0.5,-1.0), new metal(vec3(0.8, 0.2, 0.2), 0.0));
        list[3] = new box(vec3(-5.0,-0.5,-5.0), vec3(5.0,3.0,-6.0), new metal(vec3(0.8, 0.6, 0.4), 0.0));
        list[4] = new sphere(vec3(2.0,0.0,1.0), 0.5, new lambertian(vec3(0.5, 0.7, 0.6)));
        list[5] = new sphere(vec3(0.75,-0.25,5.0), 0.25, new lambertian(vec3(0.8, 0.7, 0.6)));
        list[6] = new polygon(vertexes3_2, 3, new lambertian(vec3(0.3, 0.8, 0.0)));
        hitable *world = new hitable_list(list,7);

        vec3 lookfrom(0,0,12);
        vec3 lookat(0,1,-1);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }

#elif testNumber == 21 /*21: linear equation*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "box.h"

    using namespace std;

    int main(){
#if 0
        vec3 u(3, 2, 1);
        vec3 v(2, 2, 2);
        vec3 w(4, 3, 3);
        vec3 r_before(500, 190, 90);
#endif //
#if 1
        vec3 u(0, 1, 1);
        vec3 v(-1, 1, 0);
        vec3 w(1, 2, 1);
        vec3 r_before(2, 5, 3);
#endif //
#if 0
        vec3 u(1, 0, 0);
        vec3 v(0, 1, 0);
        vec3 w(0, 0, 1);
        vec3 r_before(2, 5, 3);
#endif //

        vec3 r_after = vector_trans(r_before, u, v, w);
        vec3 r_back = vector_trans_back(r_after, u, v, w);

        std::cout << "r_before: " << r_before << endl;
        std::cout << "r_after: " << r_after << endl;
        std::cout << "r_back: " << r_back << endl;
    }
#elif testNumber == 22 /*22: get_vector_v*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>
    #include "float.h"

    #include "vec3.h"

    using namespace std;

    int main(){
#if 0
        vec3 vector_u(1, 0, 0);
        float angle = 0;
#endif
#if 0
        vec3 vector_u(1, 0, 0);
        float angle = 45;
#endif
#if 0
        vec3 vector_u(1, 0, 1);
        float angle = 30;
#endif
#if 1
        vec3 vector_u(1, 1, 0);
        float angle = 30;
#endif

        vec3 vector_v = get_vector_v(vector_u, angle);

        std::cout << "vector_u: " << vector_u << endl;
        std::cout << "vector_v: " << vector_v << endl;
        std::cout << "cos(angle): " << cos(angle*M_PI/180) << endl;
    }
#elif testNumber == 23 /*23: box2()*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>
    #include "float.h"

    #include "vec3.h"

    using namespace std;

    int main(){
#if 1
        vec3 u(1, 1, 1);
        float an = 0;
        float a = 3;
        float b = 4;
        float c = 5;
        vec3 p(0, 0, 2);
#endif // 1
        vec3 vector_u;
        vec3 vector_v;
        vec3 vector_w;
        vec3 vertex_l;
        vec3 vertex_h;
        vec3 normals[6];

        vector_u = unit_vector(u);
        vector_v = unit_vector(get_vector_v(vector_u, an));
        vector_w = unit_vector(cross(vector_u, vector_v));

        vertex_l = p;
        vertex_h = p + a*vector_u + b*vector_v - c*vector_w;

        normals[0] = vector_trans(vec3(-1, 0, 0), vector_u, vector_v, vector_w);//left
        normals[1] = vector_trans(vec3(1, 0, 0), vector_u, vector_v, vector_w);//right
        normals[2] = vector_trans(vec3(0, 1, 0), vector_u, vector_v, vector_w);;//up
        normals[3] = vector_trans(vec3(0, -1, 0), vector_u, vector_v, vector_w);;//down
        normals[4] = vector_trans(vec3(0, 0, 1), vector_u, vector_v, vector_w);;//front
        normals[5] = vector_trans(vec3(0, 0, -1), vector_u, vector_v, vector_w);;//back

        std::cout << "vector_u: " << vector_u << endl;
        std::cout << "vector_v: " << vector_v << endl;
        std::cout << "vector_w: " << vector_w << endl;
        std::cout << "vertex_l: " << vertex_l << endl;
        std::cout << "vertex_h: " << vertex_h << endl;
        for (int i=0; i<6; i++) {
            std::cout << "normals[" << i << "]:" << normals[i] << endl;
        }
    }
#elif testNumber == 24 /*24: any box*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "box.h"
    #include "box2.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/box_any.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";
/*
//triangle2, the green lambertian one
        vec3 vertexes3_2[3];
        vertexes3_2[0] = vec3(1.5,0.5,1.0);
        vertexes3_2[1] = vec3(2.5,0.5,1.0);
        vertexes3_2[2] = vec3(2.0,2.0,1.0);
*/
        hitable *list[3];
        list[0] = new sphere(vec3(0.0,-100.5,-1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[1] = new box(vec3(-1.5,-0.5,4.0), vec3(-1.0,1.5,3.5), new lambertian(vec3(0.0, 1.0, 0.5)));
/*
        list[2] = new box2(vec3(1, 0.5, -0.5), 0, 0.5, 0.5, 2, vec3(0.0,-0.5,4.0), new lambertian(vec3(0.0, 1.0, 0.5)));
        list[3] = new box2(vec3(1, 0.5, -0.5), 30, 0.5, 0.5, 2, vec3(0.0,-0.5,4.0), new lambertian(vec3(0.8, 0.1, 0.5)));
        list[4] = new box2(vec3(1, 0.5, -0.5), 45, 0.5, 0.5, 2, vec3(0.0,-0.5,4.0), new lambertian(vec3(0.6, 0.3, 0.5)));
        list[5] = new box2(vec3(1, 0.5, -0.5), 60, 0.5, 0.5, 2, vec3(0.0,-0.5,4.0), new lambertian(vec3(0.4, 0.5, 0.5)));
        list[6] = new box2(vec3(1, 0.5, -0.5), 90, 0.5, 0.5, 2, vec3(0.0,-0.5,4.0), new lambertian(vec3(0.2, 0.7, 0.5)));
        list[7] = new box2(vec3(1, 0, -0.5), 90, 8.0, 0.5, 3.5, vec3(-4.0,-0.5,0.0), new metal(vec3(0.8, 0.8, 0.8), 0.0));
        list[8] = new box2(vec3(1, 0, 0.5), -90, 8.0, 0.5, 3.5, vec3(-1.0,-0.5,-4.0), new metal(vec3(0.8, 0.8, 0.8), 0.0));
*/

//        list[2] = new box2(vec3(1, 0, -0.5), 0, 1.0, 0.5, 2, vec3(0.0,-0.5,4.0), new lambertian(vec3(0.0, 0.1, 0.5)));
        list[2] = new box2(vec3(1, 0.5, -0.5), 60, 1.0, 0.5, 2, vec3(0.0,-0.5,4.0), new lambertian(vec3(0.0, 0.1, 0.5)));

/*
        list[2] = new box(vec3(-0.25,-0.5,0.0), vec3(0.75,0.5,-1.0), new metal(vec3(0.8, 0.2, 0.2), 0.0));
        list[3] = new box(vec3(-5.0,-0.5,-5.0), vec3(5.0,3.0,-6.0), new metal(vec3(0.8, 0.6, 0.4), 0.0));
        list[4] = new sphere(vec3(2.0,0.0,1.0), 0.5, new lambertian(vec3(0.5, 0.7, 0.6)));
        list[5] = new sphere(vec3(0.75,-0.25,5.0), 0.25, new lambertian(vec3(0.8, 0.7, 0.6)));
        list[6] = new polygon(vertexes3_2, 3, new lambertian(vec3(0.3, 0.8, 0.0)));
*/
        hitable *world = new hitable_list(list,3);

        vec3 lookfrom(0,0,12);
        vec3 lookat(0,1,-1);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }

#elif testNumber == 25 /*25: ellipsoid*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "box.h"
    #include "box2.h"
    #include "quadratic_ellipsoid.h"
    #include "quadratic.h"
    #include "quadratic_paraboloid.h"
    #include "quadratic_hyperbolic_paraboloid.h"
    #include "elliptic_plane.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/ellipsoid.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";
/*
//triangle2, the green lambertian one
        vec3 vertexes3_2[3];
        vertexes3_2[0] = vec3(1.5,0.5,1.0);
        vertexes3_2[1] = vec3(2.5,0.5,1.0);
        vertexes3_2[2] = vec3(2.0,2.0,1.0);
*/
        hitable *list[4];
        list[0] = new sphere(vec3(0.0,-100.5,-1), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[1] = new quadratic(vec3(-4, 2.5, -1), 1, 1.5, 2, 1, 1, 10, new lambertian(vec3(1.0, 0.0, 0.0)));
        list[2] = new quadratic(vec3(0, 2.5, -1), 1, 1.5, 2, 1, 2, 10, new lambertian(vec3(0.0, 1.0, 0.0)));
        list[3] = new quadratic(vec3(4, 2.5, -1), 1, 1.5, 2, 1, 3, 10, new lambertian(vec3(0.0, 0.0, 1.0)));
//        list[1] = new box2(vec3(0.3, 0, -1), 0, 12.0, 0.5, 3.5, vec3(-4.0,-0.5,2.0), new lambertian(vec3(0.0, 1.0, 0.5)));
//        list[2] = new box2(vec3(-0.5, 0, -1), 0, 12.0, 0.5, 3.5, vec3(4.0,-0.5,2.0), new metal(vec3(0.8, 0.8, 0.8), 0.0));
/*
        list[1] = new quadratic_ellipsoid(vec3(-6.4, 1, 0), 0.7, 1.0, 1.5, new lambertian(vec3(0.0, 0.1, 0.5)));
        list[2] = new quadratic(vec3(-4.5, 1.5, 0), 0.5, 0.7, 1, -1, 1, 1.5, new lambertian(vec3(0.1, 0.1, 0.5)));
        list[3] = new quadratic(vec3(-2.6, 1.5, 0), 0.15, 0.2, 0.2, -1, -1, 1.5, new lambertian(vec3(0.2, 0.1, 0.5)));
        list[4] = new quadratic(vec3(-1.0, 1.5, 0), 0.5, 0.7, 1, -1, 0, 1.5, new lambertian(vec3(0.3, 0.1, 0.5)));
//        list[3] = new elliptic_plane(vec3(0.5, 1.5, 0), vec3(0, 1, 1), 0.5, 0.5, 0.5, new lambertian(vec3(0.8, 0.1, 0.5)));
        list[5] = new quadratic(vec3(0.6, 1.5, 0), 0.5, 0.7, 1, 0, 1, 1.5, new lambertian(vec3(0.5, 0.1, 0.5)));
        list[6] = new quadratic_paraboloid(vec3(2.5, 1.5, 0), 1, 2, 1.5, new lambertian(vec3(0.7, 0.4, 0.5)));
        list[7] = new quadratic_hyperbolic_paraboloid(vec3(5.0, 1.5, 0), 0.25, 2, 1.5, 1, new lambertian(vec3(0.9, 0.2, 0.5)));
        list[8] = new sphere(vec3(0, 2, -6), 2.5, new metal(vec3(0.5, 1.0, 0.5), 0.0));
*/
        hitable *world = new hitable_list(list,4);

        vec3 lookfrom(0,5,20);
        vec3 lookat(0,1.5,0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }
#elif testNumber == 26 /*26: inverse mapping*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "box.h"
    #include "box2.h"
    #include "quadratic_ellipsoid.h"
    #include "quadratic.h"
    #include "quadratic_paraboloid.h"
    #include "quadratic_hyperbolic_paraboloid.h"
    #include "elliptic_plane.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/inverse_mapping.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

//triangle
        vec3 vertexes3_1[3];
        vertexes3_1[0] = vec3(-6.75, 0.0, 0.0);
        vertexes3_1[1] = vec3(-4.25, 0.0, 0.0);
        vertexes3_1[2] = vec3(-5.25, 2.5, 0.0);

//quadrilateral1
        vec3 vertexes4_1[4];
        vertexes4_1[0] = vec3(-6.5, 3.0 ,0.0);
        vertexes4_1[1] = vec3(-1.5, 3.0 ,0.0);
        vertexes4_1[2] = vec3(-1.5, 5.5 ,0.0);
        vertexes4_1[3] = vec3(-6.5, 5.5 ,0.0);
//quadrilateral2
        vec3 vertexes4_2[4];
        vertexes4_2[0] = vec3(1.75, 3.0 ,0.0);
        vertexes4_2[1] = vec3(6.5, 3.0 ,0.0);
        vertexes4_2[2] = vec3(5.5, 5.5 ,0.0);
        vertexes4_2[3] = vec3(3.0, 4.0 ,0.0);

        hitable *list[4];
        list[0] = new sphere(vec3(0.0,-100,0), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[1] = new polygon(vertexes3_1, 3, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[2] = new polygon(vertexes4_1, 4, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[3] = new box2(vec3(-1.0, 0, -0.5), 0, 20.0, 0.5, 6.5, vec3(15.0, 0.0, 0.0), new metal(vec3(0.8, 0.8, 0.8), 0.0));
/*
        list[1] = new sphere(vec3(-4.75, 1.25, 0), 1.25, new lambertian(vec3(0.8, 0.8, 0.0)));
//        list[2] = new polygon(vertexes3_1, 3, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[2] = new polygon(vertexes4_1, 4, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[3] = new polygon(vertexes4_2, 4, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[4] = new elliptic_plane(vec3(0.0, 1.25, 0), vec3(0, 0, 1), 1.25, 1.25, 1.25, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[5] = new quadratic(vec3(4.75, 1.25, 0), 1.25, 1.25, 1.25, 0, 1, 1.25, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[6] = new box2(vec3(-1.0, 0, -0.5), 0, 20.0, 0.5, 6.5, vec3(15.0, 0.0, 0.0), new metal(vec3(0.8, 0.8, 0.8), 0.0));
//        list[7] = new quadratic(vec3(5.5, 2.5, 0), 1.25, 2.5, 1.25, -1, 0, 2.5, new lambertian(vec3(0.8, 0.8, 0.0)));
*/
        hitable *world = new hitable_list(list,4);

        vec3 lookfrom(0, 5, 20);
        vec3 lookat(0, 2.5, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }
#elif testNumber == 27 /*27: find the roots of quartic equation*/

    #include <iostream>
    #include <limits>
    #include <stdlib.h>
    #include <math.h>

    #include "vec3.h"

    using namespace std;

    int main(){
//        float *roots;
//        r = roots_quadratic_equation(1.0, -3.0, 2.0);
//        roots = roots_cubic_equation(1.0, -5.0, 4.0, 0.0);
//        r = roots_cubic_equation(1.0, -6.0, 11.0, -6.0);//1,2,3
//        r = roots_cubic_equation(1.0, 0.0, -3.0, 2.0);//1,1,-2
//        r = roots_cubic_equation(1.0, 4.0, 5.0, 2.0); //-1,-1,-2
//        r = roots_cubic_equation(1.0, -3.0, 3.0, -1.0);
//        r = roots_quartic_equation(1.0, -4.0, 6.0, -4.0, 1.0);//1,1,1,1
//        roots = roots_quartic_equation(1.0, -10.0, 35.0, -50.0, 24.0);//1,2,3,4
//        roots = roots_quartic_equation(1.0, 0.0, -5.0, 0.0, 4.0);//1,-1,2,-2
//        r = roots_quartic_equation(1.0, 0.0, 0.0, 0.0, 0.0);//0,0,0,0
//        roots = roots_quartic_equation(210.243698, -1204.56726, 2743.98389, -2917.42529, 1210.08252);
//        roots = roots_quartic_equation(150.066818, -856.130798, 2094.49243, -2491.8606, 1210.08252);
//        roots = roots_quartic_equation(150.421722, -858.00708, 2097.89111, -2494.13525, 1210.08252);
//        roots = roots_quartic_equation(179.835159, -1024.26025, 2406.73975, -2701.34473, 1210.08252);
//        roots = roots_quartic_equation(0.554273665, -2.59592271, 5.93173456, -1.98254371, 0.0000005419);//0.39751,2.73335e-007
//        float roots[8] = {7, 7, 3, 2, 1, 10, 20, -1};

        double roots[5];
//        if (roots_quartic_equation2(159.627808, -909.279724, 2056.90576, -2171.1167, 891.264038, roots)) {
//        if (roots_quartic_equation2(1.4656595, -1.43020415, 4.95183134, -2.24647474, 0.0000393861628, roots)) {
//        if (roots_quartic_equation2(double(0.530991495), double(-1.24407244), double(784.004944), double(-2238.29272), double(1599), roots)) {
//        if (roots_quartic_equation2(double(0.511606872), double(-1.2098521), double(784.000183), double(-2238.30859), double(1599), roots)) {
//        if (roots_quartic_equation2(double(1.37457093e-007), double(-1.42776053e-005), double(785.264893), double(-2239.96143), double(1599), roots)) {
//        if (roots_quartic_equation2(double(7.61253727e-009), double(-2.88689307e-006), double(9.76842785), double(6.62989712), double(0.000253806094), roots)) {
//        if (roots_quartic_equation2_rain(double(2.0219367237839236e-005), double(-0.005427479092524834), double(792.74193737627502), double(-2268.5334587097168), double(1975), roots)) {
//        if (roots_quartic_equation2_rain(double(3.2079723767212276e-011), double(-2.4263054315687427e-007), double(882.03654864296459), double(-2246.4140176773071), double(1975), roots)) {
//        if (roots_quartic_equation2_rain((long double)(0.000000000032079723767212276), (long double)(-0.000000242630543156874274263054315687427), (long double)(882.03654864296459), (long double)(-2246.4140176773071), (long double)(1975), roots)) {
//        if (roots_quartic_equation2_rain((double)(0.00000000000011500021464039682), (double)(0.0000000035546470930310822), (double)(852.36045831028866), (double)(-2245.4008877277374), (double)(1975), roots)) {
        if (roots_quartic_equation2_rain(0.0000000000013403606699015957, 0.00000002242274461418392, 866.7537054248055, -2245.2320790290833, 1975, roots)) {
            for (int i=0; i<(roots[0]+1); i++) {
                std::cout << "roots[" << i << "]=" << roots[i] << endl;
            }
        }

/*
        float temp;
        for (int i=1; i<int(roots[0]); i++) {
            for (int j=i+1; j<int(roots[0])+1; j++) {
                if (roots[i] > roots[j]) {
                    temp = roots[i];
                    roots[i] = roots[j];
                    roots[j] = temp;
                }
            }
        }

        for (int i=0; i<(roots[0]+1); i++) {
            std::cout << "roots[" << i << "]=" << roots[i] << endl;
        }
*/
    }
#elif testNumber == 28 /*28: tori*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "box.h"
    #include "box2.h"
    #include "quadratic_ellipsoid.h"
    #include "quadratic.h"
    #include "quadratic_paraboloid.h"
    #include "quadratic_hyperbolic_paraboloid.h"
    #include "elliptic_plane.h"
    #include "tori.h"
    #include "quadratic_cylinder_all.h"
    #include "tori_part.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/tori.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        hitable *list[7];

//        list[0] = new tori_part(vec3(0, 3.2, 0), 3.2, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 30, 390);


        list[0] = new tori(vec3(0, 3.2, 0), 3.2, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)));
        list[1] = new tori(vec3(0, 4.8, 0), 1.2, 0.2, new lambertian(vec3(1.0, 0.0, 0.0)));
        list[2] = new quadratic_cylinder_all(vec3(0, 3.2, 0), 0.2, 0.2, 0.20001, 2.9,
                                             new lambertian(vec3(1.0, 0.0, 0.0)), vec3(0, 1, 0), 0);
        list[3] = new quadratic_cylinder_all(vec3(0, 3.2, 0), 0.2, 0.2, 0.20001, 2.9,
                                             new lambertian(vec3(1.0, 0.0, 0.0)), vec3(1, 0, 0), 0);
        list[4] = new quadratic_cylinder_all(vec3(0, 4.8, 0), 0.2, 0.2, 0.20001, 1.3,
                                             new lambertian(vec3(1.0, 0.0, 0.0)), vec3(1, 0, 0), 0);
        list[5] = new tori_part(vec3(-1.6, 3.2, 0), 1.6, 0.2, new lambertian(vec3(1.0, 0.0, 0.0)), 240, 360);
        list[6] = new tori_part(vec3(1.6, 3.2, 0), 1.6, 0.2, new lambertian(vec3(1.0, 0.0, 0.0)), 180, 300);

/*
        list[0] = new tori_part(vec3(-1.6, 1.2, 0), 1.6, 0.2, new lambertian(vec3(1.0, 0.0, 0.0)), 0, 360);
        list[1] = new quadratic_cylinder_all(vec3(0, 4.8, 0), 1, 1, 1, 1.3,
                                             new lambertian(vec3(1.0, 0.0, 0.0)), vec3(1, 0, 0), 0);
        list[2] = new quadratic_cylinder_all(vec3(3, 4.8, 0), 1, 1, 1, 1.3,
                                             new lambertian(vec3(1.0, 0.0, 0.0)), vec3(1, 0, -0.2), 0);
        list[3] = new quadratic_cylinder_all(vec3(3, 1.8, 0), 1, 1, 1, 1.3,
                                             new lambertian(vec3(1.0, 0.0, 0.0)), vec3(1, 0, -0.4), 30);
*/
        hitable *world = new hitable_list(list,7);

        vec3 lookfrom(0, 3, 5);
        vec3 lookat(0, 3, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 80, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }

#elif testNumber == 29 /*29: cylinder all*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "box.h"
    #include "box2.h"
    #include "quadratic_ellipsoid.h"
    #include "quadratic.h"
    #include "quadratic_paraboloid.h"
    #include "quadratic_hyperbolic_paraboloid.h"
    #include "elliptic_plane.h"
    #include "tori.h"
    #include "quadratic_cylinder_all.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/cylinder_all.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        hitable *list[2];
//        list[0] = new sphere(vec3(0.0,-100,0), 100, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[0] = new quadratic(vec3(-3.75, 3.25, 0), 1.25, 1.25, 1.25, 0, 1, 1.25, new lambertian(vec3(0.8, 0.8, 0.0)));
        list[1] = new quadratic_cylinder_all(vec3(3.75, 3.25, 0), 1.25, 1.25, 1.25, 1.25,
                                             new lambertian(vec3(0.8, 0.8, 0.0)), vec3(1, 1, 0), 0);
//        list[3] = new quadratic_cylinder_all(vec3(3.75, 1.25, 0), 1.25, 1.25, 1.25, 1.25, new lambertian(vec3(0.8, 0.8, 0.0)), vec3(1, 1, 0), 30);
        hitable *world = new hitable_list(list,2);

        vec3 lookfrom(0, 5, 20);
        vec3 lookat(0, 3, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }

#elif testNumber == 30 /*30: get_vector_vw()*/

    #include <iostream>
    #include <limits>
    #include <stdlib.h>
    #include <math.h>

    #include "vec3.h"

    using namespace std;

    int main(){
    /*
        vec3 vector_u = vec3 (0.5, 1.0, 0.0);
        float angle = 30;
        vec3 vector_v;
        vec3 vector_w;

        get_vector_vw(vector_u, angle, vector_v, vector_w);
    */
        vec3 x = vec3(1.0, 0.0, 0.0);
        vec3 pc = vec3(-1.0, -1.0, 0.0);

        float cos_theta = dot(pc, x) / (pc.length() * x.length());
        float theta = acos(cos_theta)*180/(M_PI);

        std::cout << "cos_theta:" << cos_theta << endl;
        std::cout << "theta:" << theta << endl;
    }
#elif testNumber == 31 /*31: tori_part_all*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "box.h"
    #include "box2.h"
    #include "quadratic_ellipsoid.h"
    #include "quadratic.h"
    #include "quadratic_paraboloid.h"
    #include "quadratic_hyperbolic_paraboloid.h"
    #include "elliptic_plane.h"
    #include "tori.h"
    #include "quadratic_cylinder_all.h"
    #include "tori_part.h"
    #include "tori_part_all.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/tori_part_all.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        hitable *list[4];
//        list[0] = new tori_part_all(vec3(0.0, 2.0, 0), 2.0, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 90, 360, vec3(1, 1, 1), 270);

//        list[0] = new tori_part(vec3(-3.0, 1.2, 0), 0.8, 0.2, new lambertian(vec3(1.0, 0.0, 0.0)), 30, 360);

        list[0] = new tori_part_all(vec3(0.0, 8.25, -2.121), 3, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 30, 360, vec3(0, 1, 1), 0);
        list[1] = new tori_part_all(vec3(0.0, 5.129, 0), 1.4, 0.2, new lambertian(vec3(1.0, 0.0, 0.0)), 30, 360, vec3(1, 0, 1), 0);
        list[2] = new tori_part_all(vec3(0.0, 2.129, 0), 2.2, 0.2, new lambertian(vec3(1.0, 1.0, 0.0)), 30, 360, vec3(1, 0, 0.5), 0);
        list[3] = new quadratic_cylinder_all(vec3(0, 0.329, 0), 0.4, 0.4, 0.40001, 4,
                                             new lambertian(vec3(0.0, 1.0, 1.0)), vec3(1, 0, 0.25), 0);

//        list[2] = new tori_part_all(vec3(3.0, 1.2, 0), 0.8, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 30, 360, vec3(1, 0, 0), 0);
//        list[3] = new tori_part_all(vec3(-3.0, 4.2, 0), 0.8, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 30, 360, vec3(0, 0, 1), 0);
//        list[4] = new tori_part_all(vec3(0.0, 4.2, 0), 0.8, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 30, 360, vec3(0.5, 0, 1), 0);
//        list[5] = new tori_part_all(vec3(3.0, 4.2, 0), 0.8, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 30, 360, vec3(0.5, -0.5, 1), 0);
/*
        list[0] = new tori_part(vec3(-3.0, 1.2, 0), 0.8, 0.2, new lambertian(vec3(1.0, 0.0, 0.0)), 0, 360);
        list[1] = new tori_part_all(vec3(0.0, 1.2, 0), 0.8, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 0, 360, vec3(1, 0, 0), 0);
        list[2] = new tori_part_all(vec3(3.0, 1.2, 0), 0.8, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 0, 360, vec3(0, 1, 0), 0);
        list[3] = new tori_part_all(vec3(-3.0, 3.2, 0), 0.8, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 0, 360, vec3(0, 0, 1), 0);
        list[4] = new tori_part_all(vec3(0.0, 3.2, 0), 0.8, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 0, 360, vec3(0, 0, 1), 30);
        list[5] = new tori_part_all(vec3(3.0, 3.2, 0), 0.8, 0.2, new lambertian(vec3(0.0, 1.0, 0.0)), 0, 360, vec3(0, 0, 1), 60);
*/
        hitable *world = new hitable_list(list,4);

        vec3 lookfrom(0.0, 3, 5);
        vec3 lookat(0, 3, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 80, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }

#elif testNumber == 32 /*32: acos(z_n, u_xoz)*/

    #include <iostream>
    #include <limits>
    #include <stdlib.h>
    #include <math.h>

    #include "vec3.h"

    using namespace std;

    int main(){
//        vec3 u = vec3(1.0, 1.0, -1.0); //phi=45
//        vec3 u = vec3(1.0, 1.0, 0.0); //phi=90
//        vec3 u = vec3(1.0, 1.0, 1.0); //phi=135
//        vec3 u = vec3(-1.0, 1.0, -1.0); //phi=45
        vec3 u = vec3(-1.0, 1.0, 0.0); //phi=90
//        vec3 u = vec3(-1.0, 1.0, 1.0); //phi=135
        vec3 z_n = vec3(0.0, 0.0, -1.0);
        vec3 u_xoz = vec3(u.x(), 0.0, u.z());

        float cos_phi = dot(z_n, u_xoz) / (z_n.length()*u_xoz.length());
        float phi = acos(cos_phi)*180/(M_PI);

        std::cout << "cos_phi:" << cos_phi << endl;
        std::cout << "phi:" << phi << endl;
    }

#elif testNumber == 33 /*33: quartic_blend_cylinder*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "quartic_blend_cylinder.h"
    #include "quadratic_cylinder_all.h"
    #include "quadratic.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/quartic_blend_cylinder.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        hitable *list[1];
        list[0] = new quartic_blend_cylinder(vec3(0, 3, 0), 1, 1, 3, vec3(0, 3, 0), 1, 1, 3, 0.5, 0.5, new lambertian(vec3(0.0, 1.0, 0.0)));
/*
        list[0] = new quadratic_cylinder_all(vec3(0, 3, 0), 1, 1, 1.00001, 3,
                                             new lambertian(vec3(1.0, 0.0, 0.0)), vec3(0, 1, 0), 0);
        list[1] = new quadratic_cylinder_all(vec3(0, 3, 0), 1, 1, 1.00001, 3,
                                             new lambertian(vec3(1.0, 0.0, 1.0)), vec3(1, 0, 0), 0);
*/
//        list[2] = new quartic_blend_cylinder(vec3(0, 3, 0), 1, 1, 3, vec3(0, 3, 0), 1, 1, 3, 0.5, 0.5, new lambertian(vec3(0.0, 1.0, 0.0)));
/*
        list[0] = new quadratic(vec3(0, 3, 0), 1.414, 1.414, 1, 1, 1, 0.8, new lambertian(vec3(1.0, 0.0, 0.0)));
        list[1] = new quadratic(vec3(0, 3, 0), 1.414, 1.414, 1, 1, 2, 0.8, new lambertian(vec3(0.0, 1.0, 0.0)));
        list[2] = new quadratic(vec3(0, 3, 0), 1.414, 1.414, 1, 1, 3, 0.8, new lambertian(vec3(1.0, 1.0, 0.0)));
*/
        hitable *world = new hitable_list(list,1);

        vec3 lookfrom(0, 6, 3);
        vec3 lookat(0, 3, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 80, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }

#elif testNumber == 34 /*34: superellipsoid*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "quartic_blend_cylinder.h"
    #include "quadratic_cylinder_all.h"
    #include "quadratic_ellipsoid.h"
    #include "superellipsoid.h"
    #include "box.h"
    #include "superquadratic.h"
    #include "superhyperboloid.h"
    #include "supertoroid.h"
    #include "blobs.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 400;
        int ny = 200;
        int ns = 100;

//        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/superhyperboloid-two.txt", ios_base::out);
        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/blobs.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        hitable *list[5];
//        list[0] = new superellipsoid(vec3(0, 3, 0), 3, 4.5, 3, 4, 1, 6, 0.0001, new lambertian(vec3(0.0, 1.0, 0.0)));
//        list[0] = new superhyperboloid(vec3(0, 3, 0), 3, 4.5, 3, 3, 0.1, -1, 1, 4.5, 16, 0.001, new lambertian(vec3(0.0, 1.0, 0.0)));
//        list[0] = new superhyperboloid(vec3(0, 3, 0), 1, 1, 1, 1, 1, -1, -1, 4, 16, 0.001, new lambertian(vec3(0.0, 1.0, 0.0)));
//        list[0] = new supertoroid(vec3(0, 3, 0), 1, 1, 1, 6, 1, 2, 16, 0.0001, new lambertian(vec3(0.0, 1.0, 0.0)));
//        list[0] = new superquadratic(vec3(0, 3, 0), 3, 4.5, 3, 0.25, 30, 0.25, new lambertian(vec3(0.0, 1.0, 0.0)));
//        list[0] = new box(vec3(-3, 0, -4), vec3(3, 6, 4), new metal(vec3(0.0, 1.0, 0.0), 0.0));
//        list[0] = new quadratic_ellipsoid(vec3(0, 3, 0), 1, 2, 3, new lambertian(vec3(0.0, 1.0, 0.0)));
/*
        list[0] = new blobs(vec3(-5.2, 1.75, 0), vec3(-5.2, 4.25, 0), -4, 1.2, -4, 1.2, 1.55, 5, 0.0001, new lambertian(vec3(0.0, 1.0, 0.0)));
        list[1] = new blobs(vec3(-2.6, 1.75, 0), vec3(-2.6, 4.25, 0), -2, 1.2, -2, 1.2, 1.55, 5, 0.0001, new lambertian(vec3(0.0, 1.0, 0.0)));
        list[2] = new blobs(vec3(0, 1.75, 0), vec3(0, 4.25, 0), -1, 1.2, -1, 1.2, 1.55, 5, 0.0001, new lambertian(vec3(0.0, 1.0, 0.0)));
        list[3] = new blobs(vec3(2.6, 1.75, 0), vec3(2.6, 4.25, 0), -0.5, 1.2, -0.5, 1.2, 1.55, 5, 0.0001, new lambertian(vec3(0.0, 1.0, 0.0)));
        list[4] = new blobs(vec3(5.2, 1.75, 0), vec3(5.2, 4.25, 0), -0.25, 1.2, -0.25, 1.2, 1.55, 5, 0.0001, new lambertian(vec3(0.0, 1.0, 0.0)));
*/
/*
        list[0] = new blobs(vec3(-5.2, 1.75, 0), vec3(-5.2, 4.25, 0), -4, 1, -4, 1, 1.35, 5, 0.0001, new lambertian(vec3(1.0, 0.0, 0.0)));
        list[1] = new blobs(vec3(-2.6, 1.75, 0), vec3(-2.6, 4.25, 0), -2, 1, -2, 1, 1.35, 5, 0.0001, new lambertian(vec3(1.0, 0.0, 0.0)));
        list[2] = new blobs(vec3(0, 1.75, 0), vec3(0, 4.25, 0), -1, 1, -1, 1, 1.35, 5, 0.0001, new lambertian(vec3(1.0, 0.0, 0.0)));
        list[3] = new blobs(vec3(2.6, 1.75, 0), vec3(2.6, 4.25, 0), -0.5, 1, -0.5, 1, 1.35, 5, 0.0001, new lambertian(vec3(1.0, 0.0, 0.0)));
        list[4] = new blobs(vec3(5.2, 1.75, 0), vec3(5.2, 4.25, 0), -0.25, 1, -0.25, 1, 1.35, 5, 0.0001, new lambertian(vec3(1.0, 0.0, 0.0)));
*/
        list[0] = new blobs(vec3(-5.2, 1.75, 0), vec3(-5.2, 4.25, 0), -4, 1.4, -4, 1.4, 1.75, 5, 0.0001, new lambertian(vec3(0.0, 0.0, 1.0)));
        list[1] = new blobs(vec3(-2.6, 1.75, 0), vec3(-2.6, 4.25, 0), -2, 1.4, -2, 1.4, 1.75, 5, 0.0001, new lambertian(vec3(0.0, 0.0, 1.0)));
        list[2] = new blobs(vec3(0, 1.75, 0), vec3(0, 4.25, 0), -1, 1.4, -1, 1.4, 1.75, 5, 0.0001, new lambertian(vec3(0.0, 0.0, 1.0)));
        list[3] = new blobs(vec3(2.6, 1.75, 0), vec3(2.6, 4.25, 0), -0.5, 1.4, -0.5, 1.4, 1.75, 5, 0.0001, new lambertian(vec3(0.0, 0.0, 1.0)));
        list[4] = new blobs(vec3(5.2, 1.75, 0), vec3(5.2, 4.25, 0), -0.25, 1.4, -0.25, 1.4, 1.75, 5, 0.0001, new lambertian(vec3(0.0, 0.0, 1.0)));

        hitable *world = new hitable_list(list,5);

        vec3 lookfrom(0, 3, 20);
        vec3 lookat(0, 3, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }

#elif testNumber == 35 /*35: rain*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "quartic_blend_cylinder.h"
    #include "quadratic_cylinder_all.h"
    #include "quadratic_ellipsoid.h"
    #include "superellipsoid.h"
    #include "box.h"
    #include "superquadratic.h"
    #include "superhyperboloid.h"
    #include "supertoroid.h"
    #include "blobs.h"
    #include "rain.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
            if (isnan(t)) {
                t = 0;
            }
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

//        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/superhyperboloid-two.txt", ios_base::out);
        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/rain.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        hitable *list[5];
/*
        list[0] = new rain(vec3(0, 3, 0), 2, 2, 2, new lambertian(vec3(1.0, 0.0, 0.0)));
        list[1] = new rain(vec3(-3, 4, -4), 2, 2, 2, new lambertian(vec3(0.0, 1.0, 0.0)));
        list[2] = new rain(vec3(3, 4, -4), 2, 2, 2, new lambertian(vec3(0.0, 1.0, 0.0)));
        list[3] = new rain(vec3(-7, 5, -6), 2, 2, 2, new lambertian(vec3(0.0, 0.0, 1.0)));
        list[4] = new rain(vec3(7, 5, -6), 2, 2, 2, new lambertian(vec3(0.0, 0.0, 1.0)));
*/
/*
        list[0] = new rain(vec3(0, 3, 0), 2, 2, 2, new dielectric(9));
        list[1] = new rain(vec3(-3, 4, -4), 2, 2, 2, new dielectric(9));
        list[2] = new rain(vec3(3, 4, -4), 2, 2, 2, new dielectric(9));
        list[3] = new rain(vec3(-6.5, 4.5, -6), 2, 2, 2, new dielectric(9));
        list[4] = new rain(vec3(6.5, 4.5, -6), 2, 2, 2, new dielectric(9));
*/
        list[0] = new rain(vec3(0, 3, 0), 2, 2, 2, new lambertian(vec3(1.0, 0.0, 0.0)));
        list[1] = new rain(vec3(-3, 4, -4), 2, 2, 2, new lambertian(vec3(0.0, 1.0, 0.0)));
        list[2] = new rain(vec3(3, 4, -4), 2, 2, 2, new lambertian(vec3(0.0, 1.0, 0.0)));
        list[3] = new rain(vec3(-6.5, 4.5, -6), 2, 2, 2, new lambertian(vec3(0.0, 0.0, 1.0)));
        list[4] = new rain(vec3(6.5, 4.5, -6), 2, 2, 2, new lambertian(vec3(0.0, 0.0, 1.0)));

        hitable *world = new hitable_list(list,5);

        vec3 lookfrom(0, 0, 20);
        vec3 lookat(0, 3, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }

#elif testNumber == 36 /*36: parametric surface*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "parametric_surface.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
            if (isnan(t)) {
                t = 0;
            }
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 1;

//        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/superhyperboloid-two.txt", ios_base::out);
        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/parametric_surface_bezier3_teapot.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

/*
//1
        vec3 ctrl_points[16] = {vec3(-3, 2,  0), vec3(-2, 3,  0), vec3(2, 3,  0), vec3(3, 2,  0),
                                  vec3(-3, 5, -2), vec3(-2, 5, -2), vec3(2, 5, -1), vec3(3, 5, -2),
                                  vec3(-3, 5, -4), vec3(-2, 5, -4), vec3(2, 5, -3), vec3(3, 5, -4),
                                  vec3(-3, 2, -6), vec3(-2, 3, -6), vec3(2, 3, -4), vec3(3, 2, -6)};
//2
        vec3 ctrl_points[16] = {vec3(-3, 2,  0), vec3(-2, 3,  0), vec3(2, 3,  0), vec3(3, 2,  0),
                                  vec3(-3, 5, -2), vec3(-2, 5, -2), vec3(2, 4, -1), vec3(3, 3, -2),
                                  vec3(-3, 5, -4), vec3(-2, 5, -4), vec3(2, 4, -3), vec3(3, 3, -4),
                                  vec3(-3, 2, -6), vec3(-2, 3, -6), vec3(2, 3, -4), vec3(3, 2, -6)};
//3
        vec3 ctrl_points[16] = {vec3(-3, 2,  0), vec3(-2, 3,  0), vec3(2, 3,  0), vec3(3, 2,  0),
                                  vec3(-3, 5, -2), vec3(-2, 5, -2), vec3(2, 2, -1), vec3(3, 1, -2),
                                  vec3(-3, 5, -4), vec3(-2, 5, -4), vec3(2, 2, -3), vec3(3, 1, -4),
                                  vec3(-3, 2, -6), vec3(-2, 3, -6), vec3(2, 3, -4), vec3(3, 2, -6)};
//4
        vec3 ctrl_points[16] = {vec3(-3, 2,  0), vec3(-2, 3,  0), vec3(2, 5,  0), vec3(3, 5,  0),
                                  vec3(-3, 5, -2), vec3(-2, 5, -2), vec3(2, 2, -1), vec3(3, 1, -2),
                                  vec3(-3, 5, -4), vec3(-2, 5, -4), vec3(2, 2, -3), vec3(3, 1, -4),
                                  vec3(-3, 2, -6), vec3(-2, 3, -6), vec3(2, 5, -4), vec3(3, 5, -6)};
//5
        vec3 ctrl_points[16] = {vec3(-3, 2,  0), vec3(-2, 3,  0), vec3(2, 5,  0), vec3(3, 5,  0),
                                  vec3(-3, 5, -2), vec3(-2, 5, -2), vec3(2, 0, -1), vec3(3, -1, -2),
                                  vec3(-3, 5, -4), vec3(-2, 5, -4), vec3(2, 0, -3), vec3(3, -1, -4),
                                  vec3(-3, 2, -6), vec3(-2, 3, -6), vec3(2, 5, -4), vec3(3, 5, -6)};
//6
        vec3 ctrl_points[16] = {vec3(-3, 2,  0), vec3(-2, 3,  0), vec3(2, 5,  0), vec3(3, 5,  0),
                                  vec3(-3, -1, -2), vec3(-2, 5, -2), vec3(2, 0, -1), vec3(3, -1, -2),
                                  vec3(-3, -1, -4), vec3(-2, 5, -4), vec3(2, 0, -3), vec3(3, -1, -4),
                                  vec3(-3, 2, -6), vec3(-2, 3, -6), vec3(2, 5, -4), vec3(3, 5, -6)};

        int patches[32][16];
        float vertices[306][3];
        get_teapot_data(patches, vertices);
        vec3 patches_vertices[32][16];
        for (int i=0; i<32; i++) {
            for (int j=0; j<16; j++) {
                patches_vertices[i][j] = vec3(vertices[(patches[i][j])][0], vertices[(patches[i][j])][2], vertices[(patches[i][j])][1]);
            }
        }


        hitable *list[32];
//        list[0] = new parametric_surface(vec3(0.0, 3.0, 0.0), 3.0, 3.0, 3.0, 8.0, new lambertian(vec3(1.0, 0.0, 0.0)), 0.05, 20, 200);
//        list[0] = new parametric_surface(new lambertian(vec3(1.0, 0.0, 0.0)), 0.05, 20, 200, ctrl_points);
//        list[0] = new parametric_surface(new lambertian(vec3(1.0, 0.0, 0.0)), 0.05, 20, 200, rim_v);
        for (int i=0; i<32; i++) {
            list[i] = new parametric_surface(new lambertian(vec3(1.0, 0.0, 0.0)), 0.05, 20, 200, patches_vertices[i]);
        }
        hitable *world = new hitable_list(list,32);

        vec3 lookfrom(0.75, 1.5, 20);
        vec3 lookat(0.75, 1.5, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
*/
//11
        vec3 ctrl_points[16] = {vec3(-3, 2,  0), vec3(-2, 3,  0), vec3(2, 5,  0), vec3(-6, 3,  0),
                                  vec3(-3, -1, -2), vec3(-2, 5, -2), vec3(2, 0, -1), vec3(-6,  3, -2),
                                  vec3(-3, -1, -4), vec3(-2, 5, -4), vec3(2, 0, -3), vec3(-6,  3, -4),
                                  vec3(-3, 2, -6), vec3(-2, 3, -6), vec3(2, 5, -4), vec3(-6, 3, -6)};
        hitable *list[1];
        list[0] = new parametric_surface(new lambertian(vec3(1.0, 0.0, 0.0)), 0.05, 20, 200, ctrl_points);
        hitable *world = new hitable_list(list,1);

        vec3 lookfrom(10, 10, 10);
        vec3 lookat(0, 3, -2);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }

#elif testNumber == 37 /*37: matrix*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>
    #include <vec3.h>

    using namespace std;

    int main(){
/*
        vec3 v1 = vec3(1, 2, -3);
        vec3 v2 = vec3(0, 1, 2);
        vec3 v3 = vec3(1, 0, -5);

        float m1[3][1], m2[3][1], m3[3][1], m33[3][3], in[3][3];

        get_matrix_3_1(v1, m1);
        get_matrix_3_1(v2, m2);
        get_matrix_3_1(v3, m3);
        get_matrix_3_3(v1, v2, v3, m33);

        get_matrix_inverse_3_3(m33, in);

        vec3 v4 = vec3(2, -1, 0);
        vec3 v5 = vec3(1, 1, 3);
        vec3 v6 = vec3(1, 2, 1);
        vec3 v7 = vec3(1, -1, 2);

        float m4[3][1], m5[3][3], m6[3][1];

        get_matrix_3_1(v7, m4);
        get_matrix_3_3(v4, v5, v6, m5);
        matrix_3_3_multiply_3_1(m5, m4, m6);
*/
//        float matrix1[4][4] = {{1, 2, 3, 4}, {2, 3, 4, 5}, {3, 4, 5, 6}, {4, 5, 6, 7}};
//        float matrix2[4][4] = {{2, 3, 4, 5}, {3, 4, 5, 6}, {4, 5, 6, 7}, {5, 6, 7, 8}};
        float matrix1[4][4] = {{-3, -1, 1, 3}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
        float matrix2[4][4] = {{1, -3, 3, -1}, {0, 3, -6, 3}, {0, 0, 3, -3}, {0, 0, 0, 1}};
        float matrix3[4][4] = {{1, 2, 3, 4}, {3, 4, 5, 6}, {5, 6, 7, 8}, {7, 8, 9, 0}};
        float matrix4[3][3] = {{1, 0, 1}, {2, 1, 0}, {-3, 2, -5}};
        float result[4][4];
        float transpose[4][4];
        float inverse[3][3];
        matrix_4_4_multiply_4_4(matrix1, matrix2, result);
//        get_matrix_transpose_4_4(matrix3, transpose);
//        get_matrix_inverse_3_3(matrix4, inverse);


        std::cout << "matrix1:" << endl;
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                std::cout << matrix1[i][j] << "   ";
            }
            std::cout << endl;
        }
        std::cout << endl;
        std::cout << "matrix2:" << endl;
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                std::cout << matrix2[i][j] << "   ";
            }
            std::cout << endl;
        }
        std::cout << endl;
        std::cout << "matrix1 multiply matrix2:" << endl;
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                std::cout << result[i][j] << "   ";
            }
            std::cout << endl;
        }
/*
        std::cout << "matrix3:" << endl;
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                std::cout << matrix3[i][j] << "   ";
            }
            std::cout << endl;
        }
        std::cout << endl;
        std::cout << "transpose:" << endl;
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                std::cout << transpose[i][j] << "   ";
            }
            std::cout << endl;
        }

        std::cout << "matrix4:" << endl;
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                std::cout << matrix4[i][j] << "    ";
            }
            std::cout << endl;
        }
        std::cout << endl;
        std::cout << "inverse:" << endl;
        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                std::cout << inverse[i][j] << "    ";
            }
            std::cout << endl;
        }
*/
    }

#elif testNumber == 38 /*38: read file*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>
    #include <string.h>
    #include <vec3.h>

    using namespace std;

    int main(){
        int patches[32][16];
        float vertices[306][3];
        get_teapot_data(patches, vertices);

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/log.txt", ios_base::out);
        for (int i=0; i<32; i++) {
            outfile << "patches[" << i << "]:";
            for (int j=0; j<16; j++) {
                outfile <<  patches[i][j] << "   ";
            }
            outfile << endl;
        }
        for (int i=0; i<306; i++) {
            outfile << "vertices[" << i << "]:";
            for (int j=0; j<3; j++) {
                outfile <<  vertices[i][j] << "   ";
            }
            outfile << endl;
        }
        outfile.close();
        return 0;
    }
#elif testNumber == 39 /*39: translational sweeping*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "parametric_surface.h"
    #include "sweeping_translational.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
            if (isnan(t)) {
                t = 0;
            }
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

//        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/superhyperboloid-two.txt", ios_base::out);
        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/sweeping_translational.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        vec3 ctrl_points[17] = {vec3(-4.0,  0.0, -1.0), vec3(-3.0,  0.0, -4.0), vec3(-2.0,  0.0, -4.0), vec3(-1.0,  0.0,  -1.0),
                                vec3( 1.0,  0.0, -1.0), vec3( 1.0,  0.0, -5.0), vec3( 3.0,  0.0, -5.0), vec3( 4.0,  0.0,   1.0),
                                vec3( 2.0,  0.0,  2.0), vec3( 0.5,  0.0,  2.0), vec3(-2.0,  0.0,  1.0), vec3(-2.0,  0.0,  -2.0),
                                vec3(-2.5,  0.0, -3.0), vec3(-3.0,  0.0, -1.0), vec3(-4.0,  0.0, -1.0), vec3(-3.0,  0.0,  -4.0), vec3(-2.0,  0.0, -4.0)};
        hitable *list[1];
        list[0] = new sweeping_translational(ctrl_points, -2.0, 3.0, new lambertian(vec3(1.0, 0.0, 0.0)), 2);
        hitable *world = new hitable_list(list,1);

        vec3 lookfrom(0.0001, 10, 20);
        vec3 lookat(0, 1, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
/*
        list[0] = new sweeping_translational(ctrl_points, 0.0, 4.0, new lambertian(vec3(1.0, 0.0, 0.0)), 2);
        hitable *world = new hitable_list(list,1);

        vec3 lookfrom(0.0001, 12, 20);
        vec3 lookat(0, 2, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
*/
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }

#elif testNumber == 40 /*40: roots number of equation of 6th degree or 10th degree*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>
    #include <vec3.h>

    using namespace std;

    int main(){

//equation of 6th degree
//        float ee6[7] = {1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};//x^6------------1
//        float ee6[7] = {1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0};//x^5*(x-1)------------2
//        float ee6[7] = {1.0, 0.0,  -1.0,  0.0,  0.0,  0.0,  0.0};//x^4*(x-1)*(x+1)------------3
//        float ee6[7] = {1.0,  0.0,  -5.0,  0.0,  4.0,  0.0,  0.0};//x^2*(x-1)*(x+1)*(x-2)*(x+2)------------5
//        float ee6[7] = {1.0,  3.0,  -5.0, -15.0,  4.0,  12.0,  0.0};//x*(x-1)*(x+1)*(x-2)*(x+2)*(x+3)------------6
        float ee6[7] = {1.0,  0.0, -0.4236111, 0.0,  0.050347222,  0.0, -0.0017361111};//(x-1/2)*(x+1/2)*(x-1/3)*(x+1/3)*(x-1/4)*(x+1/4)------------6
        float a =  0.5;
        float b =  0.6;
        int num;
        float tol = 1e-6;
        float roots[7];
        roots_num_equation_6th(ee6, a, b, num);
        std::cout << "the coefficients of the equation of 6th degree are:" << endl;
        for (int i=0; i<7; i++) {
            std::cout << ee6[i] << "  ";
        }
        std::cout << endl;
        std::cout << "the corresponding equation is:" << endl;
        std::cout << "(x-1/2)*(x+1/2)*(x-1/3)*(x+1/3)*(x-1/4)*(x+1/4)=0" << endl;
        std::cout << endl;
        std::cout << "the number of roots of this equation in the interval [" << a << ", " << b << "] is: " << num << endl;

        if (num > 0) {
            roots_equation_6th(ee6, a, b, tol, roots);
            std::cout << "the roots are: ";
            for (int j=1; j<(num+1); j++){
                std::cout << roots[j] << " ";
            }
            std::cout << endl;
        }
/*
//equation of 10th degree
//        float ee10[11] = {1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};//x^10------------1
//        float ee10[11] = {1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};//x^9*(x-1)------------2
//        float ee10[11] = {1.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};//x^8*(x+1)*(x-1)------------3
//        float ee10[11] = {1.0,  0.0, -5.0,  0.0,  4.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};//x^6*(x+1)*(x-1)*(x+2)*(x-2)------------5
//        float ee10[11] = {1.0,  0.0, -14.0,  0.0,  49.0,  0.0,  -36.0,  0.0,  0.0,  0.0,  0.0};//x^4*(x+1)*(x-1)*(x+2)*(x-2)*(x+3)*(x-3)------------7
        float ee10[11] = {1.0,  0.0, -30.0,  0.0,  273.0,  0.0,  -820.0,  0.0,  576.0,  0.0,  0.0};//x^2*(x+1)*(x-1)*(x+2)*(x-2)*(x+3)*(x-3)*(x+4)*(x-4)------------9
//        float ee10[11] = {1.0, -5.0, -30.0,  150.0,  273.0, -1365.0,  -820.0, 4100.0,  576.0, -2880.0,  0.0};//x*(x+1)*(x-1)*(x+2)*(x-2)*(x+3)*(x-3)*(x+4)*(x-4)*(x-5)------------10
        float a = -10;
        float b =  0;
        int num;
        float tol = 1e-6;
        float roots[11];
        roots_num_equation_10th(ee10, a, b, num);
        std::cout << "the coefficients of the equation of 10th degree are:" << endl;
        for (int i=0; i<11; i++) {
            std::cout << ee10[i] << "  ";
        }
        std::cout << endl;
        std::cout << "the corresponding equation is:" << endl;
        std::cout << "x^2*(x+1)*(x-1)*(x+2)*(x-2)*(x+3)*(x-3)*(x+4)*(x-4)" << endl;
        std::cout << endl;
        std::cout << "the number of roots of this equation in the interval [" << a << ", " << b << "] is: " << num << endl;

        if (num > 0) {
            roots_equation_10th(ee10, a, b, tol, roots);
            std::cout << "the roots are: ";
            for (int j=1; j<(num+1); j++){
                std::cout << roots[j] << " ";
            }
            std::cout << endl;
        }
*/
    }

#elif testNumber == 41 /*41: rotational sweeping*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "parametric_surface.h"
    #include "sweeping_rotational.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
            if (isnan(t)) {
                t = 0;
            }
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

//        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/superhyperboloid-two.txt", ios_base::out);
        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/sweeping_rotational.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";
/*
//1
        vec3 ctrl_points[6] = {vec3(-1.0,  5.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//2
        vec3 ctrl_points[6] = {vec3(-4.0,  5.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//3
        vec3 ctrl_points[6] = {vec3( 2.0,  5.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//4
        vec3 ctrl_points[6] = {vec3( 4.0,  5.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//5
        vec3 ctrl_points[6] = {vec3( 2.0,  7.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//6
        vec3 ctrl_points[6] = {vec3(-4.0,  2.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//7
        vec3 ctrl_points[6] = {vec3(-4.0,  0.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//8--discover some problems
        vec3 ctrl_points[6] = {vec3( 4.0,  4.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//9
        vec3 ctrl_points[6] = {vec3( 4.0,  2.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//10
        vec3 ctrl_points[6] = {vec3( 4.0,  0.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//11
        vec3 ctrl_points[6] = {vec3( 3.0, -5.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//12
        vec3 ctrl_points[6] = {vec3(-2.0, -3.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
*/
//1
        vec3 ctrl_points1[6] = {vec3(-1.0,  5.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-2.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//2
        vec3 ctrl_points2[6] = {vec3(-4.0,  5.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//3
        vec3 ctrl_points3[6] = {vec3( 2.0,  5.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//4
        vec3 ctrl_points4[6] = {vec3( 4.0,  5.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//5
        vec3 ctrl_points5[6] = {vec3( 2.0,  7.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//6
        vec3 ctrl_points6[6] = {vec3(-4.0,  2.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//7
        vec3 ctrl_points7[6] = {vec3(-4.0,  0.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//8--discover some problems
        vec3 ctrl_points8[6] = {vec3( 4.0,  4.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//9
        vec3 ctrl_points9[6] = {vec3( 4.0,  2.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//10
        vec3 ctrl_points10[6] = {vec3( 4.0,  0.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//11
        vec3 ctrl_points11[6] = {vec3( 3.0, -5.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//12
        vec3 ctrl_points12[6] = {vec3(-2.0, -3.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};

        hitable *list[1];
        list[0] = new sweeping_rotational(vec3( 0.0, 0,  0.0), ctrl_points1, new lambertian(vec3(1.0, 0.0, 0.0)));
//        list[1] = new sweeping_rotational(vec3(-5.0, 0,  0.0), ctrl_points, new lambertian(vec3(1.0, 0.0, 0.0)));
//        list[2] = new sweeping_rotational(vec3( 5.0, 0,  0.0), ctrl_points, new lambertian(vec3(1.0, 0.0, 0.0)));
//        list[1] = new sweeping_rotational(vec3(-5.0, 0,  0.0), ctrl_points, new lambertian(vec3(1.0, 0.0, 0.0)));
//        list[2] = new sweeping_rotational(vec3( 5.0, 0,  0.0), ctrl_points, new lambertian(vec3(1.0, 0.0, 1.0)));
//        list[1] = new sphere(vec3(0,-103,-1), 100, new metal(vec3(0.0, 1.0, 0.5), 0.0), 0);
        hitable *world = new hitable_list(list,1);

        vec3 lookfrom(0, 5, 14);
        vec3 lookat(0, 1, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 30, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
/*
        list[0] = new sweeping_translational(ctrl_points, 0.0, 4.0, new lambertian(vec3(1.0, 0.0, 0.0)), 2);
        hitable *world = new hitable_list(list,1);

        vec3 lookfrom(0.0001, 12, 20);
        vec3 lookat(0, 2, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
*/
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }
#elif testNumber == 42 /*42: sphere sweeping*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "triangle.h"
    #include "polygon.h"
    #include "sphere.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "parametric_surface.h"
    #include "sweeping_rotational.h"
    #include "sweeping_sphere.h"

    using namespace std;

    extern int dielectric_counter;
    extern int refl_counter;
    extern int refr_counter;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
//                std::cout << "-------------color()----------------depth:" << depth << endl;
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
            if (isnan(t)) {
                t = 0;
            }
//            std::cout << "-------------color()--5----------------" << endl;
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }


    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

//        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/superhyperboloid-two.txt", ios_base::out);
        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/sweeping_sphere.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";
/*
//0
        vec3 ctrl_points1[6] = {vec3(-1.0,  5.0,  0.0), vec3( 2.0,  4.0,  0.0),
                               vec3( 2.0,  1.0,  0.0), vec3(-0.5,  1.0,  0.0),
                               vec3( 1.5, -3.0,  0.0), vec3( 3.0,  0.0,  0.0)};
//4---error
        vec3 ctrl_points_c[6] = {vec3(10.0,  3.0,  0.0), vec3(-5.0,  5.0,  0.0),
                                 vec3(-5.0, -5.0,  0.0), vec3( 1.0, -5.0,  0.0),
                                 vec3( 1.0,  0.0,  0.0), vec3(-8.0,  0.0,  0.0)};
//5
        vec3 ctrl_points_c[6] = {vec3(10.0,  3.0,  0.0), vec3(-5.0,  5.0,  0.0),
                                 vec3(-5.0, -5.0,  0.0), vec3( 1.0, -5.0,  0.0),
                                 vec3( 1.0,  0.0,  0.0), vec3(-1.0,  0.0,  0.0)};
//22,4,4,40
        vec3 ctrl_points_c[8] = {vec3(-4.0, -4.0,  8.0), vec3(-4.0,  0.0, -4.0),
                                 vec3(-4.0,  8.0, -4.0), vec3(-4.0,  8.0,  4.0),
                                 vec3( 4.0,  8.0,  4.0), vec3( 4.0,  8.0, -4.0),
                                 vec3( 4.0,  0.0, -4.0), vec3( 4.0, -4.0,  8.0)};
        vec3 ctrl_points_r[8] = {vec3( 2.0,  0.0,  0.0), vec3( 2.0,  0.0,  0.0),
                                 vec3( 1.0,  0.0,  0.0), vec3( 0.5,  0.0,  0.0),
                                 vec3( 0.5,  0.0,  0.0), vec3( 1.0,  0.0,  0.0),
                                 vec3( 2.0,  0.0,  0.0), vec3( 2.0,  0.0,  0.0)};
//1
        vec3 ctrl_points_r[8] = {vec3( 0.0,  4.0,  0.0), vec3( 1.0,  3.5,  0.0),
                                 vec3( 1.0,  3.0,  0.0), vec3( 0.2,  2.5,  0.0),
                                 vec3( 3.0,  1.5,  0.0), vec3( 3.0, -1.5,  0.0),
                                 vec3( 2.0, -4.0,  0.0), vec3( 0.0, -4.0,  0.0)};
*/
//rotational sweeping
        vec3 ctrl_points_r[8] = {vec3(-4.0,  5.0,  0.0), vec3( 1.0,  3.5,  0.0),
                                 vec3( 1.0,  3.0,  0.0), vec3( 0.2,  2.5,  0.0),
                                 vec3( 3.0,  1.5,  0.0), vec3( 3.0, -1.5,  0.0),
                                 vec3( 2.0, -4.0,  0.0), vec3(-4.0, -4.0,  0.0)};

//sphere sweeping, teapot handle
        vec3 ctrl_points_sc_h[8] = {vec3( 1.0,  0.5,  0.0), vec3( 2.5,  0.5,  0.0),
                                    vec3( 4.0,  1.0,  0.0), vec3( 5.0,  0.8,  0.0),
                                    vec3( 4.0, -1.5,  0.0), vec3( 3.0, -2.0,  0.0),
                                    vec3( 2.5, -2.5,  0.0), vec3( 1.5, -1.5,  0.0)};
        vec3 ctrl_points_sr_h[8] = {vec3( 0.4,  0.0,  0.0), vec3( 0.4,  0.0,  0.0),
                                    vec3( 0.4,  0.0,  0.0), vec3( 0.4,  0.0,  0.0),
                                    vec3( 0.4,  0.0,  0.0), vec3( 0.4,  0.0,  0.0),
                                    vec3( 0.4,  0.0,  0.0), vec3( 0.4,  0.0,  0.0)};
//sphere sweeping, teapot spout
        vec3 ctrl_points_sc_s[8] = {vec3(-0.5, -1.0,  0.0), vec3(-2.5, -1.5,  0.0),
                                    vec3(-4.0, -1.5,  0.0), vec3(-4.5, -0.5,  0.0),
                                    vec3(-5.0,  0.5,  0.0), vec3(-6.0,  1.5,  0.0),
                                    vec3(-7.0,  1.5,  0.0), vec3(-8.0,  0.0,  0.0)};
        vec3 ctrl_points_sr_s[8] = {vec3( 1.0,  0.0,  0.0), vec3( 1.0,  0.0,  0.0),
                                    vec3( 0.8,  0.0,  0.0), vec3( 0.8,  0.0,  0.0),
                                    vec3( 0.4,  0.0,  0.0), vec3( 0.4,  0.0,  0.0),
                                    vec3( 0.2,  0.0,  0.0), vec3( 0.2,  0.0,  0.0)};
//22,4,4,40
        vec3 ctrl_points_sc[8] = {vec3(-4.0, -4.0,  8.0), vec3(-4.0,  0.0, -4.0),
                                 vec3(-4.0,  8.0, -4.0), vec3(-4.0,  8.0,  4.0),
                                 vec3( 4.0,  8.0,  4.0), vec3( 4.0,  8.0, -4.0),
                                 vec3( 4.0,  0.0, -4.0), vec3( 4.0, -4.0,  8.0)};
        vec3 ctrl_points_sr[8] = {vec3( 2.0,  0.0,  0.0), vec3( 2.0,  0.0,  0.0),
                                 vec3( 1.0,  0.0,  0.0), vec3( 0.5,  0.0,  0.0),
                                 vec3( 0.5,  0.0,  0.0), vec3( 1.0,  0.0,  0.0),
                                 vec3( 2.0,  0.0,  0.0), vec3( 2.0,  0.0,  0.0)};

        hitable *list[1];
//        list[0] = new sweeping_rotational(vec3( 0.0, 0.0,  0.0), ctrl_points_r, new lambertian(vec3(1.0, 0.0, 0.0)), 0.0);
//        list[1] = new sweeping_sphere(ctrl_points_sc_h, ctrl_points_sr_h, new lambertian(vec3(1.0, 0.0, 0.0)));
//        list[2] = new sweeping_sphere(ctrl_points_sc_s, ctrl_points_sr_s, new lambertian(vec3(1.0, 0.0, 0.0)));
        list[0] = new sweeping_sphere(ctrl_points_sc, ctrl_points_sr, new lambertian(vec3(1.0, 0.0, 0.0)));
        hitable *world = new hitable_list(list,1);

        vec3 lookfrom(10, 20, 10);
        vec3 lookat(0, 4, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 40, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
/*
        list[0] = new sweeping_translational(ctrl_points, 0.0, 4.0, new lambertian(vec3(1.0, 0.0, 0.0)), 2);
        hitable *world = new hitable_list(list,1);

        vec3 lookfrom(0.0001, 12, 20);
        vec3 lookat(0, 2, 0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 20, float(nx)/float(ny), aperture, 0.7*dist_to_focus);
*/
        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    //generate a random in range [0,1]
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
//                    std::cout << "-------------main()--1----------------" << endl;
                    ray r = cam.get_ray(u, v);
//                    std::cout << "-------------main()--2----------------" << endl;

//                    vec3 p = r.point_at_parameter(2.0);
//                    std::cout << "(i,j): (" << i << "," << j << ")==" << "in call color()" << endl;
                    col += color(r, world, 0);
//                    std::cout << "-------------main()--3----------------" << endl;
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
     std::cout << "dielectric_counter: " << dielectric_counter << " " << "refl_counter: " << refl_counter << " " << "refr_counter: " << refr_counter << endl;

    }

#elif testNumber == 43 /*43: binary tree*/
     #include<iostream>
     using namespace std;

     struct tree
     {
         int data;
         tree *left,*right;
     };
     class Btree
     {
     public:
         Btree()
         {
             root=NULL;
         }
         void create_Btree(int);
         void Preorder(tree *);
         void inorder(tree *);
         void Postorder(tree *);
         void display1() {Preorder(root); std::cout<<endl;}
         void display2() {inorder(root);std::cout<<endl;}
         void display3() {Postorder(root); std::cout<<endl;}
         int count(tree *);
         int findleaf(tree *);
         int findnode(tree *);
         static int n;
         static int m;
         tree *root;
     };
     int Btree::n=0;
     int Btree::m=0;
     void Btree::create_Btree(int x)
     {
         tree *newnode=new tree;
         newnode->data=x;
         newnode->right=newnode->left=NULL;
         if(root==NULL)
             root=newnode;
         else
         {
             tree *back;
             tree *current=root;
             while(current!=NULL)
             {
                 back=current;
                 if(current->data>x)
                     current=current->left;
                 else
                     current=current->right;
             }
             if(back->data>x)
                 back->left=newnode;
             else
                 back->right=newnode;
         }
     }
     int Btree::count(tree *p)
     {
         if(p==NULL)
             return 0;
         else
             return count(p->left)+count(p->right)+1;
     }
     void Btree::Preorder(tree *temp)
     {
         if(temp!=NULL)
         {
             std::cout<<temp->data<<" ";
             Preorder(temp->left);
             Preorder(temp->right);
         }
     }
     void Btree::inorder(tree *temp)
     {
         if(temp!=NULL)
         {
             inorder(temp->left);
             std::cout<<temp->data<<" ";
             inorder(temp->right);
         }
     }
     void Btree::Postorder(tree *temp)
     {
         if(temp!=NULL)
         {
             Postorder(temp->left);
             Postorder(temp->right);
             std::cout<<temp->data<<" ";
         }
     }
     int Btree::findleaf(tree *temp)
     {
         if(temp==NULL)return 0;
         else
         {
             if(temp->left==NULL&&temp->right==NULL)return n+=1;
             else
             {
                 findleaf(temp->left);
                 findleaf(temp->right);
             }
             return n;
         }
     }
     int Btree::findnode(tree *temp)
     {
         if(temp==NULL)return 0;
         else
         {
             if(temp->left!=NULL&&temp->right!=NULL)
             {
                 findnode(temp->left);
                 findnode(temp->right);
             }
             if(temp->left!=NULL&&temp->right==NULL)
             {
                 m+=1;
                 findnode(temp->left);
             }
             if(temp->left==NULL&&temp->right!=NULL)
             {
                 m+=1;
                 findnode(temp->right);
             }
         }
         return m;
     }


     int main()
     {
         Btree A;
         int array[]={7,4,2,3,15,35,6,45,55,20,1,14,56,57,58};
         int k;
         k=sizeof(array)/sizeof(array[0]);
         std::cout<<"create binary tree: "<<endl;
         for(int i=0;i<k;i++)
         {
             std::cout<<array[i]<<" ";
             A.create_Btree(array[i]);
         }
         std::cout<<endl;
         std::cout<<"the number of nodes: "<<A.count(A.root)<<endl;
         std::cout<<"the number of leaves: "<<A.findleaf(A.root)<<endl;
         std::cout<<"the number of degree 1 nodes: "<<A.findnode(A.root)<<endl;
         std::cout<<endl<<"preorder traverse: "<<endl;
         A.display1();
         std::cout<<endl<<"inorder traverse: "<<endl;
         A.display2();
         std::cout<<endl<<"postorder traverse: "<<endl;
         A.display3();
         return 0;
     }
#elif testNumber == 44 /*44: binary tree 2*/
    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include <TreeType.h>

    using namespace std;

    int main()
    {
         TreeType binarySearchTree;
         ItemType infoArray[15] = {20, 12, 2, 14, 15, 35, 6, 45, 25, 23, 1, 24, 16, 57, 28};
         std::cout << "CreateTree: " << endl;
         binarySearchTree.CreateTree(infoArray, 15);
         std::cout << endl;

         std::cout << "PreorderTraverse: " << endl;
         binarySearchTree.PreorderTraverse(binarySearchTree.root);
         std::cout << endl;

         std::cout << "InorderTraverse: " << endl;
         binarySearchTree.InorderTraverse(binarySearchTree.root);
         std::cout << endl;

         std::cout << "PostorderTraverse: " << endl;
         binarySearchTree.PostorderTraverse(binarySearchTree.root);
         std::cout << endl;
         std::cout << endl;

         ItemType insertItems[4] = {40, 32, 42, 28};
         std::cout << "InsertItems: ";
         for (int i=0; i<4; i++) {
            std::cout << insertItems[i] << " ";
         }
         std::cout << endl;
         binarySearchTree.InsertItems(insertItems, 4);
         std::cout << endl;

         std::cout << "PreorderTraverse: " << endl;
         binarySearchTree.PreorderTraverse(binarySearchTree.root);
         std::cout << endl;

         std::cout << "InorderTraverse: " << endl;
         binarySearchTree.InorderTraverse(binarySearchTree.root);
         std::cout << endl;

         std::cout << "PostorderTraverse: " << endl;
         binarySearchTree.PostorderTraverse(binarySearchTree.root);
         std::cout << endl;

         std::cout << endl;
         ItemType item = 35;
         bool found;
         std::cout << "FindItem: " << item << endl;
         found = binarySearchTree.FindItem(binarySearchTree.root, item);
         std::cout << "found: " << found << endl;

         std::cout << endl;
         std::cout << "MinValue: " << endl;
         ItemType minItem = binarySearchTree.MinValue(binarySearchTree.root);
         std::cout << "minItem: " << minItem << endl;

         std::cout << endl;
         std::cout << "MaxValue: " << endl;
         ItemType maxItem = binarySearchTree.MaxValue(binarySearchTree.root);
         std::cout << "maxItem: " << maxItem << endl;

         std::cout << endl;
         std::cout << "NodesNum: " << endl;
         int num = binarySearchTree.NodesNum(binarySearchTree.root);
         std::cout << "NodesNum: " << num << endl;

         std::cout << endl;
         std::cout << "LeavesNum: " << endl;
         int lnum = 0;
         binarySearchTree.LeavesNum(binarySearchTree.root, lnum);
         std::cout << "LeavesNum: " << lnum << endl;

         std::cout << "PreorderTraverseLeaves: " << endl;
         binarySearchTree.PreorderTraverseLeaves(binarySearchTree.root);
         std::cout << endl;

         std::cout << "InorderTraverseLeaves: " << endl;
         binarySearchTree.InorderTraverseLeaves(binarySearchTree.root);
         std::cout << endl;

         std::cout << "PostorderTraverseLeaves: " << endl;
         binarySearchTree.PostorderTraverseLeaves(binarySearchTree.root);
         std::cout << endl;

         std::cout << endl;
         item = 45;
         std::cout << "DeleteItem: " << item << endl;
         bool deleted = binarySearchTree.DeleteItem(binarySearchTree.root, item);
         std::cout << "deleted: " << deleted << endl;
         std::cout << endl;

         std::cout << "PreorderTraverse: " << endl;
         binarySearchTree.PreorderTraverse(binarySearchTree.root);
         std::cout << endl;

         std::cout << "InorderTraverse: " << endl;
         binarySearchTree.InorderTraverse(binarySearchTree.root);
         std::cout << endl;

         std::cout << "PostorderTraverse: " << endl;
         binarySearchTree.PostorderTraverse(binarySearchTree.root);
         std::cout << endl;

         std::cout << endl;
         std::cout << "NodesNum: " << endl;
         num = binarySearchTree.NodesNum(binarySearchTree.root);
         std::cout << "NodesNum: " << num << endl;

         std::cout << endl;
         std::cout << "LeavesNum: " << endl;
         lnum = 0;
         binarySearchTree.LeavesNum(binarySearchTree.root, lnum);
         std::cout << "LeavesNum: " << lnum << endl;

         std::cout << "PreorderTraverseLeaves: " << endl;
         binarySearchTree.PreorderTraverseLeaves(binarySearchTree.root);
         std::cout << endl;

         std::cout << "InorderTraverseLeaves: " << endl;
         binarySearchTree.InorderTraverseLeaves(binarySearchTree.root);
         std::cout << endl;

         std::cout << "PostorderTraverseLeaves: " << endl;
         binarySearchTree.PostorderTraverseLeaves(binarySearchTree.root);
         std::cout << endl;
    }

#elif testNumber == 45 /*45: intervals*/
    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>
    #include <vec3.h>

    using namespace std;

    int main(){

//equation of 10th degree
//        float ee10[11] = {1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};//x^10------------1
//        float ee10[11] = {1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};//x^9*(x-1)------------2
//        float ee10[11] = {1.0,  0.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};//x^8*(x+1)*(x-1)------------3
//        float ee10[11] = {1.0,  0.0, -5.0,  0.0,  4.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0};//x^6*(x+1)*(x-1)*(x+2)*(x-2)------------5
//        float ee10[11] = {1.0,  0.0, -14.0,  0.0,  49.0,  0.0,  -36.0,  0.0,  0.0,  0.0,  0.0};//x^4*(x+1)*(x-1)*(x+2)*(x-2)*(x+3)*(x-3)------------7
//        float ee10[11] = {1.0,  0.0, -30.0,  0.0,  273.0,  0.0,  -820.0,  0.0,  576.0,  0.0,  0.0};//x^2*(x+1)*(x-1)*(x+2)*(x-2)*(x+3)*(x-3)*(x+4)*(x-4)------------9
        float ee10[11] = {1.0, -5.0, -30.0,  150.0,  273.0, -1365.0,  -820.0, 4100.0,  576.0, -2880.0,  0.0};//x*(x+1)*(x-1)*(x+2)*(x-2)*(x+3)*(x-3)*(x+4)*(x-4)*(x-5)------------10
        float a = -10;
        float b =  10;
        int num;
        float tol = 1e-6;
        float roots[11];
        roots_num_equation_10th(ee10, a, b, num);
        std::cout << "the coefficients of the equation of 10th degree are:" << endl;
        for (int i=0; i<11; i++) {
            std::cout << ee10[i] << "  ";
        }
        std::cout << endl;
        std::cout << "the corresponding equation is:" << endl;
        std::cout << "x^2*(x+1)*(x-1)*(x+2)*(x-2)*(x+3)*(x-3)*(x+4)*(x-4)" << endl;
        std::cout << endl;
        std::cout << "the number of roots of this equation in the interval [" << a << ", " << b << "] is: " << num << endl;

        if (num > 0) {
            roots_equation_10th_bt(ee10, a, b, tol, roots);
//            roots_equation_10th(ee10, a, b, tol, roots);
            std::cout << "the roots are: ";
            for (int j=1; j<(num+1); j++){
                std::cout << roots[j] << " ";
            }
            std::cout << endl;
        }
    }
#elif testNumber == 46 /*46: CSG*/

    #include <iostream>
    #include <fstream>
    #include <limits>
    #include <stdlib.h>

    #include "sphere.h"
    #include "box.h"
    #include "hitable_list.h"
    #include "float.h"
    #include "camera.h"
    #include "material.h"
    #include "lambertian.h"
    #include "metal.h"
    #include "dielectric.h"
    #include "csgTree.h"

    using namespace std;

    vec3 color(const ray& r, hitable *world, int depth) {
        hit_record rec;
        if (world->hit(r, 0.001, (numeric_limits<float>::max)(), rec)) {
            ray scattered;
            vec3 attenuation;
            if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
                return attenuation*color(scattered, world, depth+1);
            }
            else {
                return vec3(0,0,0);
            }
        }
        else {
            vec3 unit_direction = unit_vector(r.direction());
            float t = 0.5*(unit_direction.y() + 1.0);
            if (isnan(t)) {
                t = 0;
            }
            return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);//white, light blue
        }
    }

vec3 lookfrom;

    int main(){
        int nx = 200;
        int ny = 100;
        int ns = 100;

        ofstream outfile( "/Users/libingzeng/CG/AnIntroductionToRayTracing/results/csg.txt", ios_base::out);
        outfile << "P3\n" << nx << " " << ny << "\n255\n";

        std::cout << "P3\n" << nx << " " << ny << "\n255\n";

        csgTree *csTree = new csgTree();
        csgTreeNode *rootNode = new csgTreeNode;
        rootNode->info.operation = 4;
        rootNode->info.solid = NULL;

        hitable *list_csg[2];
//1
//        list_csg[0] = new sphere(vec3(2.5, 5.0, -2.5), 3.0, new lambertian(vec3(0.0, 0.0, 1.0)), 0, 1);
//2
        list_csg[0] = new sphere(vec3(0.0, 2.5, 0.0), 3.3, new lambertian(vec3(0.0, 0.0, 1.0)), 0, 1);
//        list[0] = new box(vec3(-2.5, 0.0, 2.5), vec3(2.5, 5.0, -2.5), new lambertian(vec3(1.0, 0.0, 0.0)), 1);
        list_csg[1] = new box(vec3(-2.5, 0.0, 2.5), vec3(2.5, 5.0, -2.5), new lambertian(vec3(1.0, 0.0, 0.0)), 1);
        csgTreeNode *node1 = new csgTreeNode;
        node1->info.operation = 0;
        node1->info.solid = list_csg[0];
        node1->left = NULL;
        node1->right = NULL;
        csgTreeNode *node2 = new csgTreeNode;
        node2->info.operation = 0;
        node2->info.solid = list_csg[1];
        node2->left = NULL;
        node2->right = NULL;

        rootNode->left = node1;
        rootNode->right = node2;
        csTree->root = rootNode;

        hitable *list[1];
        list[0] = csTree;
        hitable *world = new hitable_list(list,1);

        lookfrom = vec3(10, 10, 10);
        vec3 lookat(0.0, 2.5, 0.0);
        float dist_to_focus = (lookfrom - lookat).length();
        float aperture = 0.0;
        camera cam(lookfrom, lookat, vec3(0,1,0), 40, float(nx)/float(ny), aperture, 0.7*dist_to_focus);

        for (int j = ny-1; j >= 0; j--){
            for (int i = 0; i < nx; i++){
                vec3 col(0, 0, 0);
                for (int s = 0; s < ns; s++){
                    float random = rand()%(100)/(float)(100);
                    float u = float(i + random) / float(nx);
                    float v = float(j + random) / float(ny);
                    ray r = cam.get_ray(u, v);
                    col += color(r, world, 0);
                }
                col /= float(ns);
                col = vec3( sqrt(col[0]), sqrt(col[1]), sqrt(col[2]) );
                int ir = int (255.99*col[0]);
                int ig = int (255.99*col[1]);
                int ib = int (255.99*col[2]);

                outfile << ir << " " << ig << " " << ib << "\n";
                std::cout << ir << " " << ig << " " << ib << "\n";
            }
        }
    }

#endif // testNumber
