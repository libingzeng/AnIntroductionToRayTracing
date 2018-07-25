# AnIntroductionToRayTracing
This render is based on the framework of Peter Shirley's ray tracing in one weekend and is extended with tracing almost all of the surfaces mentioned in the book, An Introduction to Ray Tracing. The extended surfaces are all NON-TRIANGULATE, which include box, sphere, polygon, quadric surfaces, tori, blending and joining surface, superellipsoid, superhyperboloid, supertoroid, blobs, tear drop, bicubic Bezier patches, bicubic B-spline patches, translational sweeping surface, cone sweeping surface, rotational sweeping surface, sphere sweeping surface, and CSG surfaces.
Some of the main resultant images were uploaded in the path: results/pictures.


There is a Macro named testNumber in main.cpp.
We can set testNumber to run the corresponding test case and produce image.

#define testNumber 16
/*
1: output the first image
2: test "int &ri£¨int& ri£¨int *&pri"
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

