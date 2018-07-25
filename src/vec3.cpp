#include "vec3.h"
#include "float.h"

vec3 vector_trans(const vec3& v1, const vec3& u, const vec3& v, const vec3& w) {
/*translate vector v1 from normal space to u-v-w space*/
    int i,j;
    int h1=0;
    int h2=2;
    float k1,k2,k3;//the three unknowns
    float temp[3][4] = {0};
    float mn[3][4] = {
        {u.x(), v.x(), w.x(), v1.x()},
        {u.y(), v.y(), w.y(), v1.y()},
        {u.z(), v.z(), w.z(), v1.z()}};

    /*eliminate k1*/
    for (i=0; i<3; i++) {
        if(mn[i][0] != 0) {//choose the first row for temp
            for (j=0; j<4; j++) {
                temp[h1][j] = mn[i][j]/mn[i][0];//set the coefficient of k1 in the h1-th row of temp to 1
                if(h1 != 0) {
                    temp[h1][j] = temp[h1][j] - temp[0][j];//temp: the h1 row minus the first row
                }
            }
            h1++;
        }
        else {
            for (j=0; j<4; j++) {
                temp[h2][j] = mn[i][j];//copy the row of mn whose coefficient of k1 is 0 to the last row of temp
            }
            h2--;
        }
    }
    for (i=0; i<3; i++) {
        for (j=0; j<4; j++) {
            mn[i][j] = temp[i][j];
        }
    }
    h1 = 1;
    h2 = 2;

    /*eliminate k2*/
    for (i=1; i<3; i++) {
        if(temp[i][1] != 0) {
            for (j=1; j<4; j++) {
                mn[h1][j] = temp[i][j]/temp[i][1];
                if(h1 != 1) {
                    mn[h1][j] = mn[h1][j] - mn[1][j];
                }
            }
            h1++;
        }
        else {
            for (j=1; j<4; j++) {
                mn[h2][j] = temp[i][j];
            }
            h2--;
        }
    }

    k3 = mn[2][3] / mn[2][2];
    k2 = mn[1][3] - mn[1][2]*k3;
    k1 = mn[0][3] - mn[0][2]*k3 - mn[0][1]*k2;

    return vec3(k1, k2, k3);
}

vec3 vector_trans_back(const vec3& v1, const vec3& u, const vec3& v, const vec3& w) {
/*translate vector v1 from u-v-w space to normal space*/
    return vec3((v1.x()*u.x()+v1.y()*v.x()+v1.z()*w.x()),
                (v1.x()*u.y()+v1.y()*v.y()+v1.z()*w.y()),
                (v1.x()*u.z()+v1.y()*v.z()+v1.z()*w.z()));
}

vec3 get_vector_v(const vec3& vector_u, float angle) {
/*rectangle: determine the v axis of u-v-w space*/
/*if y coordinate of vector_u is zero, we regard the angle as theta, because in this case, there is only one certain value;*/
/*if y coordinate of vector_u is not zero, we regard the angle as phi, because in this case, the theta is limited*/
    if (vector_u.y() == 0) {
        if (angle == 0) {
            return vec3(0, 1, 0);
        }
        else {
            float theta = angle*M_PI/180;
            float A = vector_u.z();
            float B = vector_u.x();
            float tan_half_phi, sin_phi, cos_phi;
            tan_half_phi = (B-sqrt(B*B+A*A))/A;;
            sin_phi = 2*tan_half_phi/(1+tan_half_phi*tan_half_phi);
            if (sin_phi < 0) {
                tan_half_phi = (B-sqrt(B*B+A*A))/A;
                sin_phi = 2*tan_half_phi/(1+tan_half_phi*tan_half_phi);
            }

            cos_phi = (1-tan_half_phi*tan_half_phi)/(1+tan_half_phi*tan_half_phi);

            return (vec3(sin(theta)*sin_phi, cos(theta), sin(theta)*cos_phi));
        }

    }
    else {
        float phi = angle*M_PI/180;
        float A = vector_u.y();
        float B = vector_u.x()*sin(phi) + vector_u.z()*cos(phi);
        float tan_half_theta, sin_theta, cos_theta;
        tan_half_theta = (B+sqrt(B*B+A*A))/A;
        sin_theta = 2*tan_half_theta/(1+tan_half_theta*tan_half_theta);
        if (sin_theta < 0) {
            tan_half_theta = (B-sqrt(B*B+A*A))/A;
            sin_theta = 2*tan_half_theta/(1+tan_half_theta*tan_half_theta);
        }
        cos_theta = (1-tan_half_theta*tan_half_theta)/(1+tan_half_theta*tan_half_theta);

        return (vec3(sin_theta*sin(phi), cos_theta, sin_theta*cos(phi)));
    }
}

bool get_vector_vw2(vec3& u, float angle, vec3& v, vec3& w) {
/*determine the v axis and w axis of u-v-w space*/
    vec3 v1, w1;
    u = unit_vector(u);
    float theta = angle*M_PI/180;
    if ((u.x() == 0.0) && (u.z() == 0.0)) {
        v1 = vec3(1.0, 0.0, 0.0);
    }
    else {
        v1 = unit_vector(vec3(-u.z(), 0.0, u.x()));
//        theta = (angle-180)*M_PI/180;
//        theta = (angle-90)*M_PI/180;
        theta = (angle)*M_PI/180;
    }
    w1 = cross(v1, u);
    float x = cos(theta);
    float z = sin(theta);
    if (fabs(x) < 0.0001) {
        x = 0.0;
    }
    if (fabs(z) < 0.0001) {
        z = 0.0;
    }
    vec3 v_uv1w1 = vec3(x, 0, z);
    v = unit_vector(vector_trans_back(v_uv1w1, v1, u, w1));
    w = cross(v, u);
    return true;
}

bool get_vector_vw(vec3& u, float angle, vec3& v, vec3& w) {
/*determine the v axis and w axis of u-v-w space*/
    vec3 v1, w1;
    u = unit_vector(u);
    float theta = angle*M_PI/180;
    if (u.x() == 0.0) {
        v1 = vec3(1.0, 0.0, 0.0);
        theta = (angle)*M_PI/180;
    }
    else {
        v1 = unit_vector(vec3(-u.z(), 0.0, u.x()));
        vec3 z_n = vec3(0, 0, -1);
        vec3 u_xoz = vec3(u.x(), 0, u.z());
        float phi = acos(dot(z_n, u_xoz) / (z_n.length()*u_xoz.length()));
        if (u.x() > 0) {
            theta = (angle)*M_PI/180 - phi;
        }
        else {
            theta = (angle)*M_PI/180 + phi;
        }
    }
    w1 = cross(v1, u);
    float x = cos(theta);
    float z = sin(theta);
    if (fabs(x) < 0.0001) {
        x = 0.0;
    }
    if (fabs(z) < 0.0001) {
        z = 0.0;
    }
    vec3 v_uv1w1 = vec3(x, 0, z);
    v = unit_vector(vector_trans_back(v_uv1w1, v1, u, w1));
    w = cross(v, u);
    return true;
}

bool roots_quadratic_equation2(float a, float b, float c, float (&roots)[3]) {
    //the first element is the number of the real roots, and other elements are the real roots.
    if (a == 0.0) {
        if (b == 0.0) {
            roots[0] = 0.0;
        }
        else {
            roots[1] = -c/b;
            roots[0] = 1.0;
        }
    }
    else {
        float d = b*b - 4*a*c;
        if (d < 0.0) {
            roots[0] = 0.0;
        }
        else {
            roots[1] = (-b + sqrt(d)) / (2*a);
            roots[2] = (-b - sqrt(d)) / (2*a);
            roots[0] = 2.0;
        }
    }
    return true;
}

float* roots_quadratic_equation(double a, double b, double c) {
    //the first element is the number of the real roots, and other elements are the real roots.
    float *roots = new float[3];
    if (fabs(a) < DBL_EPSILON) {
        if (fabs(b) < DBL_EPSILON) {
            roots[0] = 0.0;
        }
        else {
            roots[1] = -c/b;
            roots[0] = 1.0;
        }
    }
    else {
        double d = b*b - 4*a*c;
        if (d < -DBL_EPSILON) {
            roots[0] = 0.0;
        }
        else if (fabs(d) <= DBL_EPSILON) {
            roots[1] = float((-b + 0) / (2*a));
            roots[2] = float((-b - 0) / (2*a));
            roots[0] = float(1.0);
        }
        else {
            roots[1] = float((-b + sqrt(d)) / (2*a));
            roots[2] = float((-b - sqrt(d)) / (2*a));
            roots[0] = float(2.0);
        }
    }
    return roots;
}
/*
float* roots_cubic_equation2(float a, float b, float c, float d) {
    //the first element is the number of the real roots, and other elements are the real roots.
    //Shengjin's formula
    float *roots = new float[4];
    if (a == 0) {
        roots = roots_quadratic_equation(b, c, d);
    }
    else {
        float A = b*b - 3*a*c;
        float B = b*c - 9*a*d;
        float C = c*c - 3*b*d;
        float deita = B*B - 4*A*C;

        if ((A == B) && (A == 0)) {
            //the three roots are the same
            if (a != 0) {
                roots[1] = -b/(3*a);
            }
            else {
                if (b != 0) {
                    roots[1] = -c/b;
                }
                else {
                    if (c != 0) {
                        roots[1] = -3*d/c;
                    }
                }
            }
            roots[2] = roots[1];
            roots[3] = roots[1];
            roots[0] = 3;
        }
        else if (deita > 0) {
            //only one real root
            float y1 = A*b + (3*a/2)*(-B + sqrt(deita));
            float y2 = A*b + (3*a/2)*(-B - sqrt(deita));
            float pow_y1, pow_y2;
            if (y1 < 0) {
            //for pow(a,b), when b is not int, a should not be negative.
                pow_y1 = - pow(-y1, 1.0/3.0);
            }
            else {
                pow_y1 = pow(y1, 1.0/3.0);
            }
            if (y2 < 0) {
                pow_y2 = - pow(-y2, 1.0/3.0);
            }
            else {
                pow_y2 = pow(y2, 1.0/3.0);
            }
            roots[1] = (-b - pow_y1 - pow_y2) / (3*a);
            roots[0] = 1;
        }
        else if (deita == 0) {
            //three real roots and two of them are the same
            float K = B/A;
            roots[1] = -b/a + K;
            roots[2] = -K/2;
            roots[3] = -K/2;
            roots[0] = 3;
        }
        else if (deita < 0) {
            //three different real roots
            float theta = acos((2*A*b-3*a*B) / (2*pow(A, 1.5)));
            roots[1] = (-b - 2*sqrt(A)*cos(theta/3)) / (3*a);
            roots[2] = (-b + sqrt(A) * (cos(theta/3) + sqrt(3)*sin(theta/3))) / (3*a);
            roots[3] = (-b + sqrt(A) * (cos(theta/3) - sqrt(3)*sin(theta/3))) / (3*a);
            roots[0] = 3;
        }
    }
    return roots;
}
*/
float* roots_cubic_equation(double a, double b, double c, double d) {
    //the first element is the number of the real roots, and other elements are the real roots.
    //Shengjin's formula
    float *roots = new float[4];
    if (fabs(a) < DBL_EPSILON) {
        float *roots2;
        roots2 = roots_quadratic_equation(b, c, d);
        for (int i=0; i<int(roots2[0])+1; i++) {
            roots[i] = roots2[i];
        }
        delete [] roots2;
    }
    else {
        double A = double(b*b - 3*a*c);
        double B = double(b*c - 9*a*d);
        double C = double(c*c - 3*b*d);
        double deita = double(B*B - 4*A*C);

        if ((A == B) && (fabs(A) < DBL_EPSILON)) {
            //the three roots are the same
            if (!(fabs(a) < DBL_EPSILON)) {
                roots[1] = -b/(3*a);
            }
            else {
                if (!(fabs(b) < DBL_EPSILON)) {
                    roots[1] = -c/b;
                }
                else {
                    if (!(fabs(c) < DBL_EPSILON)) {
                        roots[1] = -3*d/c;
                    }
                }
            }
            roots[2] = roots[1];
            roots[3] = roots[1];
            roots[0] = 3;
        }
        else if (deita > DBL_EPSILON) {
            //only one real root
            double y1 = double(A*b + (3*a)*(-B + sqrt(deita))/2);
            double y2 = double(A*b + (3*a)*(-B - sqrt(deita))/2);
            double pow_y1, pow_y2;
            if (y1 < -DBL_EPSILON) {
            //for pow(a,b), when b is not int, a should not be negative.
                pow_y1 = - pow(double(-y1), double(1.0/3.0));
            }
            else if (fabs(y1) <= DBL_EPSILON) {
                pow_y1 = 0;
            }
            else {
                pow_y1 = pow(double(y1), double(1.0/3.0));
            }
            if (y2 < -DBL_EPSILON) {
                pow_y2 = - pow(double(-y2), double(1.0/3.0));
            }
            else if (fabs(y2) <= DBL_EPSILON) {
                pow_y2 = 0;
            }
            else {
                pow_y2 = pow(double(y2), double(1.0/3.0));
            }
            roots[1] = float((-b - pow_y1 - pow_y2) / (3*a));
            roots[0] = float(1);
        }
        else if (fabs(deita) <= DBL_EPSILON) {
            //three real roots and two of them are the same
            double K = B/A;
            roots[1] = float(-b/a + K);
            roots[2] = float(-K/2);
            roots[3] = float(-K/2);
            roots[0] = float(3);
        }
        else if (deita < -DBL_EPSILON) {
            //three different real roots
            double theta = acos((2*A*b-3*a*B) / (2*pow(A, 1.5)));
            roots[1] = float((-b - 2*sqrt(A)*cos(theta/3)) / (3*a));
            roots[2] = float((-b + sqrt(A) * (cos(theta/3) + sqrt(3)*sin(theta/3))) / (3*a));
            roots[3] = float((-b + sqrt(A) * (cos(theta/3) - sqrt(3)*sin(theta/3))) / (3*a));
            roots[0] = float(3);
        }
    }
    return roots;
}

/*
float* roots_quartic_equation(float a, float b, float c, float d, float e) {
    //the first element is the number of the real roots, and other elements are the real roots.
    //Ferrari's Method.
    float *roots = new float[5];
    if (a == 0) {
        float *roots3;
        roots3 = roots_cubic_equation(b, c, d, e);
        for (int i=0; i<int(roots3[0])+1; i++) {
            roots[i] = roots3[i];
        }
        delete [] roots3;
    }
    else {
        float b1 = b/a;
        float c1 = c/a;
        float d1 = d/a;
        float e1 = e/a;
        if ((b1 == 0) && (c1 == 0) && (d1 == 0)) {
        //in this special case, such as a=1, b=c=d=0, e=-1, the roots should be +1 and -1
            if (e1 > 0) {
                roots[0] = 0.0;
            }
            else {
                roots[1] = sqrt(sqrt(-e1));
                roots[2] = -sqrt(sqrt(-e1));
                roots[0] = 2.0;
            }
        }
        else {
            float *roots_y;
            roots_y = roots_cubic_equation(-1.0, c1, 4*e1-b1*d1, d1*d1+e1*b1*b1-4*e1*c1);

            for (int i=0; i<(roots_y[0]+1); i++) {
                std::cout << "roots_y[" << i << "]=" << roots_y[i] << endl;
            }

            float y = roots_y[3];

            std::cout << "y = roots_y[3]:" << y << endl;

            delete [] roots_y;
            float A1, A2, B1, B2, C1, C2;
            if (b1*b1-4*c1+4*y == 0) {
                A1 = 1.0;
                A2 = 1.0;
                B1 = b1/2;
                B2 = b1/2;
                C1 = y/2;
                C2 = y/2;
            }
            else if (b1*b1-4*c1+4*y < 0) {
                A1 = 0.0;
                A2 = 1.0;
                B1 = 1.0;
                B2 = b1/2;
                C1 = (b1*y-2*d1)/(b1*b1-4*c1+4*y);
                C2 = y/2;
            }
            else {
                A1 = 1.0;
                A2 = 1.0;
                B1 = b1/2 - sqrt(b1*b1-4*c1+4*y)/2;
                B2 = b1/2 + sqrt(b1*b1-4*c1+4*y)/2;
                C1 = y/2 - (b1*y-2*d1)/(2*sqrt(b1*b1-4*c1+4*y));
                C2 = y/2 + (b1*y-2*d1)/(2*sqrt(b1*b1-4*c1+4*y));
            }
            float *roots_x1;
            float *roots_x2;
            roots_x1 = roots_quadratic_equation(A1, B1, C1);
            roots_x2 = roots_quadratic_equation(A2, B2, C2);
            if (roots_x1[0] != 0) {
                for (int i=1; i<roots_x1[0]+1; i++) {
                    roots[i] = roots_x1[i];
                }
            }
            if (roots_x2[0] != 0) {
                int roots_x1_number = int(roots_x1[0]);
                for (int j=1; j<roots_x2[0]+1; j++) {
                    roots[roots_x1_number+j] =roots_x2[j];
                }
            }
            roots[0] = roots_x1[0] + roots_x2[0];
            delete [] roots_x1;
            delete [] roots_x2;
        }
    }
    return roots;
}
*/

bool roots_quartic_equation2(double a, double b, double c, double d, double e, float (&roots)[5]) {
    //the first element is the number of the real roots, and other elements are the real roots.
    //Descartes's Method.
    if (fabs(a) < DBL_EPSILON) {
        float *roots3;
        roots3 = roots_cubic_equation(b, c, d, e);
        for (int i=0; i<int(roots3[0])+1; i++) {
            roots[i] = roots3[i];
        }
        delete [] roots3;
    }
    else {
        double a1 = b/a;
        double b1 = c/a;
        double c1 = d/a;
        double d1 = e/a;
        double p = (-3*a1*a1)/8 + b1;
        double q = (a1*a1*a1)/8 - (a1*b1)/2 + c1;
        double r = (-3*a1*a1*a1*a1)/256 + (a1*a1*b1)/16 - (a1*c1)/4 + d1;
        double sa = 2*p;
        double sb = p*p - 4*r;
        double sc = -q*q;
        double k, m, t;
        float *roots_s = roots_cubic_equation(1.0, sa, sb, sc);
#if VEC3_LOG == 1
        for (int i=0; i<(roots_s[0]+1); i++) {
            std::cout << "roots_s3[" << i << "]=" << roots_s[i] << endl;
        }
#endif // VEC3_LOG
        double temp;
        for (int i=1; i<int(roots_s[0]); i++) {
            for (int j=i+1; j<int(roots_s[0])+1; j++) {
                if (roots_s[i] > roots_s[j]) {
                    temp = roots_s[i];
                    roots_s[i] = roots_s[j];
                    roots_s[j] = temp;
                }
            }
        }
        if (roots_s[int(roots_s[0])] > DBL_EPSILON) {
            k = sqrt(double(roots_s[int(roots_s[0])]));
            delete [] roots_s;
        }
        else {
            if (fabs(q) < DBL_EPSILON) {
                k = 0.0;
                delete [] roots_s;
            }
            else {
                roots[0] = 0.0;
                delete [] roots_s;
                return true;
            }
        }
        if (fabs(k) == 0.0) {
            if (fabs(q) < DBL_EPSILON) {
                float *roots2;
                roots2 = roots_quadratic_equation(1.0, -p, r);
                if (roots2[0] == 0.0) {
                    roots[0] = 0.0;
                    delete [] roots2;
                    return true;
                }
                else {
                    m = roots2[1];
                    t = roots2[2];
                }
                delete [] roots2;
            }
            else {
                roots[0] = 0.0;
                return true;
            }
        }
        else {
            m = (k*k*k + k*p + q) / (2*k);
            t = (k*k*k + k*p - q) / (2*k);
        }

        float *roots_y1;
        float *roots_y2;
        roots_y1 = roots_quadratic_equation(1.0, k, t);
#if VEC3_LOG == 1
        for (int i=0; i<(roots_y1[0]+1); i++) {
            std::cout << "roots_y1[" << i << "]=" << roots_y1[i] << endl;
        }
#endif // VEC3_LOG
        roots_y2 = roots_quadratic_equation(1.0, -k, m);
#if VEC3_LOG == 1
        for (int i=0; i<(roots_y2[0]+1); i++) {
            std::cout << "roots_y2[" << i << "]=" << roots_y2[i] << endl;
        }
#endif // VEC3_LOG
        if (roots_y1[0] != 0.0) {
            for (int i=1; i<roots_y1[0]+1; i++) {
                roots[i] = roots_y1[i] - a1/4;
            }
        }
        if (roots_y2[0] != 0.0) {
            int roots_y1_number = int(roots_y1[0]);
            for (int j=1; j<roots_y2[0]+1; j++) {
                roots[roots_y1_number+j] =roots_y2[j] - a1/4;
            }
        }
        roots[0] = roots_y1[0] + roots_y2[0];
        delete [] roots_y1;
        delete [] roots_y2;
    }
    return true;
}

double* roots_quadratic_equation_rain(double a, double b, double c) {
    //the first element is the number of the real roots, and other elements are the real roots.
    double *roots = new double[3];
    if (fabs(a) < DBL_EPSILON) {
        if (fabs(b) < DBL_EPSILON) {
            roots[0] = 0.0;
        }
        else {
            roots[1] = -c/b;
            roots[0] = 1.0;
        }
    }
    else {
        double d = b*b - 4*a*c;
        if (d < -DBL_EPSILON) {
            roots[0] = 0.0;
        }
        else if (fabs(d) <= DBL_EPSILON) {
            roots[1] = double((-b) / (2*a));
            roots[2] = double((-b) / (2*a));
            roots[0] = double(1.0);
        }
        else {
            roots[1] = double((-b + sqrt(d)) / (2*a));
            roots[2] = double((-b - sqrt(d)) / (2*a));
            roots[0] = double(2.0);
        }
    }
    return roots;
}


double* roots_cubic_equation_rain(double a, double b, double c, double d) {
    //the first element is the number of the real roots, and other elements are the real roots.
    //Shengjin's formula
    double *roots = new double[4];
    if (fabs(a) < DBL_EPSILON) {
        double *roots2;
        roots2 = roots_quadratic_equation_rain(b, c, d);
        for (int i=0; i<int(roots2[0])+1; i++) {
            roots[i] = roots2[i];
        }
        delete [] roots2;
    }
    else {
        double A = double(b*b - 3.0*a*c);
        double B = double(b*c - 9.0*a*d);
        double C = double(c*c - 3.0*b*d);
        double deita = double(B*B - 4.0*A*C);

        if ((A == B) && (fabs(A) < DBL_EPSILON)) {
            //the three roots are the same
            if (!(fabs(a) < DBL_EPSILON)) {
                roots[1] = -b/(3.0*a);
            }
            else {
                if (!(fabs(b) < DBL_EPSILON)) {
                    roots[1] = -c/b;
                }
                else {
                    if (!(fabs(c) < DBL_EPSILON)) {
                        roots[1] = -3.0*d/c;
                    }
                }
            }
            roots[2] = roots[1];
            roots[3] = roots[1];
            roots[0] = double(3.0);
        }
        else if (deita > DBL_EPSILON) {
            //only one real root
            double y1 = double(A*b + (3*a)*(-B + sqrt(deita))/2);
            double y2 = double(A*b + (3*a)*(-B - sqrt(deita))/2);
            double pow_y1, pow_y2;
            if (y1 < -DBL_EPSILON) {
            //for pow(a,b), when b is not int, a should not be negative.
                pow_y1 = - pow(double(-y1), double(1.0/3.0));
            }
            else if (fabs(y1) <= DBL_EPSILON) {
                pow_y1 = 0;
            }
            else {
                pow_y1 = pow(double(y1), double(1.0/3.0));
            }
            if (y2 < -DBL_EPSILON) {
                pow_y2 = - pow(double(-y2), double(1.0/3.0));
            }
            else if (fabs(y2) <= DBL_EPSILON) {
                pow_y2 = 0;
            }
            else {
                pow_y2 = pow(double(y2), double(1.0/3.0));
            }
            roots[1] = double((-b - pow_y1 - pow_y2) / (3.0*a));
            roots[0] = double(1.0);
        }
        else if (fabs(deita) <= DBL_EPSILON) {
            //three real roots and two of them are the same
            double K = B/A;
            roots[1] = double(-b/a + K);
            roots[2] = double(-K/2.0);
            roots[3] = double(-K/2.0);
            roots[0] = double(3.0);
        }
        else if (deita < -DBL_EPSILON) {
            //three different real roots
            double theta = acos((2.0*A*b-3.0*a*B) / (2.0*pow(A, 1.5)));
            roots[1] = double((-b - 2.0*sqrt(A)*cos(theta/3.0)) / (3.0*a));
            roots[2] = double((-b + sqrt(A) * (cos(theta/3.0) + sqrt(double(3.0))*sin(theta/3.0))) / (3.0*a));
            roots[3] = double((-b + sqrt(A) * (cos(theta/3.0) - sqrt(double(3.0))*sin(theta/3.0))) / (3.0*a));
            roots[0] = double(3.0);
        }
    }
    return roots;
}


bool roots_quartic_equation2_rain(double a, double b, double c, double d, double e, double (&roots)[5]) {
    //the first element is the number of the real roots, and other elements are the real roots.
    //Descartes's Method.
    if (fabs(a) < DBL_EPSILON) {
        double *roots3;
        roots3 = roots_cubic_equation_rain(b, c, d, e);
        for (int i=0; i<int(roots3[0])+1; i++) {
            roots[i] = roots3[i];
        }
        delete [] roots3;
    }
    else {
        double a1 = b/a;
        double b1 = c/a;
        double c1 = d/a;
        double d1 = e/a;
        double p = (-3*a1*a1)/8 + b1;
        double q = (a1*a1*a1)/8 - (a1*b1)/2 + c1;
        double r = (-3*a1*a1*a1*a1)/256 + (a1*a1*b1)/16 - (a1*c1)/4 + d1;
        double sa = 2*p;
        double sb = p*p - 4*r;
        double sc = -q*q;
        double k, m, t;
        double *roots_s = roots_cubic_equation_rain(1.0, sa, sb, sc);
#if VEC3_LOG == 1
        for (int i=0; i<(roots_s[0]+1); i++) {
            std::cout << "roots_s3[" << i << "]=" << roots_s[i] << endl;
        }
//        std::cout << "DBL_EPSILON:" << DBL_EPSILON << endl;
#endif // VEC3_LOG
        double temp;
        for (int i=1; i<int(roots_s[0]); i++) {
            for (int j=i+1; j<int(roots_s[0])+1; j++) {
                if (roots_s[i] > roots_s[j]) {
                    temp = roots_s[i];
                    roots_s[i] = roots_s[j];
                    roots_s[j] = temp;
                }
            }
        }
        if (roots_s[int(roots_s[0])] > DBL_EPSILON) {
            k = sqrt(double(roots_s[int(roots_s[0])]));
            delete [] roots_s;
        }
        else {
            if (fabs(q) < DBL_EPSILON) {
                k = 0.0;
                delete [] roots_s;
            }
            else {
                roots[0] = 0.0;
                delete [] roots_s;
                return true;
            }
        }
        if (fabs(k) == 0.0) {
            if (fabs(q) < DBL_EPSILON) {
                double *roots2;
                roots2 = roots_quadratic_equation_rain(1.0, -p, r);
                if (roots2[0] == 0.0) {
                    roots[0] = 0.0;
                    delete [] roots2;
                    return true;
                }
                else {
                    m = roots2[1];
                    t = roots2[2];
                }
                delete [] roots2;
            }
            else {
                roots[0] = 0.0;
                return true;
            }
        }
        else {
            m = (k*k*k + k*p + q) / (2*k);
            t = (k*k*k + k*p - q) / (2*k);
        }

        double *roots_y1;
        double *roots_y2;
        roots_y1 = roots_quadratic_equation_rain(1.0, k, t);
#if VEC3_LOG == 1
        for (int i=0; i<(roots_y1[0]+1); i++) {
            std::cout << "roots_y1[" << i << "]=" << roots_y1[i] << endl;
        }
#endif // VEC3_LOG
        roots_y2 = roots_quadratic_equation_rain(1.0, -k, m);
#if VEC3_LOG == 1
        for (int i=0; i<(roots_y2[0]+1); i++) {
            std::cout << "roots_y2[" << i << "]=" << roots_y2[i] << endl;
        }
#endif // VEC3_LOG
        if (roots_y1[0] != 0.0) {
            for (int i=1; i<roots_y1[0]+1; i++) {
                roots[i] = roots_y1[i] - a1/4;
            }
        }
        if (roots_y2[0] != 0.0) {
            int roots_y1_number = int(roots_y1[0]);
            for (int j=1; j<roots_y2[0]+1; j++) {
                roots[roots_y1_number+j] =roots_y2[j] - a1/4;
            }
        }
        roots[0] = roots_y1[0] + roots_y2[0];
        delete [] roots_y1;
        delete [] roots_y2;
    }
    return true;
}


bool get_matrix_3_1(const vec3 a, float (&m)[3][1]) {
        m[0][0] = a.x();
        m[1][0] = a.y();
        m[2][0] = a.z();
        return true;
}

bool get_matrix_3_3(const vec3 a, const vec3 b, const vec3 c, float (&m)[3][3]) {
        m[0][0] = a.x();
        m[1][0] = a.y();
        m[2][0] = a.z();
        m[0][1] = b.x();
        m[1][1] = b.y();
        m[2][1] = b.z();
        m[0][2] = c.x();
        m[1][2] = c.y();
        m[2][2] = c.z();
        return true;
}

bool get_matrix_inverse_3_3(const float m[3][3], float (&inverse)[3][3]) {
        float det_m = m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1] -
                      m[0][2]*m[1][1]*m[2][0] - m[0][0]*m[1][2]*m[2][1] - m[0][1]*m[1][0]*m[2][2];
        if (fabs(det_m) < 1e-6) {
            return false;
        }
        else {
            vec3 a = (1/det_m)*vec3(m[1][1]*m[2][2] - m[1][2]*m[2][1],
                                    m[1][2]*m[2][0] - m[1][0]*m[2][2],
                                    m[1][0]*m[2][1] - m[1][1]*m[2][0]);
            vec3 b = (1/det_m)*vec3(m[0][2]*m[2][1] - m[0][1]*m[2][2],
                                    m[0][0]*m[2][2] - m[0][2]*m[2][0],
                                    m[0][1]*m[2][0] - m[0][0]*m[2][1]);
            vec3 c = (1/det_m)*vec3(m[0][1]*m[1][2] - m[0][2]*m[1][1],
                                    m[0][2]*m[1][0] - m[0][0]*m[1][2],
                                    m[0][0]*m[1][1] - m[0][1]*m[1][0]);
            get_matrix_3_3(a, b, c, inverse);
            return true;
        }
}

bool matrix_3_3_multiply_3_1(const float m1[3][3], const float m2[3][1], float (&m)[3][1]) {
        m[0][0] = m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0] + m1[0][2]*m2[2][0];
        m[1][0] = m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0] + m1[1][2]*m2[2][0];
        m[2][0] = m1[2][0]*m2[0][0] + m1[2][1]*m2[1][0] + m1[2][2]*m2[2][0];
        return true;
}

bool matrix_3_1_minus_3_1(const float m1[3][1], const float m2[3][1], float (&m)[3][1]) {
        m[0][0] = m1[0][0] - m2[0][0];
        m[1][0] = m1[1][0] - m2[1][0];
        m[2][0] = m1[2][0] - m2[2][0];
        return true;
}

bool get_bezier_matrix_xyz_4_4(const vec3 point[4][4], float (&bezier_x)[4][4], float (&bezier_y)[4][4], float (&bezier_z)[4][4]) {
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                bezier_x[i][j] = point[i][j].x();
                bezier_y[i][j] = point[i][j].y();
                bezier_z[i][j] = point[i][j].z();
            }
        }
        return true;
}

bool matrix_4_4_multiply_4_4(const float matrix1[4][4], const float matrix2[4][4], float (&result)[4][4]) {
        for (int k=0; k<4; k++) {
            for (int i=0; i<4; i++) {
                result[i][k] = 0.0;
                for (int j=0; j<4; j++) {
                    result[i][k] = result[i][k] + matrix1[i][j]*matrix2[j][k];
                }
            }
        }
        return true;
}
bool matrix_4_4_multiply_4_1(const float matrix1[4][4], const float matrix2[4][1], float (&result)[4][1]) {
        for (int k=0; k<1; k++) {
            for (int i=0; i<4; i++) {
                result[i][k] = 0.0;
                for (int j=0; j<4; j++) {
                    result[i][k] = result[i][k] + matrix1[i][j]*matrix2[j][k];
                }
            }
        }
        return true;
}

bool get_matrix_transpose_4_4(const float matrix[4][4], float (&result)[4][4]) {
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            result[j][i] = matrix[i][j];
        }
    }
    return true;
}

bool get_teapot_data(int (&patches)[32][16], float (&vertices)[306][3]) {
        ifstream infile( ".\\teaset\\teapot");
        char str[100];
        char *token;
        int item_num;
        int patch_num = 0;
        int vertex_num = 0;
        int flag = 0;
        int tokens_i[16];
        float tokens_f[3];
        while (infile >> str) {
            item_num = 0;
            token = strtok(str, ",");
            if ((flag == 0) || (flag == 1)) {
                sscanf(token, "%d", &(tokens_i[item_num]));
                if (tokens_i[0] == 32) {
                    flag = 1;
                }
                else if (tokens_i[0] == 306) {
                    flag = 2;
                }
                else {
                    patches[patch_num][item_num] = tokens_i[item_num];
                    item_num ++;
                }
                token = strtok(NULL, ",");
            }

            while (token != NULL) {
                if (flag == 1) {
                    sscanf(token, "%d", &(tokens_i[item_num]));
                    patches[patch_num][item_num] = tokens_i[item_num];
                }
                if (flag == 2) {
                    sscanf(token, "%f", &(tokens_f[item_num]));
                    vertices[vertex_num][item_num] = tokens_f[item_num];
                }
                item_num ++;
                token = strtok(NULL, ",");
            }
            if ((flag == 1) && (tokens_i[0] != 32)) {patch_num ++;}
            if (flag == 2) {
                if ((tokens_i[0] == 306)) {tokens_i[0] = 0;}
                else {vertex_num ++;}
            }
        }
        infile.close();
        return true;
}

bool roots_num_equation_6th(float ee6[7], double a, double b, int &num) {
        double temp1, temp2, aaa, bbb;
        double ee77[7][7] = {0.0};
        double fun_a[7] = {0.0};
        double fun_b[7] = {0.0};
        int change_num_a = 0;
        int change_num_b = 0;
        if (fabs(a)<1e-6) {
            aaa = 1e-6;
        }
        else {
            aaa = a;
        }
        if (fabs(b)<1e-6) {
            bbb = 1e-6;
        }
        else {
            bbb = b;
        }

        for (int i=0; i<6; i++) {
        //determine f0, f1
            ee77[0][i] = ee6[i];
            ee77[1][i] = (6-i)*ee6[i];
        }
        ee77[0][6] = ee6[6];

        for (int i=2; i<6; i++) {
        //determine f2, f3, f4, f5
            temp1 = ee77[i-2][0] / ee77[i-1][0];
            temp2 = (ee77[i-2][1] - temp1*ee77[i-1][1]) / ee77[i-1][0];
            for (int j=0; j<(6-i); j++) {
                ee77[i][j] = temp1*ee77[i-1][j+2] + temp2*ee77[i-1][j+1] - ee77[i-2][j+2];
                if (fabs(ee77[i][j])<1e-5) {
                    ee77[i][j] = 0.0;
                }
            }
            ee77[i][6-i] = temp2*ee77[i-1][6-i+1] - ee77[i-2][6-i+2];
            if (fabs(ee77[i][6-i])<1e-5) {
                ee77[i][6-i] = 0.0;
            }

            if ((ee77[i][0] == 0.0) && (ee77[i][1] == 0.0) && (ee77[i][2] == 0.0) &&
                (ee77[i][3] == 0.0) && (ee77[i][4] == 0.0) && (ee77[i][5] == 0.0) &&
                (ee77[i][6] == 0.0)) {
                break;
            }
        }
        if (fabs(ee77[5][0])>1e-5) {
        //determine f6
            temp1 = ee77[4][0] / ee77[5][0];
            temp2 = (ee77[4][1] - temp1*ee77[5][1]) / ee77[5][0];
            ee77[6][0] = temp2*ee77[5][1] - ee77[4][2];
            if (fabs(ee77[6][0])<1e-5) {
                ee77[6][0] = 0.0;
            }
        }
        else {
            ee77[6][0] = 0.0;
        }

        for (int j=0; j<7; j++) {
            fun_a[0] = fun_a[0] + ee77[0][j]*pow(aaa, (6-j));
            fun_b[0] = fun_b[0] + ee77[0][j]*pow(bbb, (6-j));
        }
        for (int i=1; i<7; i++) {
            for (int j=0; j<(7-i); j++) {
                fun_a[i] = fun_a[i] + ee77[i][j]*pow(aaa, ((6-i)-j));
                fun_b[i] = fun_b[i] + ee77[i][j]*pow(bbb, ((6-i)-j));
            }
            if ((fun_a[i]*fun_a[i-1]) < 0) {
                change_num_a ++;
            }
            if ((fun_b[i]*fun_b[i-1]) < 0) {
                change_num_b ++;
            }
        }
        num = abs(change_num_a - change_num_b);
        return true;
}

bool get_equation_6th_function_and_derivative(float ee6[7], double xxx, double& fnc, double& drv) {
    /*determine function value and derivative value at xxx*/
        fnc = 0.0;
        drv = 0.0;
        for (int i=0; i<7; i++) {
            fnc = fnc + ee6[i]*pow(xxx, (6-i));
            if (i<6) {
                drv = drv + (6-i)*ee6[i]*pow(xxx, (5-i));
            }
        }
        return true;
}

bool get_root_by_newton_iteration_for_equation_6th(float ee6[7], double x0[2], float tol, double (&root)[2]) {
    /*find the single root in the interval [x0[0], x0[1]] by newton iteration */
        double t_k, t_k1, ft_k, ft_d_k;
        int get = 0;

        for (int i=0; i<2; i++) {
            t_k = x0[i];
            for (int k=0; k<20; k++) {
                if (!(isnan(t_k))) {
                    get_equation_6th_function_and_derivative(ee6, t_k, ft_k, ft_d_k);
                    if ((ft_d_k != 0) && !(isnan(ft_k)) && !(isnan(ft_d_k))) {
                        t_k1 = t_k - ft_k/ft_d_k;
                        if (((fabs(t_k1)>=1e-2)&&(fabs((t_k1 - t_k)/t_k1) < tol)) ||
                            ((fabs(t_k1)<1e-2)&&(fabs((t_k1 - t_k)) < tol))) {
                            if ((t_k1 > x0[0]) && (t_k1<=x0[1])) {
                                root[1] = (fabs(t_k1)<(tol*10))? 0.0:t_k1;
                                root[0] = 1.0;
                                get ++;
                            }
                            else {
                                x0[1] = (x0[0] + x0[1])/2;
                            }
                            break;
                        }
                        else {
                            t_k = t_k1;
                        }
                    }
                    else {
                        break;
                    }
                }
                else {
                    break;
                }
            }
            if (get == 1) {
                break;
            }
        }
        root[0] = double(get);
        return true;
}

bool roots_equation_6th(float ee6[7], float a, float b, float tol, float (&roots)[7]) {
        double root[2], xxx0[2], temp;
        double left = a;
        double right = b;
        int total_num, interval_num;
        int offset = 0;
        roots[0] = 0.0;
        roots_num_equation_6th(ee6, double(a), double(b), total_num);
        if (total_num == 0) {
            return true;
        }
        else {
            interval_num = total_num;
            for (int i=1; i<(total_num+1+offset); i++) {
                while(interval_num != 1) {
                    if (left > right) {
                        temp = left;
                        left = right;
                        right = temp;
                    }
                    right = (left + right)/2;

                    roots_num_equation_6th(ee6, left, right, interval_num);
                    if (interval_num == 0) {
                        left = right;
                        right = b;
                        if (left == right) {
                        /*s2: handle the situation that both left and right value reach the right end of the whole interval*/
                            break;
                        }
                    }
                    else if (fabs(left-right)<tol) {
                    /*s1: left and right value are are so close that their values even don't change under sub-dividing*/
                        break;
                    }
                }
                if (interval_num == 0) {
                /*s2: handle the situation that both left and right value reach the right end of the whole interval*/
                    break;
                }
                xxx0[0] = left;
                xxx0[1] = right;
                get_root_by_newton_iteration_for_equation_6th(ee6, xxx0, tol, root);
                if (root[0] != 0.0) {
                    roots[0] = roots[0] + 1.0;
                    roots[int(roots[0])] = root[1];
                    left = right;
                    right = b;
                }
                else if ((root[0] == 0.0)&&(fabs(left-right)<tol)&&(interval_num>=1)) {
                /*s1: left and right value are are so close that their values even don't change under sub-dividing,
                the interval_num may be more than one, so we modify "interval_num==1" to "interval_num>=1"*/
                    roots[0] = roots[0] + 1.0;
                    roots[int(roots[0])] = left;
                    left = right;
                    right = b;
                }
                else {
                    offset ++;
                    right = (left + right)/2;
                }

                roots_num_equation_6th(ee6, left, right, interval_num);
            }
        }
        return true;
}

bool roots_num_equation_10th(float ee10[11], double a, double b, int &num) {
        double temp1, temp2, aaa, bbb;
        double ee1111[11][11] = {0.0};
        double fun_a[11] = {0.0};
        double fun_b[11] = {0.0};
        int change_num_a = 0;
        int change_num_b = 0;
        if (fabs(a)<1e-6) {
            aaa = 1e-6;
        }
        else {
            aaa = a;
        }
        if (fabs(b)<1e-6) {
            bbb = 1e-6;
        }
        else {
            bbb = b;
        }

        for (int i=0; i<10; i++) {
        //determine f0, f1
            ee1111[0][i] = ee10[i];
            ee1111[1][i] = (10-i)*ee10[i];
        }
        ee1111[0][10] = ee10[10];

        for (int i=2; i<10; i++) {
        //determine f2, f3, f4, f5, f6, f7, f8, f9
            temp1 = ee1111[i-2][0] / ee1111[i-1][0];
            temp2 = (ee1111[i-2][1] - temp1*ee1111[i-1][1]) / ee1111[i-1][0];
            for (int j=0; j<(10-i); j++) {
                ee1111[i][j] = temp1*ee1111[i-1][j+2] + temp2*ee1111[i-1][j+1] - ee1111[i-2][j+2];
                if (fabs(ee1111[i][j])<1e-5) {
                    ee1111[i][j] = 0.0;
                }
            }
            ee1111[i][10-i] = temp2*ee1111[i-1][10-i+1] - ee1111[i-2][10-i+2];
            if (fabs(ee1111[i][10-i])<1e-5) {
                ee1111[i][10-i] = 0.0;
            }

            if ((ee1111[i][0] == 0.0) && (ee1111[i][1] == 0.0) && (ee1111[i][2] == 0.0) &&
                (ee1111[i][3] == 0.0) && (ee1111[i][4] == 0.0) && (ee1111[i][5] == 0.0) &&
                (ee1111[i][6] == 0.0) && (ee1111[i][7] == 0.0) && (ee1111[i][8] == 0.0) &&
                (ee1111[i][9] == 0.0) && (ee1111[i][10] == 0.0)) {
                break;
            }
        }
        if (fabs(ee1111[9][0])>1e-5) {
        //determine f10
            temp1 = ee1111[8][0] / ee1111[9][0];
            temp2 = (ee1111[8][1] - temp1*ee1111[9][1]) / ee1111[9][0];
            ee1111[10][0] = temp2*ee1111[9][1] - ee1111[8][2];
            if (fabs(ee1111[10][0])<1e-5) {
                ee1111[10][0] = 0.0;
            }
        }
        else {
            ee1111[10][0] = 0.0;
        }

        for (int j=0; j<11; j++) {
            fun_a[0] = fun_a[0] + ee1111[0][j]*pow(aaa, (6-j));
            fun_b[0] = fun_b[0] + ee1111[0][j]*pow(bbb, (6-j));
        }
        for (int i=1; i<11; i++) {
            for (int j=0; j<(11-i); j++) {
                fun_a[i] = fun_a[i] + ee1111[i][j]*pow(aaa, ((6-i)-j));
                fun_b[i] = fun_b[i] + ee1111[i][j]*pow(bbb, ((6-i)-j));
            }
            if ((fun_a[i]*fun_a[i-1]) < 0) {
                change_num_a ++;
            }
            if ((fun_b[i]*fun_b[i-1]) < 0) {
                change_num_b ++;
            }
        }
        num = abs(change_num_a - change_num_b);
        return true;
}

bool get_equation_10th_function_and_derivative(float ee10[11], double xxx, double& fnc, double& drv) {
    /*determine function value and derivative value at xxx*/
        fnc = 0.0;
        drv = 0.0;
        for (int i=0; i<11; i++) {
            fnc = fnc + ee10[i]*pow(xxx, (10-i));
            if (i<10) {
                drv = drv + (10-i)*ee10[i]*pow(xxx, (9-i));
            }
        }
        return true;
}

bool get_root_by_newton_iteration_for_equation_10th(float ee10[11], double x0[2], float tol, double (&root)[2]) {
    /*find the single root in the interval [x0[0], x0[1]] by newton iteration */
        double t_k, t_k1, ft_k, ft_d_k;
        int get = 0;

        for (int i=0; i<2; i++) {
            t_k = x0[i];
            for (int k=0; k<20; k++) {
                if (!(isnan(t_k))) {
                    get_equation_10th_function_and_derivative(ee10, t_k, ft_k, ft_d_k);
                    if ((ft_d_k != 0) && !(isnan(ft_k)) && !(isnan(ft_d_k))) {
                        t_k1 = t_k - ft_k/ft_d_k;
                        if (((fabs(t_k1)>=1e-2)&&(fabs((t_k1 - t_k)/t_k1) < tol)) ||
                            ((fabs(t_k1)<1e-2)&&(fabs((t_k1 - t_k)) < tol))) {
                            if ((t_k1 > x0[0]) && (t_k1<=x0[1])) {
                                root[1] = (fabs(t_k1)<(tol*10))? 0.0:t_k1;
                                root[0] = 1.0;
                                get ++;
                            }
                            else {
                                x0[1] = (x0[0] + x0[1])/2;
                            }
                            break;
                        }
                        else {
                            t_k = t_k1;
                        }
                    }
                    else {
                        break;
                    }
                }
                else {
                    break;
                }
            }
            if (get == 1) {
                break;
            }
        }
        root[0] = double(get);
        return true;
}

bool roots_equation_10th(float ee10[11], float a, float b, float tol, float (&roots)[11]) {
        double root[2], xxx0[2], temp;
        double left = a;
        double right = b;
        int total_num, interval_num;
        int offset = 0;
        roots[0] = 0.0;
        roots_num_equation_10th(ee10, double(a), double(b), total_num);
        if (total_num == 0) {
            return true;
        }
        else {
            interval_num = total_num;
            for (int i=1; i<(total_num+1+offset); i++) {
                while(interval_num != 1) {
                    if (left > right) {
                        temp = left;
                        left = right;
                        right = temp;
                    }
                    right = (left + right)/2;

                    roots_num_equation_10th(ee10, left, right, interval_num);
                    if (interval_num == 0) {
//                        left = right;
//                        right = b;
                        right = 2*right - left;
                        left = (left + right) / 2;
                        if (left == right) {
                        /*s2: handle the situation that both left and right value reach the right end of the whole interval*/
                            break;
                        }
                    }
                    else if (fabs(left-right)<tol) {
                    /*s1: left and right value are are so close that their values even don't change under sub-dividing*/
                        break;
                    }
                }
                if (interval_num == 0) {
                /*s2: handle the situation that both left and right value reach the right end of the whole interval*/
                    break;
                }
                xxx0[0] = left;
                xxx0[1] = right;
                get_root_by_newton_iteration_for_equation_10th(ee10, xxx0, tol, root);
                if (root[0] != 0.0) {
                    roots[0] = roots[0] + 1.0;
                    roots[int(roots[0])] = root[1];
//                        left = right;
//                        right = b;
                    right = 2*right - left;
                    left = (left + right) / 2;
                }
                else if ((root[0] == 0.0)&&(fabs(left-right)<tol)&&(interval_num>=1)) {
                /*s1: left and right value are are so close that their values even don't change under sub-dividing,
                the interval_num may be more than one, so we modify "interval_num==1" to "interval_num>=1"*/
                    roots[0] = roots[0] + 1.0;
                    roots[int(roots[0])] = left;
//                        left = right;
//                        right = b;
                    right = 2*right - left;
                    left = (left + right) / 2;
                }
                else {
                    offset ++;
                    right = (left + right)/2;
                }

                roots_num_equation_10th(ee10, left, right, interval_num);
            }
        }
        return true;
}

struct IntervalStruct
{
    double key;
    double leftValue;
    double rightValue;
    int rootNum;
};
typedef IntervalStruct ItemType;
struct IntervalTreeNode
{
    ItemType info;
    IntervalTreeNode* left;
    IntervalTreeNode* right;
};

void CreateIntervalTree(float ee10[11], double a, double b, int num, IntervalTreeNode *rootNode) {
    double left = a;
    double right = b;
    int leftNum, rightNum;
    rootNode->info.key = (left+right)/2;
    rootNode->info.leftValue = left;
    rootNode->info.rightValue = right;
    rootNode->info.rootNum = num;
    rootNode->left = NULL;
    rootNode->right = NULL;
    if ((num > 1) && (fabs(left-right)>1e-6)) {
        IntervalTreeNode *leftNode = new IntervalTreeNode;
        rootNode->left = leftNode;
        left = rootNode->info.leftValue;
        right = (rootNode->info.leftValue + rootNode->info.rightValue) / 2;
        roots_num_equation_10th(ee10, left, right, leftNum);
        CreateIntervalTree(ee10, left, right, leftNum, leftNode);

        IntervalTreeNode *rightNode = new IntervalTreeNode;
        rootNode->right = rightNode;
        left = (rootNode->info.leftValue + rootNode->info.rightValue) / 2;
        right = rootNode->info.rightValue;
//        roots_num_equation_10th(ee10, left, right, rightNum);
        rightNum = num-leftNum;
        CreateIntervalTree(ee10, left, right, rightNum, rightNode);
    }
}

void TraverseLeavesForValidIntervals(IntervalTreeNode *rootNode, int &counter, ItemType (&leaves)[11]) {
    if(rootNode!=NULL) {
        TraverseLeavesForValidIntervals(rootNode->left, counter, leaves);
        if ((rootNode->left == NULL) && (rootNode->right == NULL) && (rootNode->info.rootNum != 0)) {
            if (counter == 0) {counter = 1;}
            leaves[counter].leftValue = rootNode->info.leftValue;
            leaves[counter].rightValue = rootNode->info.rightValue;
            leaves[counter].rootNum = rootNode->info.rootNum;
            leaves[0].rootNum = counter;
            std::cout << "[" << leaves[counter].leftValue << ", " << leaves[counter].rightValue << "]: rootNum " << leaves[counter].rootNum << endl;
            counter ++;
        }
        TraverseLeavesForValidIntervals(rootNode->right, counter, leaves);
    }
}

/*as to single real root intervals, the interval ends may be not close enough to the root,
so we use this function to narrow the intervals*/
void NarrowValidIntervals(float ee10[11], double a, double b, double &na, double &nb) {
    int num = 0;
    na = a;
    nb = (a+b)/2;
    roots_num_equation_10th(ee10, na, nb, num);
    if (num == 0) {
        na = (a+b)/2;
        nb = b;
    }
}

bool roots_equation_10th_bt(float ee10[11], float a, float b, float tol, float (&roots)[11]) {
    IntervalTreeNode *rootNode = new IntervalTreeNode;
    ItemType intervals[11];
    double root[2], xxx0[2];
    double left = a;
    double right = b;
    double left_v, right_v;
    int totalNum, num;
    bool found;
    num = 0;
    roots[0] = 0.0;
    roots_num_equation_10th(ee10, left, right, totalNum);
    CreateIntervalTree(ee10, left, right, totalNum, rootNode);
    TraverseLeavesForValidIntervals(rootNode, num, intervals);
    for (int i=1; i<num; i++) {
        left = intervals[i].leftValue;
        right = intervals[i].rightValue;
        found = false;
        while (!found) {
            xxx0[0] = left;
            xxx0[1] = right;
            get_root_by_newton_iteration_for_equation_10th(ee10, xxx0, tol, root);
            if (root[0] != 0.0) {
                roots[0] = roots[0] + 1.0;
                roots[int(roots[0])] = root[1];
                found = true;
            }
            else if ((root[0] == 0.0)&&(fabs(left-right)<tol)) {
                roots[0] = roots[0] + 1.0;
                roots[int(roots[0])] = left;
                found = true;
            }
            else {
                NarrowValidIntervals(ee10, left, right, left_v, right_v);
                left = left_v;
                right = right_v;
            }
        }
    }
    return true;
}


