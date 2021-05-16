// Bismillahir Rahmanir Raheem
#include<iostream>
#define PI (2*acos(0.0))
#include "bitmap_image.hpp"

using namespace std;

///////// class ///////////

struct point {
    double x,y,z,w;

    point() {
        x = 0; y = 0; z = 0; w = 1;
    }

    point(double xx,double yy,double zz) {
        x = xx; y = yy; z = zz; w = 1;
    }
};



struct matrix {
    double val[4][4];

    matrix() {
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++) {
                val[i][j] = 0;
            }
        }

        val[3][3] = 1;
    }

    matrix(double x11,double x12,double x13,double x14,
           double x21,double x22,double x23,double x24,
           double x31,double x32,double x33,double x34,
           double x41,double x42,double x43,double x44) {

        val[0][0] = x11; val[0][1] = x12; val[0][2] = x13; val[0][3] = x14;
        val[1][0] = x21; val[1][1] = x22; val[1][2] = x23; val[1][3] = x24;
        val[2][0] = x31; val[2][1] = x32; val[2][2] = x33; val[2][3] = x34;
        val[3][0] = x41; val[3][1] = x42; val[3][2] = x43; val[3][3] = x44;

    }

};

//////// func /////////////

void print_matrix(matrix m,string s) {
    cout << " >> "+s << endl;
    for(int i=0;i<4;i++) {
        for(int j=0;j<4;j++) {
            cout << m.val[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}


void print_vect(point p,string s) {
    cout << s + " : " ;
    cout << p.x << " " << p.y << " " << p.z << endl;
}

struct point transformPoint(struct matrix T, struct point P) {

    point ret;

    ret.x = T.val[0][0]*P.x + T.val[0][1]*P.y + T.val[0][2]*P.z + T.val[0][3]*P.w;
    ret.y = T.val[1][0]*P.x + T.val[1][1]*P.y + T.val[1][2]*P.z + T.val[1][3]*P.w;
    ret.z = T.val[2][0]*P.x + T.val[2][1]*P.y + T.val[2][2]*P.z + T.val[2][3]*P.w;
    ret.w = T.val[3][0]*P.x + T.val[3][1]*P.y + T.val[3][2]*P.z + T.val[3][3]*P.w;

    return ret;
}


struct matrix product(struct matrix L,struct matrix R) {
    matrix ret;

    for(int i = 0; i < 4; ++i)
    {
        for(int j = 0; j < 4; ++j)
        {
            for(int k=0; k<4; ++k)
            {
                ret.val[i][j] += L.val[i][k] * R.val[k][j];
            }
        }
    }

    return ret;

}

struct point normalize(struct point p) {
    struct point ret;

    double div = sqrt(  p.x*p.x + p.y*p.y + p.z*p.z );

    ret.x = p.x / div;
    ret.y = p.y / div;
    ret.z = p.z / div;

    return ret;
}

struct point cross_product(struct point p1,struct point p2)
{
    struct point ret;

    ret.x = p1.y * p2.z - p1.z * p2.y;
    ret.y = p1.z * p2.x - p1.x * p2.z;
    ret.z = p1.x * p2.y - p1.y * p2.x;


    return ret;

}

double dot_product (struct point a,struct point b) {
    return (a.x*b.x + a.y*b.y + a.z*b.z);
}

struct point R(point L,point R,double angle) {

    angle = (PI * angle) / 180.0;
    double cosO = cos(angle);
    double sinO = sqrt( 1-cosO*cosO );

    double dot = dot_product(R,L);
    point cross = cross_product(R,L);

    struct point c;
    c.x = R.x*cosO + L.x*dot*(1-cosO)  + cross.x*sinO;
    c.y = R.y*cosO + L.y*dot*(1-cosO)  + cross.y*sinO;
    c.z = R.z*cosO + L.z*dot*(1-cosO)  + cross.z*sinO;


    return c;
}

struct matrix rotate(double angle,point a)
{


    a = normalize(a);

    point i (1,0,0);
    point j (0,1,0);
    point k (0,0,1);

    point c1 = R(i,a,angle);
    point c2 = R(j,a,angle);
    point c3 = R(k,a,angle);


    matrix ret(
            c1.x, c2.x, c3.x, 0,
            c1.y, c2.y, c3.y, 0,
            c1.z, c2.z, c3.z, 0,
            0,    0,    0,    1
    );

    return ret;
}

struct point add(struct point one,struct point two) {

    return point( one.x+two.x, one.y+two.y, one.z+two.z );
}

struct point getVector(struct point a,struct point b) {
    point r;
    r.x = b.x-a.x;
    r.y = b.y-a.y;
    r.z = b.z-a.z;
    return r;
}

struct point negVector(struct point v) {

    return point(-v.x,-v.y,-v.z);
}




///////// var /////////////

double eyeX, eyeY, eyeZ, lookX, lookY, lookZ, upX, upY, upZ, fovY, aspectRatio, near, f;
string command;

stack<matrix> S;
stack<int> L;

int main()
{

    ///////////// init ////////////

    matrix I4(
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1
    );

    S.push(I4);

    std::ios::sync_with_stdio(false);


    ////////////////////////////////////////// STAGE 1 ////////////////////////////////////////////////


    ifstream in1;
    in1.open("scene.txt");


    ofstream out1;
    out1.open("stage1.txt");

    out1 << fixed << setprecision(2);


    //////////// load data ////////////

    in1 >> eyeX >> eyeY >> eyeZ >> lookX >> lookY >> lookZ >> upX >> upY >> upZ >> fovY >> aspectRatio >> near >> f;



    while (true)
    {
        in1 >>  command;


        if(command.compare("triangle")==0) {

            point p1, p2, p3;

            in1>> p1.x >> p1.y >> p1.z;
            in1>> p2.x >> p2.y >> p2.z;
            in1>> p3.x >> p3.y >> p3.z;

            point pt1 = transformPoint(S.top(),p1);
            point pt2 = transformPoint(S.top(),p2);
            point pt3 = transformPoint(S.top(),p3);

            ////////////// print ////////////////////



            out1 << pt1.x << "\t" << pt1.y << "\t" << pt1.z << endl;
            out1 << pt2.x << "\t" << pt2.y << "\t" << pt2.z << endl;
            out1 << pt3.x << "\t" << pt3.y << "\t" << pt3.z << endl;
            out1 << endl;


        }

        else if(command.compare("translate")==0) {
            point T;
            in1 >> T.x >> T.y >> T.z;


            matrix Tr(
                    1, 0, 0, T.x,
                    0, 1, 0, T.y,
                    0, 0, 1, T.z,
                    0, 0, 0, 1
            );

            S.push( product(S.top(),Tr) );

        }

        else if(command.compare("scale")==0) {
            point s;
            in1 >> s.x >> s.y >> s.z;

            matrix Sc(
                    s.x, 0, 0, 0,
                    0, s.y, 0, 0,
                    0, 0, s.z, 0,
                    0, 0, 0, 1
            );

            S.push( product(S.top(),Sc) );


        }

        else if(command.compare("rotate")==0) {

            double angle;
            point axis;

            in1 >> angle >> axis.x >> axis.y >> axis.z;

            matrix Rt = rotate(angle,axis);

            S.push( product(S.top(),Rt) );

        }

        else if(command.compare("push")==0) {
            L.push( S.size() );
        }

        else if(command.compare("pop")==0) {

            int s = L.top();
            L.pop();

            while (S.size() > s )
            {
                S.pop();
            }

        }

        else if(command.compare("end")==0) {
            break;
        }

        else continue;



    }

    in1.close();
    out1.close();


    ////////////////////////////////////////// STAGE 2 ////////////////////////////////////////////////

    ifstream in2;
    in2.open("stage1.txt");


    ofstream out2;
    out2.open("stage2.txt");

    out2 << fixed << setprecision(2);



    point look(lookX,lookY,lookZ);
    point up(upX,upY,upZ);
    point eye(eyeX,eyeY,eyeZ);

    point l = add( look, negVector(eye) );
    l = normalize(l);

    point r = cross_product(l,up);
    r = normalize(r);

    point u = cross_product(r,l);

    matrix T(
            1, 0, 0, -eyeX,
            0, 1, 0, -eyeY,
            0, 0, 1, -eyeZ,
            0, 0, 0, 1
    );



    matrix R(
            r.x, r.y, r.z, 0,
            u.x, u.y, u.z, 0,
            -l.x, -l.y, -l.z, 0,
            0, 0, 0, 1
    );

    print_matrix(R,"R");
    print_matrix(T,"T");

    matrix V = product(R,T);

    print_matrix(V,"V");


    while (true)
    {
        matrix m;


        for(int i=0;i<3;i++) {
            for(int j=0;j<3;j++) {
                if(in2>>m.val[i][j]) {
                    //////
                }
                else
                    goto out;
            }

        }

        print_matrix(m,"input ");

        //// found a new matrix; ////////

        matrix mm = product(V,m);


        for(int i=0;i<3;i++) {
            for(int j=0;j<3;j++) {
                out2 << mm.val[i][j] << "\t";
            }
            out2 << endl;
        }
        out2 << endl;



    }

    out:





    return 0;
}