// Efaz, 1605034

#include <iostream>
#include <stack>
#include <cmath>
#include <fstream>
#include <vector>

#define PI (2*acos(0.0))
#include "bitmap_image.hpp"

using namespace std;

///////// class ///////////


class point {
    

    public:
    double x,y,z,w;

    point() {
        x = 0; y = 0; z = 0; w = 1;
    }

    point(double xx,double yy,double zz) {
        x = xx; y = yy; z = zz; w = 1;
    }

    void setVal(double xx,double yy,double zz) {
        x = xx; y = yy; z = zz; 
    }

     void normalize()
    {
        double r = sqrt((x * x) + (y * y) + (z * z));
        x = x / r;
        y = y / r;
        z = z / r;
    }

    point operator+(point v)
    {
        point temp;
        temp.x = this->x + v.x;
        temp.y = this->y + v.y;
        temp.z = this->z + v.z;

        return temp;
    }

    point operator-(point v)
    {
        point temp;
        temp.x = this->x - v.x;
        temp.y = this->y - v.y;
        temp.z = this->z - v.z;

        return temp;
    }

    point &operator=(point v)
    {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;

        return *this;
    }

    point operator*(double d)
    {
        point temp;
        temp.x = x * d;
        temp.y = y * d;
        temp.z = z * d;
        return temp;
    }
    void print()
    {
        cout << x << " " << y << " " << z << endl;
    }

    double dot(point p)
    {
        return (this->x * p.x) + (this->y * p.y) + (this->z * p.z);
    }

    point cross(point p)
    {
        point temp;
        temp.x = (y * p.z - z * p.y);
        temp.y = (z * p.x - x * p.z);
        temp.z = (x * p.y - y * p.x);
        return temp;
    }

    void homogeneous() {
        if(w != 1) {
            x = x/w;
            y = y/w;
            z = z/w;
            w = 1;
        }
    }

    


};

struct triangle{
    point points[3];
    int color[3];
};


struct matrix {
    double val[4][4];

    matrix() {
        for(int i=0;i<4;i++){
            for(int j=0;j<4;j++) {
                val[i][j] = 0;
               // if(i==j) val[i][j] = 1;
            }
        }
        //val[3][3] = 1;
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

double getX(int x,double left_X,double dx) {
    return left_X + x*dx;
}

double getY(int y,double Top_Y,double dy) {
    return Top_Y - y*dy;    
}


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




struct matrix rotate(double angle,point a)
{


    a.normalize();
    angle = (acos(-1.0) / 180.0) * angle;

    point i (1,0,0);
    point j (0,1,0);
    point k (0,0,1);

    point c1 = i * cos(angle) + a * ((a.dot(i)) * (1 - cos(angle))) + a.cross(i) * sin(angle);

    point c2 = j * cos(angle) + a * ((a.dot(j)) * (1 - cos(angle))) + a.cross(j) * sin(angle);

    point c3 = k * cos(angle) + a * ((a.dot(k)) * (1 - cos(angle))) + a.cross(k) * sin(angle);

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

struct point transformPoint(struct matrix T, struct point P) {

    point ret;
    
    ret.x = T.val[0][0]*P.x + T.val[0][1]*P.y + T.val[0][2]*P.z + T.val[0][3]*P.w;
    ret.y = T.val[1][0]*P.x + T.val[1][1]*P.y + T.val[1][2]*P.z + T.val[1][3]*P.w;
    ret.z = T.val[2][0]*P.x + T.val[2][1]*P.y + T.val[2][2]*P.z + T.val[2][3]*P.w;
    ret.w = T.val[3][0]*P.x + T.val[3][1]*P.y + T.val[3][2]*P.z + T.val[3][3]*P.w; 

    return ret;
}




///////// var /////////////

double eyeX, eyeY, eyeZ, lookX, lookY, lookZ, upX, upY, upZ, fovY, aspectRatio, near, far;
string command;

stack<matrix> S;

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

    in1 >> eyeX >> eyeY >> eyeZ >> lookX >> lookY >> lookZ >> upX >> upY >> upZ >> fovY >> aspectRatio >> near >> far;



    while (true)
    {
        in1 >>  command;


        if(command.compare("triangle")==0) {

            point p1, p2, p3;

            in1>> p1.x >> p1.y >> p1.z;
            in1>> p2.x >> p2.y >> p2.z;
            in1>> p3.x >> p3.y >> p3.z;

            point pt1 = transformPoint(S.top(),p1) ;
            point pt2 = transformPoint(S.top(),p2) ;
            point pt3 = transformPoint(S.top(),p3) ;

            pt1.homogeneous();
            pt2.homogeneous();
            pt3.homogeneous();

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

            matrix temp = S.top();
            S.pop();

            S.push( product(temp,Tr) );

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

            matrix temp = S.top();
            S.pop();

            S.push( product(temp,Sc) );


        }

        else if(command.compare("rotate")==0) {

            double angle;
            point axis;

            in1 >> angle >> axis.x >> axis.y >> axis.z;

            matrix Rt = rotate(angle,axis);

            matrix temp = S.top();
            S.pop();

            S.push( product(temp,Rt) );

        }

        else if(command.compare("push")==0) {
            matrix temp = S.top();
            S.push( temp );
        }

        else if(command.compare("pop")==0) {

            S.pop();
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

    point l = (look - eye);
    l.normalize();
    point r = l.cross(up);
    r.normalize();
    point u = r.cross(l);


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

    matrix V = product(R,T);




    while (true)
    {
        

        double val[3][3];


        for(int i=0;i<3;i++) {
            for(int j=0;j<3;j++) {
                if(in2>>val[i][j]) {
                    //////
                }
                else
                    goto out;
            }

        }

        point p1(val[0][0],val[0][1],val[0][2]), p2(val[1][0],val[1][1],val[1][2]), p3(val[2][0],val[2][1],val[2][2]);
        

        point pt1 = transformPoint(V,p1) ;
        point pt2 = transformPoint(V,p2) ;
        point pt3 = transformPoint(V,p3) ;

        pt1.homogeneous();
        pt2.homogeneous();
        pt3.homogeneous();
        ////////////// print ////////////////////


        out2 << pt1.x << "\t" << pt1.y << "\t" << pt1.z << endl;
        out2 << pt2.x << "\t" << pt2.y << "\t" << pt2.z << endl;
        out2 << pt3.x << "\t" << pt3.y << "\t" << pt3.z << endl;
        out2 << endl;      



    }

    out:

    in2.close();
    out2.close();


    //////////////////////////////// STAGE 3 //////////////////////////////////////////////////////////////////

    ifstream in3;
    in3.open("stage2.txt");


    ofstream out3;
    out3.open("stage3.txt");

    out3 << fixed << setprecision(2);

    double fovX = fovY * aspectRatio;
    double tt = near * tan( (fovY*PI) / (2*180) );
    double rr = near * tan( (fovX*PI) / (2*180) );

    matrix P(
        near/rr, 0, 0, 0, 
        0, near/tt, 0, 0, 
        0, 0, -(far+near)/(far-near), -(2*far*near)/(far-near),
        0, 0, -1, 0
    );


    while (true)
    {
        double val[3][3];


        for(int i=0;i<3;i++) {
            for(int j=0;j<3;j++) {
                if(in3>>val[i][j]) {
                    //////
                }
                else
                    goto out2;
            }

        }



        point p1(val[0][0],val[0][1],val[0][2]), p2(val[1][0],val[1][1],val[1][2]), p3(val[2][0],val[2][1],val[2][2]);
        

        point pt1 = transformPoint(P,p1) ;
        point pt2 = transformPoint(P,p2) ;
        point pt3 = transformPoint(P,p3) ;

        pt1.homogeneous();
        pt2.homogeneous();
        pt3.homogeneous();

        ////////////// print ////////////////////


        out3 << pt1.x << "\t" << pt1.y << "\t" << pt1.z << endl;
        out3 << pt2.x << "\t" << pt2.y << "\t" << pt2.z << endl;
        out3 << pt3.x << "\t" << pt3.y << "\t" << pt3.z << endl;
        out3 << endl; 


    }

    out2:
    in3.close();
    out3.close();

    ////////////////////////////////////////////////////////// STAGE 4 //////////////////////////////////////////////////////////////////

    
    ifstream con;
    con.open("config.txt");

    ifstream in4;
    in4.open("stage3.txt");


    ofstream out4;
    out4.open("z_buffer.txt");

    int Screen_Width, Screen_Height;
    double left_X, bottom_Y, lim_Front,lim_Rear;

    con >> Screen_Width >> Screen_Height >> left_X >> bottom_Y >> lim_Front >> lim_Rear;
    int right_X = -left_X;
    int top_Y = -bottom_Y;

    vector<triangle> objects;
    

    while(true) {


        double val[3][3];

        for(int i=0;i<3;i++) {
            for(int j=0;j<3;j++) {
                if(in4>>val[i][j]) {
                    //////
                }
                else
                    goto here;
            }

        }

        triangle t;
        for(int i=0;i<3;i++) {
            t.points[i].setVal(val[i][0],val[i][1],val[i][2]);
            t.color[i] = 100 + rand()%155;
        }
        objects.push_back(t);


        
    }
    here: {}

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    double dx = (right_X-left_X) / Screen_Width;
    double dy = (top_Y-bottom_Y) / Screen_Height;


    double Top_Y = top_Y - dy/2;
    double Left_X = left_X + dx/2;
    


    double** z_buffer = new double*[Screen_Width];
 
    for (int i = 0; i < Screen_Width; i++) {
        z_buffer[i] = new double[Screen_Height];
    }
 
    
    for (int i = 0; i < Screen_Width; i++) {
        for (int j = 0; j < Screen_Height; j++) {
            z_buffer[i][j] = lim_Rear;
        }
    }
 
    bitmap_image* image = new bitmap_image(Screen_Width, Screen_Height);

    // background image color set to black

    for (int i = 0; i < Screen_Width; i++) { 
        for(int j=0;j<Screen_Height;j++) {
            image->set_pixel(i,j,0,0,0);
        }
    }
    
    for (int i=0; i<objects.size(); i++) {

        triangle t = objects[i];

        double t_minY = 9;
        double t_maxY = -9;
        for(int ii=0;ii<3;ii++) {
            t_maxY = max(t_maxY,t.points[ii].y);
            t_minY = min(t_minY,t.points[ii].y);
        }
        
        double top_scanline = min(t_maxY,Top_Y);
        double bottom_scanline = max( t_minY,bottom_Y);

        for(double Yp = top_scanline; Yp >= bottom_scanline; Yp -= dy) {
            
            double Xa=9, Xb=9; // just a random unusual init
            double Za,Zb;
            for(int j=0;j<3;j++) {
                int ii = j; int jj = (j+1) % 3; // two points to check, ii and jj

                // check
                // will get into this loop for two times 
                if ((Yp > t.points[ii].y && Yp < t.points[jj].y) || (Yp > t.points[jj].y && Yp < t.points[ii].y) ) {
                    
                    double ans =  t.points[ii].x - ((t.points[ii].x-t.points[jj].x)*(t.points[ii].y-Yp) / (t.points[ii].y-t.points[jj].y));
                    if(Xa == 9) {
                        Xa = ans;
                        Za = t.points[ii].z - ((t.points[ii].z-t.points[jj].z)*(t.points[ii].y-Yp) / (t.points[ii].y-t.points[jj].y));
                    }
                    else {
                        Xb = ans;
                        Zb = t.points[ii].z - ((t.points[ii].z-t.points[jj].z)*(t.points[ii].y-Yp) / (t.points[ii].y-t.points[jj].y));
                    };
                }
            }


             
            int left_col = round( (min(Xa,Xb) - Left_X) / dx );
            int right_col = round( (max(Xa,Xb) - Left_X) / dx );
            
            int row = (Top_Y-Yp) / dy;
            
            int z_inc = modulusx / (right_col-left_col)
            
            
        }



    }











    image->save_image("output.bmp");

    //Delete the array created
    for(int i=0;i<Screen_Width;i++)    //To delete the inner arrays
    delete [] z_buffer[i];   
    delete [] z_buffer; 
    delete image;

    cout << "done!" << endl;
    return 0;
}