//
//  main.cpp
//  assignment2
//
//  Created by Simon Åkesson on 2018-10-22.
//  Copyright © 2018 Simon Åkesson. All rights reserved.
//

#include <iostream>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include "point.h"
#include <vector>
#include <math.h>
#include "connection.h"
#include <fstream>

using namespace std;
string inputFile = "/Users/simonakesson/Document/ECS 175/assignment2/assignment2/test.txt";
string outputFile = "/Users/simonakesson/Document/ECS 175/assignment2/assignment2/test.txt";
const int INSIDE = 0;
const int LEFT = 0x1;
const int RIGHT = 0x2;
const int BELOW = 0x4;
const int TOP = 0x8;
bool toggleRot = false;
float *PixelBuffer;

GLfloat angle = 1.0f;
int refresh = 30;
int height=600,width=600;
float rx1,rx2,ry1,ry2,rz1,rz2;
float x_1,y_1,z_1,x_2,y_2,z_2;
int select;
bool rot = false,rot1 = false,rot2 = false,rot3 = false;

class polygon{
public:
    vector<point> points2;
    vector<connection> connections;
    void addPoint(point p){points2.push_back(p);}
    void addConnection(int a,int b){connections.push_back(*new connection(a,b));}
    int size = (int)points2.size();
};


vector<polygon> polygons;
vector<polygon> resetPolygons;
int getRegion(point p, int x_max, int y_max, int x_min, int y_min){
    int code = INSIDE;
    if (p.x < x_min)
        code |= LEFT;
    else if (p.x > x_max)
        code |= RIGHT;
    if (p.y < y_min)
        code |= BELOW;
    else if (p.y > y_max)
        code |= TOP;
    
    return code;
}
void drawLines(point p1, point p2){ //XY plane.
    glLineWidth(4.0);
    glBegin(GL_LINES);
    glColor3f(0.0f, 0.5f, 1.0f);
    glVertex2f(p1.x, p1.y);
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex2f(p2.x, p2.y);
    glEnd();
}
void makePix(int x,int y, int red, int green, int blue){
    int pixLoc = (x * 3) + (y * width * 3);
    if(x<width-2 && y < height-2 && x > 0 && y > 0){
        PixelBuffer[pixLoc + 0] = red;
        PixelBuffer[pixLoc + 1] = green;
        PixelBuffer[pixLoc + 2] = blue;
    }
  
}
void drawDDA(point p1, point p2,int r, int g, int b){
    float deltax = p2.x-p1.x, deltay = p2.y-p1.y, x = p1.x, y = p1.y, dist=0, incy, incx;
    if(fabs(deltax) > fabs(deltay)){
        dist = fabs(deltax);
    } else{
        dist  = fabs(deltay);
    }
    incx = deltax/dist;
    incy = deltay/dist;
    makePix(p1.x, p1.y, r, g, b);
    makePix(p2.x, p2.y, r, g, b);
    for(int i=0; i<dist; i++){
        x+=incx;
        y+=incy;
        makePix(x, y, r, g, b);
    }
}
void drawDDAXY(point p1, point p2,int r, int g, int b){
    float deltax = p2.x-p1.x, deltay = p2.y-p1.y, x = p1.x, y = p1.y, dist=0, incy, incx;
    if(fabs(deltax) > fabs(deltay)){
        dist = fabs(deltax);
    } else{
        dist  = fabs(deltay);
    }
    incx = deltax/dist;
    incy = deltay/dist;
     if(x>0 && x<width/2 && y>0 && y<height/2){
    makePix(p1.x, p1.y, r, g, b);
    makePix(p2.x, p2.y, r, g, b);
     }
    for(int i=0; i<dist; i++){
        x+=incx;
        y+=incy;
        if(x>0 && x<width/2 && y>0 && y<height/2){
        makePix(x, y, r, g, b);
        }
    }
}
void drawDDAYZ(point p1, point p2,int r, int g, int b){
    float deltax = p2.y-p1.y, deltay = p2.z-p1.z, x = p1.y, y = p1.z, dist=0, incy, incx;
    if(fabs(deltax) > fabs(deltay)){
        dist = fabs(deltax);
    } else{
        dist  = fabs(deltay);
    }
    incx = deltax/dist;
    incy = deltay/dist;
    if(x>0 && x<width/2 && y>0 && y<height/2){
    makePix(p1.y+width/2, p1.z+height/2, r, g, b);
    makePix(p2.y+width/2, p2.z+height/2, r, g, b);
    }
    for(int i=0; i<dist; i++){
        x+=incx;
        y+=incy;
        if(x>0 && x<width/2 && y>0 && y<height/2){
        makePix(x+height/2, y+height/2, r, g, b);
        }
    }
}
void drawDDAXZ(point p1, point p2,int r, int g, int b){
    float deltax = p2.x-p1.x, deltay = p2.z-p1.z, x = p1.x, y = p1.z, dist=0, incy, incx;
    if(fabs(deltax) > fabs(deltay)){
        dist = fabs(deltax);
    } else{
        dist  = fabs(deltay);
    }
    incx = deltax/dist;
    incy = deltay/dist;
    if(x>0 && x<width/2 && y>0 && y<height/2){
    makePix(p1.x+width/2, p1.z, r, g, b);
    makePix(p2.x+width/2, p2.z, r, g, b);
    }
    for(int i=0; i<dist; i++){
        x+=incx;
        y+=incy;
         if(x>0 && x<width/2 && y>0 && y<height/2){
             
        makePix(x+width/2, y, r, g, b);
         }
    }
}

void translate(vector<point> points,float dx, float dy, float dz, int shape){
    vector <point> temp;
    for(point p : points){
        p.x += dx;
        p.y += dy;
        p.z += dz;
        temp.push_back(*new point(p.x,p.y,p.z));
    }
    polygons[shape].points2 = temp;
}

point findCentre(vector<point> point, int size){
    float minX=point[0].x,minY=point[0].y,maxX=point[0].x,maxY=point[0].y,minZ=point[0].z, maxZ=point[0].z;
    for(int i = 0; i < size; i++){
        if(point[i].x>maxX) {maxX = point[i].x;}
        if(point[i].y>maxY) {maxY = point[i].y;}
        if(point[i].z>maxZ) {maxZ = point[i].z;}
        if(point[i].x<minX) {minX = point[i].x;}
        if(point[i].y<minY) {minY = point[i].y;}
        if(point[i].z<minZ) {minZ = point[i].z;}

    }
    float centreX = (maxX + minX)/2;
    float centreY = (maxY + minY)/2;
    float centreZ = (maxZ + minZ)/2;
    return *new class point(centreX,centreY,centreZ);
}

void scale(vector<point> points,float factor,int shape){
    point centre = findCentre(points, (int)points.size());
    vector <point> temp;
    for(point p : points){
        p.x = (p.x-centre.x)*factor + centre.x;
        p.y = (p.y-centre.y)*factor + centre.y;
        p.z = (p.z-centre.z)*factor + centre.z;
        temp.push_back(*new point(p.x,p.y,p.z));
    }
    polygons[shape].points2 = temp;
}

void rotateX(vector<point> points, float factor, int shape){
    point centre = findCentre(points, (int)points.size());
    double angleDegrees = factor * M_PI / 180.0;
    vector<point>temp;
    for(point p : points){
        float dx = p.x, dy = p.y, dz = p.z;
        float tempY = p.y;
        tempY = (cos(angleDegrees)*(dy-centre.y) - sin(angleDegrees)*(dz-centre.z));
        dz = (sin(angleDegrees)*(dy-centre.y) + cos(angleDegrees)*(dz-centre.z));
        dy = tempY + centre.y;
        dz += centre.z;
        temp.push_back(*new point(dx,dy,dz));
    }
    polygons[shape].points2 = temp;
}

void rotateY(vector<point> points,float factor, int shape){
    point centre = findCentre(points, (int)points.size());
    double angleDegrees = factor * M_PI / 180.0;
    vector<point>temp;
    for(point p : points){
        float dx = p.x, dy = p.y, dz = p.z;
        float tempX = p.x;
        tempX = (cos(angleDegrees)*(dx-centre.x)) + (sin(angleDegrees)*(dz-centre.z));
        dz = (cos(angleDegrees)*(dz-centre.z))-(sin(angleDegrees)*(dx-centre.x));
        dx = tempX + centre.x;
        dz += centre.z;
        temp.push_back(*new point(dx,dy,dz));
    }
    polygons[shape].points2 = temp;
}

void rotateZ(vector<point> points, float factor, int shape){
    point centre = findCentre(points, (int)points.size());
    double angleDegrees = factor * M_PI / 180.0;
    vector<point>temp;
    for(point p : points){
        float dx = p.x, dy = p.y, dz = p.z;
        float tempX = p.x;
        tempX = (cos(angleDegrees)*(dx-centre.x) - sin(angleDegrees)*(dy-centre.y));
        dy = (sin(angleDegrees)*(dx-centre.x) + cos(angleDegrees)*(dy-centre.y));
        dx = tempX + centre.x;
        dy += centre.y;
        temp.push_back(*new point(dx,dy,dz));
    }
    polygons[shape].points2 = temp;
}

void cohenSutherland(point p1, point p2,int x_max, int y_max, int x_min, int y_min, int r,int g, int b){
        int code1 = getRegion(p1, x_max,  y_max,  x_min,  y_min);
        int code2 = getRegion(p2, x_max,  y_max,  x_min,  y_min);
        bool done = false;
        while(true){
            if((code1 == 0) && (code2 == 0)){
                done = true;
                break;
            }else if (code1&code2){
                break;}
            else{
                int tempCode = 0;
                float x = 0.0 ,y = 0.0;
                (code1 != 0) ? tempCode = code1 : tempCode = code2;
                if(tempCode & LEFT){
                    x = x_min;
                    y = (p2.y-p1.y)*(x_min-p1.x)/(p2.x-p1.x) + p1.y;
                }else if (tempCode & RIGHT){
                    x = x_max;
                    y = (p2.y-p1.y)*(x_max-p1.x)/(p2.x-p1.x) + p1.y;
                }else if (tempCode & BELOW){
                    y = y_min;
                    x = (p2.x-p1.x)*(y_min-p1.y)/(p2.y-p1.y) + p1.x;
                }else if (tempCode & TOP){
                    x = (p2.x-p1.x)*(y_max-p1.y)/(p2.y-p1.y) + p1.x;
                    y = y_max;
                }
                if(code1 == tempCode){
                    p1.y = y;
                    p1.x = x;
                    code1 = getRegion(p1, x_max,  y_max,  x_min,  y_min);
                }else{
                    p2.y = y;
                    p2.x = x;
                    code2 = getRegion(p2, x_max,  y_max,  x_min,  y_min);
                }
            }
        }
    if(done){
        drawDDA(p1, p2, r, g, b);
    }
}

void drawAxis(point p1, point p2, bool erase){
    float k1 = (p2.y-p1.y)/(p2.x-p1.x);
    float m1 = p1.y - k1*p1.x;
    float y1 = k1*width/2 +m1;
    point tp1 = *new point(0,m1,0);
    point tp2 = *new point(width/2,y1,0);
    cohenSutherland(tp1, tp2,width/2,height/2, 0,0,1,1,0);
    if(erase){
        cohenSutherland(tp1, tp2,width/2,height/2, 0,0,0,0,0);
    }
    float k2 = (p2.z-p1.z)/(p2.x-p1.x);
    float m2 = p1.z - k2*p1.x;
    float z2 = k2*width/2 +m2;
    point tp1xz = *new point(width/2,m2,0);
    point tp2xz = *new point(width,z2,0);
    cohenSutherland(tp1xz, tp2xz,width,width/2,width/2, 0,0,0,1);
    if(erase){
         cohenSutherland(tp1xz, tp2xz,width,width/2,width/2, 0,0,0,0);
    }
    float k3 = (p2.z-p1.z)/(p2.y-p1.y);
    float m3 = p1.z - k3*p1.y;
    float z3 = k3*width/2 +m3;
    point tp1yz = *new point(width/2,m3+height/2,0);
    point tp2yz = *new point(width,z3+height/2,0);
    cohenSutherland(tp1yz, tp2yz,width,width,width/2, width/2,1,0,1);
    if(erase){
        cohenSutherland(tp1yz, tp2yz,width,width,width/2, width/2,0,0,0);
    }
}
void rotationArb(point p1, point p2, float factor,vector<point> points, int shape){
     point centre = findCentre(points, (int)points.size());
    drawAxis(p1, p2,false);
    float denom = sqrt(pow((p2.x-p1.x), 2) + pow((p2.y-p1.y), 2) + pow((p2.z-p1.z), 2));
    float dx = (p2.x-p1.x)/denom;
    float dy = (p2.y-p1.y)/denom;
    float dz = (p2.z-p1.z)/denom;
    centre.x -=p1.x; centre.y-=p1.y; centre.z-=p1.z;
    double angleDegrees = factor * M_PI / 180.0;
    float l = sqrt(pow(dy, 2) + pow(dz, 2)); //m21
    vector<point>temp;
    for(point p : points){
        float tempX=p.x,tempY=p.y,tempZ=p.z;
        tempX -=p1.x; tempY-=p1.y; tempZ-=p1.z; //translate to p0
        float tempMult = tempY;   // m21
        tempY = tempY*dz/l - tempZ*dy/l;    //m21
        tempZ = tempMult*dy/l + tempZ*dz/l; //m21
        
        tempMult = tempX; //m22
        tempX = tempX*l - tempZ*dx; // m22
        tempZ = tempMult*dx + tempZ*l; // m22
        
        tempMult = tempX;  //m3
        tempX = (cos(angleDegrees)*(tempX) - sin(angleDegrees)*(tempY)); //m3
        tempY = (sin(angleDegrees)*(tempMult) + cos(angleDegrees)*(tempY)); //m3
    
        tempMult = tempX;
        tempX = tempX*l + tempZ*dx; // m22(-1)
        tempZ = l*tempZ-tempMult*dx; //m22(-1)
        
        tempMult = tempY;
        tempY = tempY*dz/l + tempZ*dy/l;    //m21
        tempZ =  tempZ*dz/l - tempMult*dy/l; //m21
        
        tempX +=p1.x; tempY+=p1.y; tempZ+=p1.z;
        temp.push_back(*new point(tempX,tempY,tempZ));
    }
        
    polygons[shape].points2 = temp;
    
}
void timer (int value){
    glutPostRedisplay();
    glutTimerFunc(refresh, timer, 0);
}
void setup(){
    drawDDA(*new point(0,height/2,0), *new point(width,height/2,0), 1, 1, 1);
    drawDDA(*new point(width/2,0,0), *new point(width/2,height,0), 1, 1, 1);
}
void clearPolygonXY(){
    for(polygon poly : polygons){
        for(connection conn : poly.connections){
            drawDDAXY(poly.points2[conn.a], poly.points2[conn.b], 0, 0, 0);
            drawDDAXZ(poly.points2[conn.a], poly.points2[conn.b], 0, 0, 0);
            drawDDAYZ(poly.points2[conn.a], poly.points2[conn.b], 0, 0, 0);
        }
    }
}
void drawPolygonXY(){
    for(polygon poly : polygons){
        for(connection conn : poly.connections){
            drawDDAXY(poly.points2[conn.a], poly.points2[conn.b], 1, 0, 1);
            drawDDAXZ(poly.points2[conn.a], poly.points2[conn.b], 0, 1, 1);
            drawDDAYZ(poly.points2[conn.a], poly.points2[conn.b], 0, 1, 0);
        }
    }
}
//reads the file of the original input.
void readFile(){
    
    std::ifstream file(inputFile);
    std::string str;
    int nbrPoly;
    int nbrConn;
    if (!file) {
        cout << "Unable to open file";
        exit(1);
    }
    int counter = 0;
    polygon poly;
    std::string delimiter = " ";
    while(std::getline(file, str)){
        if(counter==0) {
            nbrPoly = stoi(str);
            counter++;
        }
        else if(str == ""){
            for(int i = 0; i<poly.points2.size(); i++){
            }
            if(poly.points2.size() !=0){
                polygons.push_back(poly);
            }
            poly = *new polygon();
        }
        else if(str.length() < 3){
            if(counter==1){
                poly.size=stoi(str);
                counter++;
            }else{
                nbrConn = stoi(str);
            }
        }
        else if(str.length() == 3){
                std::string con1 = str.substr(0, str.find(delimiter));
                std::string con2 = str.substr(con1.length()+1);
                poly.addConnection(stoi(con1)-1,stoi(con2)-1);
            }else{
                //std::cout<<str.length() << " \n";
                //std::cout<<str.find(" ") << " \n";
            std::string xcord = str.substr(0, str.find(delimiter));
            //std::string ycord = str.substr(xcord.length()+1,str.find(" "));
            std::string temps = str.substr(xcord.length()+1);
            //    std::cout<<temps << " \n";
            std::string ycord = temps.substr(0, temps.find(delimiter));
               // std::cout << zcord << " \n";
            std::string zcord = str.substr(ycord.length() + xcord.length()+1);
                std::cout << zcord << " \n";
            int winDelta = max(width,height);
            float ndcX = (((stof(xcord))*winDelta/2)+1);
            float ndcY = (((stof(ycord))*winDelta/2)+1);
            float ndcZ = (((stof(zcord))*winDelta/2)+1);
            point p = *new point(ndcX,ndcY,ndcZ);
            poly.addPoint(p);
        }
    }
    resetPolygons = polygons;
}
void keyboard(int key, int x, int y){
    int selection;
    switch(key){
        case '0' : {
            ofstream file;
            file.open (outputFile);
            file << polygons.size() << "\n\n";
            for(polygon p : polygons){
                file << p.points2.size()<< "\n";
                for(point points : p.points2){
                    file << (float)(points.x-1)/(width/2) << " "<< (float)(points.y-1)/(width/2) << " " << (float)(points.z-1)/(width/2) <<  "\n" ;
                }
                file << p.connections.size()<< "\n";
                for(connection conn : p.connections){
                    file << conn.a+1<< " "<< conn.b+1 << "\n" ;
                }
                file << "\n";
            }
            file.close();
            exit(0);
            break;
        }
        case '1' :{
            float xcord;
            float ycord;
            float zcord;
            std::cout <<"Enter which polygon \n";
            cin>>selection;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>selection;
            }
            std::cout <<"Enter new x y z coordinates to move in each direction \n";
            cin>>xcord >> ycord >> zcord;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>xcord>>ycord>>zcord;
            }
            int winDelta = max(width,height);
            float ndcX = ((xcord*winDelta/2)+1);
            float ndcY = ((ycord*winDelta/2)+1);
            float ndcZ = ((zcord*winDelta/2)+1);
            clearPolygonXY();
            translate(polygons[selection-1].points2, ndcX, ndcY, ndcZ, selection-1);
            break;
    }
        case '2' : {
            float factor;
            std::cout <<"Enter which polygon \n";
            cin>>selection;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>selection;
            }
            std::cout <<"Enter what factor to scale with \n";
            cin>>factor;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>factor;
            }
            clearPolygonXY();
            scale(polygons[selection-1].points2, factor, selection-1);
            break;
        }
        case '3' :{
            float angle;
            float x1,y1,z1,x2,y2,z2;
            if(rot==true){
                drawAxis(*new point(x_1,y_1,z_1), *new point(x_2,y_2,z_2), true);
                rot = false;
                break;
            }
            std::cout <<"Enter which polygon \n";
            cin>>selection;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>selection;
            }
            std::cout <<"Enter the first point x1 y1 z1 \n";
            cin>>x1>>y1>>z1;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>x1>>y1>>z1;
            }
            std::cout <<"Enter the second point x2 y2 z2 \n";
            cin>>x2>>y2>>z2;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>x2>>y2>>z2;
            }
            std::cout <<"Enter angle of rotation \n";
            cin>>angle;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>angle;
            }
            int winDelta = max(width,height);
            x_1 = ((x1*winDelta/2)+1);
            y_1 = ((y1*winDelta/2)+1);
            z_1 = ((z1*winDelta/2)+1);
            x_2 = ((x2*winDelta/2)+1);
            y_2 = ((y2*winDelta/2)+1);
            z_2 = ((z2*winDelta/2)+1);
            clearPolygonXY();
            rotationArb(*new point(x_1,y_1,z_1), *new point(x_2,y_2,z_2), angle, polygons[selection-1].points2, selection-1);
            rot=true;
            break;
        }
        case '4' : {
            float tx1,ty1,tz1,tx2,ty2,tz2;
            if(toggleRot==true){
                toggleRot = false;
                drawAxis(*new point(rx1,ry1,rz1), *new point(rx2,ry2,rz2), true);
                break;
            }
            std::cout <<"Enter which polygon \n";
            cin>>select;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>select;
            }
            std::cout <<"Enter the first point x1 y1 z1 \n";
            cin>>tx1>>ty1>>tz1;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>tx1>>ty1>>tz1;
            }
            std::cout <<"Enter the second point x2 y2 z2 \n";
            cin>>tx2>>ty2>>tz2;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>tx2>>ty2>>tz2;
            }
            int winDelta = max(width,height);
            rx1 = ((tx1*winDelta/2)+1);
            ry1 = ((ty1*winDelta/2)+1);
            rz1 = ((tz1*winDelta/2)+1);
            rx2 = ((tx2*winDelta/2)+1);
            ry2 = ((ty2*winDelta/2)+1);
            rz2 = ((tz2*winDelta/2)+1);
            clearPolygonXY();
            select = select-1;
            toggleRot=true;
            break;
        }
        case '5' : {
            drawAxis(*new point(x_1,y_1,z_1), *new point(x_2,y_2,z_2), true);
            clearPolygonXY();
            polygons = resetPolygons;
            break;
        }
        case '6': {
            float ang;
            int axis;
            std::cout <<"Enter which polygon \n";
            cin>>selection;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>selection;
            }
            std::cout <<"which axis (choose number) \n 1: x \n 2: y \n 3: z \n";
            cin>>axis;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>axis;
            }
            std::cout <<"Enter angle of rotation \n";
            cin>>ang;
            if(!cin){
                cout << "REENTER!  \n";
                cin.clear();
                cin.ignore((numeric_limits<streamsize>::max)(), '\n');
                cin>>ang;
            }
            switch(axis){
                case 1 : clearPolygonXY(); rotateX(polygons[selection-1].points2, ang, selection-1); break;
                case 2 : clearPolygonXY(); rotateY(polygons[selection-1].points2, ang, selection-1); break;
                case 3 : clearPolygonXY(); rotateZ(polygons[selection-1].points2, ang, selection-1); break;
                default : std::cout << "invalid entry \n"; break;
            }
            break;
        }
        case '7' :{
            rot1 ? rot1 = false : rot1 = true;
            break;
        }
        case '8' :{
            rot2 ? rot2 = false : rot2 = true;
            break;
        }
        case '9' :{
            rot3 ? rot3 = false : rot3 = true;
            break;
        }
    }
    glutPostRedisplay();
     std::cout<<"Select option 1-6 but writing number. \n 1: Translate \n 2: Scale \n 3: Rotate polygon by arbitrary axis \n 4: Continuously rotate the polygon around arbitrary axis. \n 5: Reset all transforms since starting program \n 6: Spin w.r.t. x-,y-,z axis of your choice \n 7: Spin all w.r.t. x-axis  \n 8: Spin all w.r.t. y-axis \n 9: Spin all w.r.t. z-axis  \n 0: Save and exit program \n";
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    clearPolygonXY();
    if(toggleRot){
        rotationArb(*new point(rx1,ry1,rz1), *new point(rx2,ry2,rz2), 1, polygons[select].points2,select);
        }
    if(rot1){
        for(int i = 0; i<polygons.size(); i++){
        rotateX(polygons[i].points2, 1, i);
        }
    }
    if(rot2){
        for(int i = 0; i<polygons.size(); i++){
            rotateY(polygons[i].points2, 1, i);
        }
    }
    if(rot3){
        for(int i = 0; i<polygons.size(); i++){
            rotateZ(polygons[i].points2, 1, i);
        }
    }
    
    drawPolygonXY();
    setup();
    
    glDrawPixels(height, width, GL_RGB, GL_FLOAT, PixelBuffer);
    glFlush();  // Render now
    glutSwapBuffers ();
}
/* Main function: GLUT runs as a console application starting at main()  */
int main(int argc, char** argv) {
    PixelBuffer = new float[height * width * 3];
    readFile();
    std::cout<<"Select option 1-6 but writing number. \n 1: Translate \n 2: Scale \n 3: Rotate polygon by arbitrary axis \n 4: Continuously rotate the polygon around arbitrary axis.\n 5: Reset all transforms since starting program \n 6: Spin w.r.t. x-,y-,z axis of your choice \n 7: Spin all w.r.t. x-axis  \n 8: Spin all w.r.t. y-axis \n 9: Spin all w.r.t. z-axis  \n 0: Save and exit program \n";

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE);
    glutInitWindowSize(width, height);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Assignment 2 ");
    glClearColor(0, 0, 0, 0);
    glutDisplayFunc(display);
    glutSpecialFunc(keyboard);
    glutTimerFunc(0, timer, 0);
    glutMainLoop();
    return 0;
}


