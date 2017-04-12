/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: <Your name here>
 * *************************
 */

#ifdef WIN32
#include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
#include <GL/gl.h>
#include <GL/glut.h>
#elif defined(__APPLE__)
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <math.h>
#include <algorithm>
#include "../external/glm/glm/gtc/type_ptr.hpp"


#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

using namespace std;

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
};

struct Triangle
{
    Vertex v[3];
};

struct Sphere
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
};

struct Light
{
    double position[3];
    double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

double x_min;
double y_min;
double x_max;
double y_max;

//this stores all the ray emitted from the camera.
vector<glm::vec3> directions;
double points[WIDTH*HEIGHT][3];//pixels to depict
double colors[WIDTH*HEIGHT][3];//colors of the pixels;
//this function return a normalized vector.
//vector<double> normalized(double x, double y, double z);
glm::vec3 normalized(double x, double y, double z);
//functions by alexis

//this function returns a vector of center of pixel positions:
vector<glm::vec3 > returnRayDirections(){
    
    double aspectRatio = 1.0*WIDTH/HEIGHT;
    double z = -1.0f;
    double angle = M_PI/180*fov;
    x_min = - aspectRatio * tan(angle/2);
    y_min = -tan(angle/2);
    x_max = aspectRatio * tan(angle/2);
    y_max = tan(angle/2);
    
    //    cout<<fov/2<<endl;
    //    cout<<tan(M_PI/180*fov/2)<<endl;
    cout<<x_min<<" "<<x_max<<" "<<y_min<<" "<<y_max<<endl;
    //    cout<<aspectRatio<<endl;
    //    cout<<z<<endl;
    
    vector<glm::vec3 > pixelCenter;
    
    double increment_x = 1.0*(x_max-x_min)/(WIDTH-1);
    double increment_y = 1.0*(y_max-y_min)/(HEIGHT-1);
    
//    int counterx = 0;
    for(double i = x_min; i < x_max+0.5*increment_x; i+=increment_x){
        //        counterx++;
        for(double j = y_min; j < y_max+0.5*increment_y; j+=increment_y){
            //            array<double, 3> eachPixel = {i, j, z};
            
            glm::vec3 eachPixel = glm::normalize(glm::vec3((float)i, (float)j, (float)z));
//            glm::vec3 eachPixel = normalized(i, j, z);
            pixelCenter.push_back(eachPixel);
        }
    }
    //    cout<<counterx<<" "<<increment_x<<" "<<increment_y<<" look at this"<<endl;
    cout<<"this is the length: "<<pixelCenter.size()<<" "<<HEIGHT*WIDTH<<endl;
    //    exit(0);
    return pixelCenter;
    
}


glm::vec3 normalized(double x, double y, double z){
    double length = sqrt( pow(x, 2) + pow(y, 2) + pow(z, 2));
    glm::vec3 result;
    if(length == 0){
        cout<<"the length is 0"<<endl;
        result = glm::vec3(0,0,0);
        return result;
    }else{
        result = glm::vec3(x/length,y/length,z/length);
        return result;
    }
}

//returns a normal vector of two given vectors.
vector<double> crossProduct(double v1_x, double v1_y, double v1_z, double v2_x, double v2_y, double v2_z){
    //    vector<double> result = normalized( v1_y*v2_z - v1_z*v2_y, v1_z*v2_x - v1_x*v2_z, v1_x*v2_y - v1_y*v2_x );
    vector<double> result;
    result.push_back(v1_y*v2_z - v1_z*v2_y);
    result.push_back(v1_z*v2_x - v1_x*v2_z);
    result.push_back(v1_x*v2_y - v1_y*v2_x);
    return result;
}

//returns the area given three vertex of a triangle.
double calculateArea(double v1_x, double v1_y, double v1_z, double v2_x, double v2_y, double v2_z, double v3_x, double v3_y, double v3_z ){
//    vector<double> cross_product = crossProduct(v2_x-v1_x, v2_y-v1_y, v2_z-v1_z, v3_x-v1_x, v3_y-v1_y, v3_z-v1_z);
//    //    cout<<"this is the result of a cross product"<< cross_product[0]<<" "<< cross_product[1]<<" "<< cross_product[2]<<" "<<endl;
//    double area = 0.5*sqrt( pow(cross_product[0], 2) + pow(cross_product[1], 2) + pow(cross_product[2], 2) );
//    //    cout<<"this is area from inside"<<area<<endl;
//    //    if(area < 0) {
//    //        area *= -1.0f;
//    //
//    //    }
//    return area;
        //projecting onto 2D plane:
        double area;
        //step1: calculating normal of the plan
        vector<double> normal = crossProduct(v2_x-v1_x, v2_y-v1_y, v2_z-v1_z, v3_x-v1_x, v3_y-v1_y, v3_z-v1_z);
        //step2: compare the longest vector
        double longest_value = 0;
        int longest_axis_of_normal = 0;
        for(int i = 0; i < 3; i++){
            if(abs(normal[i])>longest_value){
                longest_axis_of_normal = i;
                longest_value = abs(normal[i]);
            }
        }
        //step3: project onto the other two planes:
    //    int axis_x = (i+i)%3;
    //    int axis_y = (i+2)%3;
        if(longest_axis_of_normal==0){
            //project onto yz plane:
    //        area = abs( 0.5 * ( (v2_y-v1_y)*(v3_z-v1_z) - (v3_y-v1_y)*(v2_z-v1_z) ) );
            area =  0.5 * ( (v2_y-v1_y)*(v3_z-v1_z) - (v3_y-v1_y)*(v2_z-v1_z) ) ;
    
        }else if(longest_axis_of_normal==1){
            //project onto xz plane:
    //        area = abs( 0.5 * ( (v2_x-v1_x)*(v3_z-v1_z) - (v3_x-v1_x)*(v2_z-v1_z) ) );
            area =  0.5 * ( (v2_x-v1_x)*(v3_z-v1_z) - (v3_x-v1_x)*(v2_z-v1_z) ) ;
    
        }else if(longest_axis_of_normal==2){
            //project onto xy plane:
    //        area = abs( 0.5 * ( (v2_x-v1_x)*(v3_y-v1_y) - (v3_x-v1_x)*(v2_y-v1_y) ) );
            area =  0.5 * ( (v2_x-v1_x)*(v3_y-v1_y) - (v3_x-v1_x)*(v2_y-v1_y) ) ;
    
        }
        return area;
}
//this function takes in p0,  unit light direction, and index of the sphere
//returns -1 if there is no intersection, min(t0, t1) if intersection occurs.
int intersectSphere(double x0,double y0, double z0, double xd, double yd, double zd, int i){
    //    spheres stores all the spheres;
    //    cout<<"this is the size of spheres"<<num_spheres<<endl;
        double center[3] = {spheres[i].position[0], spheres[i].position[1], spheres[i].position[2]};
        //        cout<< "here is the sphere center"<<center[0]<<" "<<center[1]<<" "<<center[2]<<endl;
        double r = spheres[i].radius;
        double xc = spheres[i].position[0];
        double yc = spheres[i].position[1];
        double zc = spheres[i].position[2];

        double a =  1;
        double b = 2*(xd *(x0-xc) + yd*(y0-yc) + zd*(z0-zc) );
        double c = pow(x0-xc,2)+pow(y0-yc, 2)+pow(z0-zc, 2)-pow(r,2);
        
        double delta = pow(b,2)-4*c;
        //if intersection occurs
        if(delta>=0){
            double t0 = (-b-sqrt(delta))/2;
            double t1 = (-b+sqrt(delta))/2;
            if(t0 > 0){
//                cout<<"this is t0"<<t0<<endl;
                return t0;
            }
        }
        //if there are no intersection then ignore
        return -1;

}
//takes in p0, unit direction vector, and triangle index i;
//returns (0,0,0) if there is no intersection, returns the position of intersection if intersectino occurs
vector<double> intersectTriangle(double x0, double y0, double z0, double xd, double yd, double zd, int i){
    vector<double> result;
    
        vector<double> v1;//c0
        vector<double> v2;//c1
        vector<double> v3;//c2
        for(int j = 0; j < 3; j++){
            v1.push_back(triangles[i].v[0].position[j]);
            v2.push_back(triangles[i].v[1].position[j]);
            v3.push_back(triangles[i].v[2].position[j]);
        }
        
        double vec1[3] = {v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2]};
        double vec2[3] = {v3[0]-v1[0], v3[1]-v1[1], v3[2]-v1[2]};
        vector<double> temp = crossProduct(vec1[0], vec1[1], vec1[2], vec2[0], vec2[1], vec2[2]);
        glm::vec3 normal = normalized(temp[0], temp[1], temp[2]);
        
        //    vector<double> normal = normalized( vec1[1]*vec2[2]-vec1[2]*vec2[1], vec1[2]*vec2[0]-vec1[0]*vec2[2], vec1[0]*vec2[1]-vec1[1]*vec2[0] );
        
        double a = normal[0];
        double b = normal[1];
        double c = normal[2];
        
        double p_d = -(a*v1[0] + b*v1[1] + c*v1[2]); //d = -ax-by-cz;
        double t = -(a*x0+b*y0+c*z0+p_d)/(a*xd + b*yd + c*zd);
        
        //when intersected
        if(t>0){
            //            cout<<"intersecting with plane"<<endl;
            vector<double> c;
            c.push_back(x0 + xd*t);
            c.push_back(y0 + yd*t);
            c.push_back(z0 + zd*t);
            
            
            //cc1c2/c0c1c2
            double alpha = calculateArea(c[0], c[1], c[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2])/calculateArea(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2]);
            //c0cc2/c0c1c2
            double beta = calculateArea(v1[0], v1[1], v1[2], c[0], c[1], c[2], v3[0], v3[1], v3[2])/calculateArea(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2]);
            //c0c1c/c0c1c2
            double gamma = calculateArea(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], c[0], c[1], c[2])/calculateArea(v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2]);
            
            //cout << alpha << ' ' << beta << ' ' << gamma << endl;
            //cin.get();
            
            if(gamma  <= 1.0 &&  alpha <= 1.0 && beta <= 1.0 &&
               gamma >= 0.0 &&  alpha >= 0.0 && beta >= 0.0){
                result.push_back(alpha);
                result.push_back(beta);
                result.push_back(gamma);
//                result.push_back(alpha*v1[0]+beta*v2[0]+gamma*v3[0]);
//                result.push_back(alpha*v1[1]+beta*v2[1]+gamma*v3[1]);
//                result.push_back(alpha*v1[2]+beta*v2[2]+gamma*v3[2]);
            }else{
                for(int j = 0; j <3; j++){
                    result.push_back(0);
                }
            }
        }else{
            for(int j = 0; j <3; j++){
                result.push_back(0);
            }
        }
    return result;
    
    
}

//fire ray from the intersection point to the light source
//takes in a position of the intersection point, the index of the current sphere
//the index of light, a boolean value of if it is a sphere
//return true if is a shadow,
//return flase if is not a shadow
//to become a shadow: intersect with a geometry, the intersection point is not itself
bool isShadow(double px, double py, double pz, int current, int i, bool sphere){
    //back light direction:
    double lx = lights[i].position[0];
    double ly = lights[i].position[1];
    double lz = lights[i].position[2];

    glm::vec3 light_position = glm::vec3(lx, ly, lz);
    glm::vec3 intersection_position = glm::vec3(px, py,pz);
    glm::vec3 l = light_position-intersection_position;
    //getting normalized l vecotr
    l = normalized(l.x, l.y, l.z);
    for(int j = 0; j < num_spheres; j++){
        double block = intersectSphere(px, py, pz, l.x, l.y, l.z, j);
        //if intersectin occurs:
        //it will be a shadow if:
        //it intersects with other object;
        //the nearest intersection is not itself:
        if(block != -1){
            if(!sphere || current != j) return true;
            glm::vec3 block_point = intersection_position + (float)block * l;
            if(block_point!=intersection_position)return true;
        }
        //if no intersection, it will not be a shadow
    }
    for(int j = 0; j < num_triangles; j++){
        vector<double> block = intersectTriangle(px, py, pz,  l.x, l.y, l.z, j);
        if(block[0] != 0 && block[1] != 0 && block[2] != 0){
            if(sphere || current != j) return true;
//            glm::vec3 block_point = intersection_position + (float)block * l;
//            if(block_point!=intersection_position)return true;
        }
    }
    //here should calculate phone model
    return false;

}

glm::vec3 getTriangleIntersection(int index_of_object, glm::vec3 p){
    glm::vec3 v1 = glm::vec3(triangles[index_of_object].v[0].position[0],triangles[index_of_object].v[0].position[1],triangles[index_of_object].v[0].position[2]);
    glm::vec3 v2 = glm::vec3(triangles[index_of_object].v[1].position[0],triangles[index_of_object].v[1].position[1],triangles[index_of_object].v[1].position[2]);
    glm::vec3 v3 = glm::vec3(triangles[index_of_object].v[2].position[0],triangles[index_of_object].v[2].position[1],triangles[index_of_object].v[2].position[2]);
    glm::vec3 pos = p.x*v1 + p.y*v2 + p.z*v3; //position = alpha*v1 +beta*v2+gamma*v3
    return pos;
}
//this function calculate the color of a point
//input: glm::vec3: the intersection point(if sphere)OR alpha_beta_gamma value(if triangle), the light direction the shape of the object, the index of object and light
//output: glm::vec3: calculated color
glm::vec3 phongModel(glm::vec3 p, string shape, int index_of_object, int index_of_light){
    //getting the unit vector l
    glm::vec3 light_position = glm::vec3(lights[index_of_light].position[0], lights[index_of_light].position[1], lights[index_of_light].position[2]);
    glm::vec3 light_color = glm::vec3(lights[index_of_light].color[0],lights[index_of_light].color[1], lights[index_of_light].color[2]);
    glm::vec3 result_color;
    
    //calculate sphere color
    if(shape=="sphere"){
//        double position[3];
//        double color_diffuse[3];
//        double color_specular[3];
//        double shininess;
//        double radius;
   
//        double position[3];
//        double color[3];
        glm::vec3 l = glm::normalize(light_position-p);
        glm::vec3 v = glm::normalize(-1.0f * p);
        glm::vec3 sphere_center = glm::vec3(spheres[index_of_object].position[0], spheres[index_of_object].position[1], spheres[index_of_object].position[2]);
        glm::vec3 n =glm::normalize( p-sphere_center);
        
        //cout << n.x << ' ' << n.y << ' ' << n.z << endl;
        if(n.z<0){
            cout<<"something wrong"<<endl;
        }
        float dot_l_n = glm::dot(l,n);
        if(dot_l_n < 0) dot_l_n = 0;
        if(dot_l_n > 1) dot_l_n = 1;

        glm::vec3 r = glm::normalize(2*glm::dot(l,n)*n-l);

        glm::vec3 kd = glm::vec3(spheres[index_of_object].color_diffuse[0], spheres[index_of_object].color_diffuse[1], spheres[index_of_object].color_diffuse[2]);
        glm::vec3 ks = glm::vec3(spheres[index_of_object].color_specular[0],spheres[index_of_object].color_specular[1], spheres[index_of_object].color_specular[2]);
        
        float dot_r_v = glm::dot(r, v);
        if(dot_r_v < 0) dot_r_v =0;
        if(dot_r_v > 1) dot_r_v =1;
//        glm::vec3 object_color = glm::vec3(1, 1, 1) * p.z;// glm::dot(glm::vec3(0.0f, 0.0f, 1.0f), n);//kd * dot_l_n;// + ks * (float)pow(dot_r_v, spheres[index_of_object].shininess);
//        glm::vec3 object_color = kd * dot_l_n + ks * (float)(pow(dot_r_v, spheres[index_of_object].shininess));
        double obj_x = kd.x * dot_l_n + ks.x * (pow(dot_r_v, spheres[index_of_object].shininess));
        double obj_y = kd.y * dot_l_n + ks.y * (pow(dot_r_v, spheres[index_of_object].shininess));
        double obj_z = kd.z * dot_l_n + ks.z * (pow(dot_r_v, spheres[index_of_object].shininess));
        if(obj_x*light_color.x>1||obj_x*light_color.x<0||obj_y*light_color.y>1||obj_y*light_color.y<0||obj_z*light_color.z>1||obj_z*light_color.z<0){
            cout<<"what the fuck"<<endl;
        }

//        glm::vec3 object_color = kd * dot_l_n + ks * (float)(pow(dot_r_v, spheres[index_of_object].shininess));
        glm::vec3 object_color = glm::vec3((float)obj_x, (float)obj_y, (float)obj_z);
        result_color = glm::vec3(light_color.x * object_color.x, light_color.y * object_color.y, light_color.z * object_color.z);
        
    }
    //calculate color for a triangle; now p is the alpha beta gamma value
    else if(shape =="triangle"){
        
//        struct Vertex
//        {
//            double position[3];
//            double color_diffuse[3];
//            double color_specular[3];
//            double normal[3];
//            double shininess;
//        };
//        
//        struct Triangle
//        {
//            Vertex v[3];
//        };


        glm::vec3 pos = getTriangleIntersection(index_of_object, p);
       
        glm::vec3 l = light_position-pos;
        l = normalized(l.x, l.y, l.z);//normalized light vector
        glm::vec3 v = -1.0f * pos;
        v = normalized(v.x, v.y, v.z);//normalized view vector
//        glm::vec3 vec1 = v2-v1;
//        glm::vec3 vec2 = v3-v1;
//        glm::vec3 n = cross(vec1, vec2);
        glm::vec3 n1 = glm::vec3(triangles[index_of_object].v[0].normal[0],triangles[index_of_object].v[0].normal[1],triangles[index_of_object].v[0].normal[2]);
        glm::vec3 n2 = glm::vec3(triangles[index_of_object].v[1].normal[0],triangles[index_of_object].v[1].normal[1],triangles[index_of_object].v[1].normal[2]);
        glm::vec3 n3 = glm::vec3(triangles[index_of_object].v[2].normal[0],triangles[index_of_object].v[2].normal[1],triangles[index_of_object].v[2].normal[2]);
        
        glm::vec3 n = p.x*n1+p.y*n2+p.z*n3;
        n = normalized(n.x, n.y, n.z);
//        if(dot(v, n)<0) n = -1.0f*n;//normalized n vector
        glm::vec3 r = 2*glm::dot(l, n)*n-l;
        r = normalized(r.x, r.y, r.z);
        //the diffuse color of the three points
        glm::vec3 v1_kd = glm::vec3(triangles[index_of_object].v[0].color_diffuse[0], triangles[index_of_object].v[0].color_diffuse[1], triangles[index_of_object].v[0].color_diffuse[2]);
        glm::vec3 v2_kd = glm::vec3(triangles[index_of_object].v[1].color_diffuse[0], triangles[index_of_object].v[1].color_diffuse[1], triangles[index_of_object].v[1].color_diffuse[2]);
        glm::vec3 v3_kd = glm::vec3(triangles[index_of_object].v[2].color_diffuse[0], triangles[index_of_object].v[2].color_diffuse[1], triangles[index_of_object].v[2].color_diffuse[2]);
        //the specular color of the three points
        glm::vec3 v1_ks = glm::vec3(triangles[index_of_object].v[0].color_specular[0], triangles[index_of_object].v[0].color_specular[1], triangles[index_of_object].v[0].color_specular[2]);
        glm::vec3 v2_ks = glm::vec3(triangles[index_of_object].v[1].color_specular[0], triangles[index_of_object].v[1].color_specular[1], triangles[index_of_object].v[1].color_specular[2]);
        glm::vec3 v3_ks = glm::vec3(triangles[index_of_object].v[2].color_specular[0], triangles[index_of_object].v[2].color_specular[1], triangles[index_of_object].v[2].color_specular[2]);
        //the shininess of the three points
        double shininess = p.x*triangles[index_of_object].v[0].shininess + p.y*triangles[index_of_object].v[1].shininess + p.z*triangles[index_of_object].v[2].shininess;
        
        glm::vec3 kd = p.x*v1_kd + p.y*v2_kd + p.z*v3_kd;
        glm::vec3 ks = p.x*v1_ks + p.y*v2_ks + p.z*v3_ks;
        
        float dot_r_v = glm::dot(r, v);
        if(dot_r_v < 0) dot_r_v =0;
        float dot_l_n = glm::dot(l,n);
        if(dot_l_n < 0) dot_l_n = 0;
        glm::vec3 object_color = kd * dot_l_n + ks * (float)pow(dot_r_v, shininess);
        
//        glm::vec3 object_color = kd * glm::dot(l,n) + ks * (float)pow(glm::dot(r, v), shininess);
        
        result_color = glm::vec3(light_color.x * object_color.x, light_color.y * object_color.y, light_color.z * object_color.z);
        
        
    }
    return result_color;
}



//end of functions by alexis


//MODIFY THIS FUNCTION
void draw_scene()
{
    directions = returnRayDirections();

    
    for(int i = 0; i < HEIGHT*WIDTH; i++){
        for(int j = 0; j <3; j++){
            colors[i][j] = 255;
        }
        
    }
    
    
    for(int d = 0; d < directions.size(); d++){
        
        double xd = directions[d][0];
        double yd = directions[d][1];
        double zd = directions[d][2];
        bool sphere_intersected = false;
        glm::vec3 current_intersection = glm::vec3(0,0,0);
        bool intersected = false;
        glm::vec3 total_color = glm::vec3(0, 0, 0);
        for(int i = 0; i < num_spheres; i++){
            double sphere_result = intersectSphere(0, 0, 0, xd, yd, zd, i);
            //if intersect with sphere, the return value should be t0:
            if(sphere_result >= 0){
                sphere_intersected = true;
                intersected = true;
                glm::vec3 p = glm::vec3(directions[d][0]*sphere_result, directions[d][1]*sphere_result, directions[d][2]*sphere_result);//intersection position
                current_intersection = p;
//                cout << p.x << ' ' << p.y << ' ' << p.z << endl;
                //cin.get();
                //for each of the light, get the correspondent color on that point:
                //i is the current sphere, j is the current light source
                for(int j = 0; j < num_lights; j++){

                    //if(!isShadow( p.x,  p.y,  p.z, i, j,  true));
                    glm::vec3 current_color = phongModel(p, "sphere", i, j);//intersection position, light direction, shape, object index, light index
                    total_color += current_color;
                }

//                colors[d][0] = 255;
//                colors[d][1] = 255;
//                colors[d][2] = 255;
            }
        }
        //office hour ask about triangle blocked by sphere
        
        
        for(int i = 0; i < num_triangles; i++){
            vector<double> alpha_beta_gamma = intersectTriangle(0, 0, 0, xd, yd, zd, i);
            //if ifntersection occurs alpha_beta_gamma is the intersecting point
            if(alpha_beta_gamma[0]!=0 && alpha_beta_gamma[1]!= 0 && alpha_beta_gamma[2]!= 0){
                glm::vec3 abg = glm::vec3(alpha_beta_gamma[0],alpha_beta_gamma[1],alpha_beta_gamma[2]);
                intersected = true;
                glm::vec3 pos = getTriangleIntersection(i, abg);
                if(!sphere_intersected || glm::length(pos)<=glm::length(current_intersection)){
                    for(int j = 0; j < num_lights; j++){
                        glm::vec3 current_color = phongModel(abg, "triangle", i, j);
                        //                    glm::vec3 current_color = glm::vec3(1,1,1);
                        total_color = glm::vec3(0,0,0);
                        total_color +=current_color;
                    }
                    
                    //                colors[d][0] = 255;
                    //                colors[d][1] = 255;
                    //                colors[d][2] = 255;

                }else{
                    sphere_intersected = false;
                    current_intersection = glm::vec3(0,0,0);
                }
                
            }
        }
        if(!intersected){
            total_color = glm::vec3(1,1,1);
            intersected = false;
        }
        
    
        colors[d][0] = fmin(total_color.x*255, 255);
        colors[d][1] = fmin(total_color.y*255, 255);
        colors[d][2] = fmin(total_color.z*255,255);

    }
    
    
    
    
    
    
    
    //    exit(0);
    
    //a simple test output
    for(unsigned int x=0; x<WIDTH; x++)
    {
        glPointSize(2.0);
        glBegin(GL_POINTS);
        for(unsigned int y=0; y<HEIGHT; y++)
        {
            unsigned int index = x*(HEIGHT)+y;
            plot_pixel(x, y, colors[index][0], colors[index][1], colors[index][2]);
        }
        glEnd();
        glFlush();
    }
    
    printf("Done!\n"); fflush(stdout);
    
    //ray:
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    glColor3f(((double)r) / 255.0f, ((double)g) / 255.0f, ((double)b) / 255.0f);
    glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    buffer[y][x][0] = r;
    buffer[y][x][1] = g;
    buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
    plot_pixel_display(x,y,r,g,b);
    if(mode == MODE_JPEG)
        plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
    printf("Saving JPEG file: %s\n", filename);
    
    ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
    if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
        printf("Error in Saving\n");
    else
        printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
    if(strcasecmp(expected,found))
    {
        printf("Expected '%s ' found '%s '\n", expected, found);
        printf("Parse error, abnormal abortion\n");
        exit(0);
    }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check(check,str);
    fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
    printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check("rad:",str);
    fscanf(file,"%lf",r);
    printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
    char s[100];
    fscanf(file,"%s",s);
    parse_check("shi:",s);
    fscanf(file,"%lf",shi);
    printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
    FILE * file = fopen(argv,"r");
    int number_of_objects;
    char type[50];
    Triangle t;
    Sphere s;
    Light l;
    fscanf(file,"%i", &number_of_objects);
    
    printf("number of objects: %i\n",number_of_objects);
    
    parse_doubles(file,"amb:",ambient_light);
    
    for(int i=0; i<number_of_objects; i++)
    {
        fscanf(file,"%s\n",type);
        printf("%s\n",type);
        if(strcasecmp(type,"triangle")==0)
        {
            printf("found triangle\n");
            for(int j=0;j < 3;j++)
            {
                parse_doubles(file,"pos:",t.v[j].position);
                parse_doubles(file,"nor:",t.v[j].normal);
                parse_doubles(file,"dif:",t.v[j].color_diffuse);
                parse_doubles(file,"spe:",t.v[j].color_specular);
                parse_shi(file,&t.v[j].shininess);
            }
            
            if(num_triangles == MAX_TRIANGLES)
            {
                printf("too many triangles, you should increase MAX_TRIANGLES!\n");
                exit(0);
            }
            triangles[num_triangles++] = t;
        }
        else if(strcasecmp(type,"sphere")==0)
        {
            printf("found sphere\n");
            
            parse_doubles(file,"pos:",s.position);
            parse_rad(file,&s.radius);
            parse_doubles(file,"dif:",s.color_diffuse);
            parse_doubles(file,"spe:",s.color_specular);
            parse_shi(file,&s.shininess);
            
            if(num_spheres == MAX_SPHERES)
            {
                printf("too many spheres, you should increase MAX_SPHERES!\n");
                exit(0);
            }
            spheres[num_spheres++] = s;
        }
        else if(strcasecmp(type,"light")==0)
        {
            printf("found light\n");
            parse_doubles(file,"pos:",l.position);
            parse_doubles(file,"col:",l.color);
            
            if(num_lights == MAX_LIGHTS)
            {
                printf("too many lights, you should increase MAX_LIGHTS!\n");
                exit(0);
            }
            lights[num_lights++] = l;
        }
        else
        {
            printf("unknown type in scene description:\n%s\n",type);
            exit(0);
        }
    }
    return 0;
}

void display()
{
}

void init()
{
    glMatrixMode(GL_PROJECTION);
    glOrtho(0,WIDTH,0,HEIGHT,1,-1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
    //hack to make it only draw once
    static int once=0;
    if(!once)
    {
        draw_scene();
        if(mode == MODE_JPEG)
            save_jpg();
    }
    once=1;
}

void keyboardFunc(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 27: // ESC key
            exit(0); // exit the program
            break;
            
        case 'x':
            // take a screenshot
            //            saveScreenshot("screenshot.jpg");
            break;
    }
}


int main(int argc, char ** argv)
{
    if ((argc < 2) || (argc > 3))
    {
        printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
        exit(0);
    }
    if(argc == 3)
    {
        mode = MODE_JPEG;
        filename = argv[2];
    }
    else if(argc == 2)
        mode = MODE_DISPLAY;
    
    glutInit(&argc,argv);
    loadScene(argv[1]);
    
    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(WIDTH,HEIGHT);
    int window = glutCreateWindow("Chen_Ming_hw3");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboardFunc);
    
    init();
    glutMainLoop();
}

