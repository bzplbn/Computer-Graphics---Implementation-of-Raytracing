//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;

int g_width;
int g_height;
int step = 0;
int num_sphere;
int num_light;
vec3 backGround;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

// TODO: add structs for spheres, lights and anything else you may need.
struct Sphere
{
    string name;
    float pos_x;
    float pos_y;
    float pos_z;
    float scl_x;
    float scl_y;
    float scl_z;
    float r;
    float g;
    float b;
    float Ka;
    float Kd;
    float Ks;
    float Kr;
};

struct Light
{
    string name;
    float pos_x;
    float pos_y;
    float pos_z;
    float Ir;
    float Ig;
    float Ib;

};

vector<vec4> g_colors;

float  g_left;
float  g_right;
float  g_top;
float  g_bottom;
float  g_near;
vector<Sphere> sphere;
vector<Light>  light;
string output;
// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    //TODO: add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.
    const int num_labels = 11;// 0    1        2       3         4      5      6         7        8       9          10
    const string labels[] = {"NEAR", "LEFT", "RIGHT", "BOTTOM", "TOP", "RES", "SPHERE", "LIGHT", "BACK", "AMBIENT", "OUTPUT"};
    long label_id = find(labels, labels + num_labels, vs[0]) - labels;
    int i =0 , j = 0;
    
    switch (label_id) {
        case 0:      g_near   = toFloat ( vs[1] );    break;
        case 1:      g_left   = toFloat ( vs[1] );    break;
        case 2:      g_right  = toFloat ( vs[1] );    break;
        case 3:      g_bottom = toFloat ( vs[1] );    break;
        case 4:      g_top    = toFloat ( vs[1] );    break;
        case 6:      sphere[i++]   = {vs[1],toFloat(vs[2]),toFloat(vs[3]),toFloat(vs[4]),toFloat(vs[5]),
                     toFloat(vs[6]),toFloat(vs[7]),toFloat(vs[8]),toFloat(vs[9]),toFloat(vs[10]),toFloat(vs[11]),toFloat(vs[12]),toFloat(vs[13]),toFloat(vs[14])};
                     num_sphere = i;                  break;
        case 7:      light[j++]    = {vs[1],toFloat(vs[2]),toFloat(vs[3]),toFloat(vs[4]),toFloat(vs[5]),
                     toFloat(vs[6]),toFloat(vs[7])};
                     num_light = j;                  break;
        case 8:      backGround = vec3(toFloat(vs[1]), toFloat(vs[2]), toFloat(vs[3]));
        case 9:
        case 10:     output  =  vs[1];                break;

    }
    
    if (vs[0] == "RES")
    {
        g_width = (int)toFloat(vs[1]);
        g_height = (int)toFloat(vs[2]);
        g_colors.resize(g_width * g_height);
    }
    
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine

// TODO: add your ray-sphere intersection routine here.
//[IN] ray
//[IN] sphere
//[OUT] t
//[OUT] intersection
bool intersect(const Ray& ray, const Sphere&  sphere, float t, vec4 intersection)
{
   
    float delta;
    mat4 M, M_inverse;
    vec4 s = ray.origin, s_new;
    vec4 c = ray.dir, c_new;
    float r = 1;
    
    M = Scale(sphere.scl_x, sphere.scl_y, sphere.scl_z) * Translate(sphere.pos_x, sphere.pos_y, sphere.pos_z);
    InvertMatrix(M, M_inverse);
    s_new = M_inverse * s;
    c_new = M_inverse * c;
    
    //caculate delta
    delta = pow(dot(s_new,c_new),2) -pow(length(c_new),2)*(pow(length(s_new),2)-pow(r, 2));
    
    if (delta < 0)
        return false;
    else
    {
        t = -(pow(dot(s_new,c_new),2)/pow(length(c_new),2)) - sqrt(delta/pow(length(c_new),2));
        intersection = s + c * t;
        return true;
    }
    
    
}

// -------------------------------------------------------------------
// Ray tracing
vec3 recursive_trace(const Ray& ray)
{
    vec3 color;
    Ray ray_rf;
    vec4 normal;
    vec3 local , reflected;
    vec4 closest_point;
    vec4 point;
    float t;
    float closest_t = 1000;
    bool intersection;
    int maxStep = 5;
    
    step++;
    
    if(step == maxStep)
    {
        return color;
    }
    
    for(int i = 0; i < num_sphere; i++)
    {
        if((intersection = intersect(ray, sphere[i], t, point)) == true)
        {
            if (t < closest_t) {
                closest_t = t;
            }
        }
    }
    closest_point = ray.origin + closest_t * ray.dir;
    normal = normalize(vec4(closest_point[0], closest_point[1], closest_point[2], 0));
    ray_rf.origin = closest_point;
    ray_rf.dir = -2 * dot(normal,ray.dir) * normal + ray.dir;
    
    color = recursive_trace(ray_rf);
    return color;
    
}

vec4 trace(const Ray& ray)
{
    // TODO: implement your ray tracing routine here.
    
//    raytrace( ray )
//    P = closest intersection
//    color_local = ShadowRay(light1, P) + ...+ ShadowRay(lightN, P)
//    color_reflect = raytrace(reflected_ray )
//    color_refract = raytrace(refracted_ray )
//    color = color_local +
//    + krfl* color_reflect + krfa* color_refract
//    return( color )
    

//    function trace(p,d,step)
//    {
//        color local, reflected;
//        point q;
//        normal n;
//        
//        if(step > max){
//            return(backgroundColor);
//        }
//        q = intersect(p,d,status);
//        if(status == light_source){
//            return(backgroundColor);
//        }
//        
//        n = normal(q);
//        r = reflect(q,n);
//
//        
//        local = phong(q, n, r);
//        reflected = trace(q, r, step+1);
//
//        
//        return(local + reflected);
//    }
    
    Ray ray_rf;
    //Ray newRay;
    vec4 normal;
    vec4 color;
    vec3 local , reflected;
    vec4 closest_point;
    vec4 point;
    vec3 temp;
    float t;
    float closest_t = 1000;
    bool intersection;
    
    //find closest intersection
    for(int i = 0; i < num_sphere; i++)
    {
        if((intersection = intersect(ray, sphere[i], t, point)) == true)
        {
            if (t < closest_t) {
                closest_t = t;
            }
        }
    }
    closest_point = ray.origin + closest_t * ray.dir;
    
    //caculate shadow
    normal = normalize(vec4(closest_point[0], closest_point[1], closest_point[2], 0));
    ray_rf.origin = closest_point;
    ray_rf.dir = -2 * dot(normal,ray.dir) * normal + ray.dir;
    
    for(int j = 0; j < num_light; j++)
    {
        point = vec4(light[j].pos_x, light[j].pos_y , light[j].pos_z, 1);
        if ((point - ray_rf.origin) != ray_rf.dir) {
            local += backGround;
        }
        else
            local += vec3(light[j].Ir,light[j].Ig, light[j].Ib);
    }
    
    reflected = recursive_trace(ray_rf);
    temp = local + reflected;
    color = vec4(temp[0], temp[1], temp[2], 1);
    return color;
}

vec4 getDir(int ix, int iy)  //Turn screen space pixels into world space rays
{
    // TODO: modify this. This should return the direction from the origin
    // to pixel (ix, iy), normalized.
    
    vec4 dir;
    float x, y, z;
    //x = left + right * (2 * ix) / RES
    x =  g_left + ix * g_right * 2 / g_width;
    //iy = bottom + top * (2 * iy) / RES
    y =  g_bottom + iy * g_top * 2 / g_height;
    z =  - g_near;
    //dir = vec4(0.0f, 0.0f, -1.0f, 0.0f);
    dir = vec4(x, y, z, 1) - vec4(0, 0, 0, 1);
    return normalize(dir);
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, const char* fname, unsigned char* pixels)
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // TODO: clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++){
                if (((float*)g_colors[y*g_width+x])[i] > 1)
                {
                    ((float*)g_colors[y*g_width+x])[i] = 1;
                }
                buf[y*g_width*3+x*3+i] = (unsigned char)(((float*)g_colors[y*g_width+x])[i] * 255.9f);
                
            }
    // TODO: change file name based on input file name.
    const char* fname = output.c_str();
    savePPM(g_width, g_height, fname, buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 1)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
//    loadFile("testBackground.txt");
    render();
    saveFile();
	return 0;
}

