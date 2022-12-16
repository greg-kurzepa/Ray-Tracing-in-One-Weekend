// Z is defined as vertically upwards
// Y forward from camera
// X sideways and to the right

#include <vector>
#include <limits>
#include <time.h>

#include <sstream>
#include <string>

#include "gmath.h"
#include "gpng.h"

using namespace gmath;
// UnitVec3 has no extra functionality on top of Vec3, it is simply a semantic way of keeping track of vectors that are *meant* to be unit vectors.
using UnitVec3 = Vec3;

int inside_count{0};

// can currently put it wherever you want but cannot change angle or position of viewport relative to camera origin
class Camera {
    public:
        Vec3 origin;
        Vec3 centre;
        double aspect_ratio;
        double viewport_height;
        double viewport_width;
        const double viewport_dist = 1; // viewport a distance of 1 unit in y direction from camera origin

        Camera() : origin(Vec3(0,0,0)), aspect_ratio(16.0 / 9.0), viewport_height(2) {
            viewport_width = viewport_height * aspect_ratio;
            centre = origin + Vec3(0, viewport_dist, 0);
        }
        Camera(Vec3 origin, double aspect_ratio, double height) : origin(origin), aspect_ratio(aspect_ratio), viewport_height(height) {
            viewport_width = viewport_height * aspect_ratio;
            centre = origin + Vec3(0, viewport_dist, 0);
        }
};

class Line3 {
    public:
        // line defined by a point it intesects and a direction
        Vec3 p; // point
        Vec3 d; // direction

        Line3() : p(Vec3(0,0,0)), d(Vec3(0,0,0)) {}
        Line3(Vec3 point, Vec3 direction) : p(point), d(direction) {}

        Vec3 operator()(const double t) const {  // get position vector at a point t along line
            return p + t*d;
        }
};

enum class Material {
    matte,
};

class Hittable {
    public:
        // returns value of t along input ray that causes intersection with the Hittable
        virtual double intersects(const Line3& ray) const = 0;
        // returns vector normal to surface at specified point
        virtual UnitVec3 get_normal(const Line3& ray, const double t, bool& intersects_outside) const = 0;
};

// Array of all hittable objects
std::vector<Hittable*> hittables;

class Sphere3 : public Hittable {
    // sphere defined by position of its centre and its radius
    public:
        Vec3 p; // centre
        double r; // radius
        Material material = Material::matte;

        Sphere3() : p(Vec3(0,0,0)), r(0) {}
        Sphere3(Vec3 centre, double radius) : p(Vec3(centre)), r(radius) {}

        // quadratic equation for intersection has at least one solution (i.e. ray hits sphere) if discriminant >= 0
        double intersects(const Line3& ray) const override {
            // simplified form of quadratic
            Vec3 oc = ray.p - p;
            double a = ray.d.abs2();
            double half_b = dot(oc, ray.d);
            double c = oc.abs2() - r*r;
            
            double discriminant = half_b * half_b - a*c;
            
            if (discriminant < 0) {
                return -1.0;
            } else {
                return (-half_b - sqrt(discriminant)) / a; // - not + because disc is +ve and we want the closest point of intersection (i.e. smallest t along line)
            }
        }

        // returns a vector facing outwards from sphere
        UnitVec3 get_normal(const Line3& ray, const double t, bool& intersects_outside) const override {
            UnitVec3 normalvec = (ray(t) - p).unit();
            if (dot(normalvec, ray.d) > 0) { // if intersection is on inside of sphere instead of outside...
                intersects_outside = false;
                normalvec = -normalvec;
                // std::cout << "INSIDE!";
                inside_count += 1;
                // std::cout << dot(normalvec, intersect_ray) / (normalvec.abs() * intersect_ray.abs());
                if (inside_count < 0) { std::cout << "OVERFLOW\n"; }
            }
            return normalvec;
        }
};

/// @brief allows for setting pixel colours using instances of Colour (Vec3) class
class ImageVec : public gpng::Image {
    public:
        ImageVec(double width, double height) : gpng::Image(width, height) {}

        void set_pixel(int column, int row, const Vec3& u) {
            int start_idx = row * width * 3 + column * 3;
            image[start_idx] = u.x;
            image[start_idx + 1] = u.y;
            image[start_idx + 2] = u.z;
        }
};

/// @brief 
/// @param n 
/// @return number of collisions before hit background (light source)
Colour ray_recur(int n, int recur_max, Line3& ray) {

    // Find closest intersection
    double t {std::numeric_limits<double>::infinity()}; // change this to Tmax if you want 0 < t < Tmax instead of 0 < t < infinity
    int smallest_idx {-1};
    for (int i = 0; i < hittables.size(); i++) {
        double intersect_t = (*hittables[i]).intersects(ray);
        if (intersect_t < t && intersect_t > 0) {
            t = intersect_t;
            smallest_idx = i;
        }
    }
    Hittable* closest_item = hittables[smallest_idx];

    if (smallest_idx != -1) { // if at least one object intersects with the ray...
        if (n == recur_max - 1) { return Colour(0,0,0); } // return black if recursion count limit recur_max reached

        // Find normal unit ray reflection vector & whether the ray hit the object on its outside or inside
        bool intersects_outside {true};
        UnitVec3 normal_unit = (*closest_item).get_normal(ray, t, intersects_outside);

        if (!intersects_outside) {
            return Colour(0,0,255);
        }

        // New ray for next iteration, selected randomly from a unit sphere tangential to the intersected surface
        // As in https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
        Vec3 X{normal_double(), normal_double(), normal_double()};
        // std::cout << (X * std::pow(random_double(), 1.0/3.0) / X.abs());
        Line3 next_ray{ray(t), normal_unit + (X * std::pow(random_double(), 1.0/3.0) / X.abs())};

        return 0.5 * ray_recur(n+1, recur_max, next_ray);
    }
    // 255*((1.0-z_pixel/img.height)*Colour(1.0, 1.0, 1.0) + (z_pixel/img.height)*Colour(0.5, 0.7, 1.0))
    else { return Colour(255,255,255); }
}

int main() {
    Camera cam;
    ImageVec img(1920, 1080);

    Sphere3 sphere1(Vec3(0,10,0), 2);
    Sphere3 sphere2(Vec3(-4,10,0), 2);
    Sphere3 sphere3(Vec3(4 ,10,0), 2);
    Sphere3 sphere4(Vec3(0,10,-102), 100);

    // Sphere3 sphere1(Vec3(0,1,0), 0.5);
    // Sphere3 sphere2(Vec3(0,1,-100.5), 100);
    hittables.push_back(&sphere1);
    hittables.push_back(&sphere2);
    hittables.push_back(&sphere3);
    hittables.push_back(&sphere4);

    // render!
    for (double z_pixel = 0; z_pixel < img.height; z_pixel++) {
        for (double x_pixel = 0; x_pixel < img.width; x_pixel++) {

            Colour running_colour{0,0,0};
            // n_antialias rays for antialiasing
            int n_antialias = 4;
            for (int i = 0; i < n_antialias; i++) {
                // Make ray
                Line3 ray;
                ray.p = cam.origin;
                ray.d.x = cam.viewport_width*((x_pixel + random_double())/img.width - 0.5);
                ray.d.y = cam.viewport_dist;
                ray.d.z = cam.viewport_height*((z_pixel + random_double())/img.height - 0.5);

                // number of bounces of ray before it hit the background
                int recur_max = 5;
                running_colour += ray_recur(0, recur_max, ray);
                if (running_colour.z == 255.0/2.0 && running_colour.abs() == 255.0/2.0) {
                    std::cout << x_pixel << " " << z_pixel << "\n";
                }
            }

            img.set_pixel(x_pixel, img.height - z_pixel - 1, running_colour / (n_antialias + 1));
            // std::cout << inside_count << "\n";
            // inside_count = 0;
        }
    }

    std::stringstream string_stream;
    string_stream << "images/" << time(NULL) << ".png";
    img.save(string_stream.str());
    // img.save("images/test2.png");

    std::cout << inside_count << "\n";
}