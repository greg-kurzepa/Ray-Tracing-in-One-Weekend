// Z is defined as vertically upwards
// Y forward from camera
// X sideways and to the right

#include <vector>
#include <limits>
#include <time.h>
#include <iostream>
#include <cmath>
#include <vector>

#include <sstream>
#include <string>

#include "gmath.h"
#include "gpng.h"

using namespace gmath;
double min_dist_threshold{0.001}; // minimum distance a point of intersection must be from start of a line to be counted
int inside_count{0};

class Camera {
    public:
        Vec3 up{0,0,1}; // sets the rotation of the field of view box, keeping it viewing "horizontally"
        double viewport_height{2}; // 2 due to legacy reasons; consider an arbitrary magic number. viewport size is always the same. to change zoom, vary focal_length.

        Vec3 viewport_centre;
        Vec3 lookat;
        double aspect_ratio;
        double focal_length;
        Vec3 look_direction;
        double viewport_width;
        Vec3 origin;

        // orthogonal unit vectors to traverse viewport
        Vec3 d_right;
        Vec3 d_up;

        // default setting
        Camera(double aspect_ratio) :
            viewport_centre(Vec3(0,1,0)),
            lookat(Vec3(0,2,0)),
            aspect_ratio(aspect_ratio),
            focal_length(1),
            look_direction((lookat - viewport_centre).unit()),
            viewport_width(viewport_height * aspect_ratio),
            origin(Vec3(0,0,0)),
            d_right(Vec3(1,0,0)),
            d_up(Vec3(0,0,1))
        {}

        // any setting
        Camera(Vec3 viewport_centre, Vec3 lookat, double aspect_ratio, double focal_length) :
            viewport_centre(viewport_centre),
            lookat(lookat),
            aspect_ratio(aspect_ratio),
            focal_length(focal_length),
            look_direction((lookat - viewport_centre).unit()),
            viewport_width(viewport_height * aspect_ratio),
            origin(viewport_centre - look_direction * focal_length),
            d_right(cross(look_direction, up).unit()),
            d_up(-cross(look_direction, d_right).unit())
        {}
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
    metal,
    glass
};

class Hittable {
    public:
        Material material{Material::matte};
        Colour reflectance{0.5, 0.5, 0.5};
        double fuzz{0}; // for metals, should be between 0 and 1
        double refractive_index{1.5}; // for glass, should be >= 1 (1 for air, 1.5 for glass)

        Hittable() {}
        Hittable(Material material) : material(material) { setup(); }
        Hittable(Material material, Colour reflectance, double fuzz) : material(material), reflectance(reflectance), fuzz(fuzz) { setup(); }
        
        void setup() {
            if (material == Material::glass) { reflectance = Colour(1.0, 1.0, 1.0); }
        }

        // returns value of t along input ray that causes intersection with the Hittable
        virtual double intersects(const Line3& ray) const = 0;
        // returns vector normal to surface at specified point
        // virtual UnitVec3 get_normal(const Line3& ray, const double t, bool& intersects_outside) const = 0;
        //
        virtual Line3 get_next_ray(const Line3& ray, const double t) const = 0;
};

// Array of all hittable objects
std::vector<Hittable*> hittables;

class Sphere3 : public Hittable {
    // sphere defined by position of its centre and its radius
    public:
        Vec3 p; // centre
        double r; // radius
        bool is_hollow; // whether the norm_vector should be inverted

        Sphere3() : p(Vec3(0,0,0)), r(0) {}
        Sphere3(Vec3 centre, double radius) : p(Vec3(centre)), r(radius) {}

        // these allow for setting the materials and, if desired, its reflectance
        Sphere3(Vec3 centre, double radius, Material material) : Hittable(material), p(Vec3(centre)), r(radius) {}
        Sphere3(Vec3 centre, double radius, Material material, Colour reflectance, double fuzz, bool is_hollow=false) : Hittable(material, reflectance, fuzz), p(Vec3(centre)), r(radius), is_hollow(is_hollow) {}

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
                double smaller = (-half_b - sqrt(discriminant)) / a;
                if (smaller > min_dist_threshold) { // return the closest value as long as it is in the +ve direction. otherwise, return further value.
                    return smaller;
                } else {
                    return (-half_b + sqrt(discriminant)) / a;
                }
            }
        }

        Line3 get_next_ray(const Line3& ray, const double t) const override {
            // Find normal unit ray reflection vector
            Vec3 normal_unit = (ray(t) - p).unit(); // normal unit vector to sphere pointing out of sphere surface
            Line3 ret_ray;

            switch (material) {
                case Material::matte: {
                    // New ray for next iteration, selected randomly from a unit sphere tangential to the intersected surface
                    // As in https://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
                    Vec3 X{normal_double(), normal_double(), normal_double()};
                    // Line3 next_ray{ray(t), normal_unit + (X * std::pow(random_double(), 1.0/3.0) / X.abs())}; // random point *in* sphere
                    ret_ray = Line3(ray(t), (normal_unit + X / X.abs()));
                    break;
                }
                case Material::metal: {
                    Vec3 X{normal_double(), normal_double(), normal_double()};
                    ret_ray = Line3(ray(t), ray.d - 2*normal_unit*dot(ray.d, normal_unit) + fuzz * X / X.abs());
                    break;
                }
                case Material::glass: {
                    // possibilities: 
                    // entering always -ve to normal
                    // -> Normal sphere: normal_unit correct; refraction ratio 1/1.5
                    // -> Hollow section: normal_unit correct; refraction ratio 1.5
                    // leaving always +ve to normal
                    // -> Normal sphere: normal_unit inverted; refraction ratio 1.5
                    // -> Hollow section: normal_unit inverted; refraction ratio 1/1.5
                    
                    double refraction_ratio;
                    if (dot(normal_unit, ray.d) > 0) { // if ray going from inside to outside...
                        normal_unit = -normal_unit;
                        refraction_ratio = is_hollow ? 1.0/refractive_index : refractive_index;
                    } else { // if ray going from outside to inside...
                        refraction_ratio = is_hollow ? refractive_index : 1.0/refractive_index;
                    }

                    double cos_theta = -dot(normal_unit, ray.d.unit());
                    if (cos_theta < 0) {
                        std::cout << "COS NEGATIVE: " << cos_theta << "\n"; 
                    }

                    if (refraction_ratio * sqrt(1.0 - cos_theta*cos_theta) > 1.0 || schlick_reflectance(cos_theta, refraction_ratio) > random_double()) { // if total internal reflection or schlick reflection...
                        ret_ray = Line3(ray(t), ray.d - 2*normal_unit*dot(ray.d, normal_unit)); // reflect
                    } else {
                        Vec3 refracted_ray_perpendicular = refraction_ratio * (ray.d.unit() + normal_unit * cos_theta);
                        Vec3 refracted_ray_parallel = normal_unit * -sqrt(fabs(1.0 - refracted_ray_perpendicular.abs2()));
                        ret_ray = Line3(ray(t), refracted_ray_perpendicular + refracted_ray_parallel);
                    }
                    break;
                }
                default:
                    std::cerr << "Error in Sphere3.get_next_ray(): No material match found";
                    return Line3(Vec3(0,0,0), Vec3(0,0,0));
            }
            ret_ray.d = ret_ray.d.unit();
            return ret_ray;
        }
    
    private:
        static double schlick_reflectance(double cos_theta, double reflection_ratio) {
            double r0 = (1 - reflection_ratio) / (1 + reflection_ratio);
            r0 = r0*r0;
            return r0 + (1-r0)*pow((1 - cos_theta), 5);
        }
};

/// @brief allows for setting pixel colours using instances of Colour (Vec3) class
class ImageVec : public gpng::Image {
    public:
        ImageVec(double width, double height) : gpng::Image(width, height) {}

        void set_pixel(int column, int row, const Colour& u) {
            int start_idx = row * width * 3 + column * 3;
            image[start_idx] = u.x;
            image[start_idx + 1] = u.y;
            image[start_idx + 2] = u.z;
        }
};

/// @brief 
/// @param n 
/// @return number of collisions before hit background (light source)
Colour ray_recur(int n, Line3& ray, bool do_trace) {
    if (do_trace) {
        std::cout << "Level: " << n << ", Position: " << ray.p << ", Vector: " << ray.d << "\n";
    }

    // Find closest intersection
    double t {std::numeric_limits<double>::infinity()}; // change this to Tmax if you want 0 < t < Tmax instead of 0 < t < infinity
    int smallest_idx {-1};
    for (std::vector<Hittable*>::size_type i = 0; i < hittables.size(); i++) {
        double intersect_t = hittables[i]->intersects(ray);
        if (intersect_t < t && intersect_t > min_dist_threshold) {
            t = intersect_t;
            smallest_idx = i;
        }
    }
    if (do_trace) { std::cout << "smallest_idx: " << smallest_idx << "\n"; }

    if (smallest_idx != -1) { // if at least one object intersects with the ray...
        if (n == 1) { return Colour(0,0,0); } // return black if recursion count limit recur_max reached
        Hittable* closest_item_ptr = hittables[smallest_idx];
        Line3 next_ray = closest_item_ptr->get_next_ray(ray, t);
        return closest_item_ptr->reflectance * ray_recur(n-1, next_ray, do_trace);

    } else { // if hit 'sky' (i.e. if nothing else was hit)...
        // rtow colour scheme
        double t = 0.5*(ray.d.unit().z + 1.0);
        return (1.0-t)*Colour(1.0, 1.0, 1.0) + t*Colour(0.5, 0.7, 1.0);

        // quadrants colour scheme
        // if (ray.d.unit().x >=0) {
        //     if (ray.d.unit().z >= 0) {
        //         return Colour(0.0, 0.0, 0.0); // top-right black
        //     } else {
        //         return Colour(1.0, 0.0, 0.0); // bottom-right red
        //     }
        // } else {
        //     if (ray.d.unit().z >= 0) {
        //         return Colour(0.0, 0.0, 1.0); // top-left blue
        //     } else {
        //         return Colour(0.0, 1.0, 0.0); // bottom-left green
        //     }
        // }
    }
}

int main() {
    int width{1600};
    int height{900};

    ImageVec img(width, height);
    double aspect_ratio = static_cast<double>(width) / static_cast<double>(height);

    // default camera
    // Camera cam(aspect_ratio);

    // Camera(Vec3 viewport_centre, Vec3 lookat, double aspect_ratio, double focal_length)  <- definition
    // camera from above, behind and to the left
    Camera cam(Vec3(-2,-1,2), Vec3(0,1,0), aspect_ratio, 10);

    // 2 spheres scene
    // Sphere3 sphere1(Vec3(0,4,0), 2, Material::glass, Colour(0.7, 0.3, 0.3), 0);
    // Sphere3 sphere2(Vec3(0,4,-102), 100, Material::matte, Colour(0.7, 0.3, 0.3), 0);

    // 3 spheres scene
    // Sphere3 sphere1(Vec3(0,1,0), 0.5, Material::matte, Colour(0.7, 0.3, 0.3), 0);
    // Sphere3 sphere2(Vec3(-1,1,0), 0.5, Material::metal, Colour(0.8, 0.8, 0.8), 0.0);
    // Sphere3 sphere3(Vec3(1,1,0), 0.5, Material::metal, Colour(0.8, 0.6, 0.2), 0.0);
    // Sphere3 sphere4(Vec3(0,1,-100.5), 100, Material::matte, Colour(0.8, 0.8, 0.0), 0);

    // Ray Tracing in One Weekend scene
    Sphere3 sphere1(Vec3(0,1,-100.5), 100, Material::matte, Colour(0.8,0.8,0.0), 0);
    Sphere3 sphere2(Vec3(0,1,0), 0.5, Material::matte, Colour(0.1,0.2,0.5), 0);
    Sphere3 sphere3(Vec3(1,1,0), 0.5, Material::metal, Colour(0.8,0.6,0.2), 0);
    // hollow glass sphere
    Sphere3 sphere4(Vec3(-1,1,0), 0.5, Material::glass, Colour(0.8,0.8,0.8), 0);
    Sphere3 sphere5(Vec3(-1,1,0), 0.4, Material::glass, Colour(0.8,0.8,0.8), 0, true);

    // refraction test scene
    // Sphere3 sphere1(Vec3(0,1,0), 0.5, Material::glass, Colour(0.8,0.8,0.8), 0);
    // Sphere3 sphere2(Vec3(0,1,0), 0.4, Material::glass, Colour(0.8,0.8,0.8), 0, true);

    hittables.push_back(&sphere1);
    hittables.push_back(&sphere2);
    hittables.push_back(&sphere3);
    hittables.push_back(&sphere4);
    hittables.push_back(&sphere5);

    // render!
    std::cout << "start loop";
    for (double y_pixel = 0; y_pixel < img.height; y_pixel++) {
        for (double x_pixel = 0; x_pixel < img.width; x_pixel++) {
            
            // select a ray and print out its coordinates
            double x_trace = 800;
            double y_trace = 708;
            bool do_trace = false;
            if (x_pixel == x_trace && y_pixel == y_trace) {
                do_trace = true;
                std::cout << "Doing trace!\n";
            }

            Colour running_colour{0,0,0};

            // n_antialias rays for antialiasing
            int n_antialias = 4;
            for (int i = 0; i < n_antialias; i++) {
                // Make ray
                Line3 ray;
                ray.p = cam.origin;
                ray.d += cam.look_direction * cam.focal_length;
                ray.d += cam.d_right * cam.viewport_width*((x_pixel + random_double())/img.width - 0.5);
                ray.d += cam.d_up * cam.viewport_height*((y_pixel + random_double())/img.height - 0.5);
                ray.d = ray.d.unit();

                // first argument is maximum recur depth
                running_colour += ray_recur(10, ray, do_trace); // start with refractive_index=1 (air)
            }
            // average
            running_colour /= n_antialias + 1;

            // gamma correction, gamma 2 (colour to the power of 1/2)
            running_colour = pow(running_colour, 0.5);
            // output colour (if tracing this pixel, give it a colour in the image to verify what we are tracing)
            if (do_trace) {
                img.set_pixel(x_pixel, img.height - y_pixel - 1, Colour(0,255,0));
            } else {
                img.set_pixel(x_pixel, img.height - y_pixel - 1, 255*running_colour);
            }
        }
    }
    std::cout << "end loop";

    std::stringstream string_stream;
    string_stream << "images/" << time(NULL) << ".png";
    img.save(string_stream.str());
    // img.save("images/test2.png");

    std::cout << inside_count << "\n";
}