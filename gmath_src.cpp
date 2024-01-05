#include <cmath>
#include <iostream>
#include <cstdlib>
#include <random>
#include "gmath.h"

namespace gmath {

    // Constants

    const double pi = 3.14159265359;

    // Functions

    /// @brief 
    /// @return random uniform double in range [0,1), i.e. including 0 but not 1
    double random_double() {
        return std::rand() / (RAND_MAX + 1.0);
    }

    /// @brief 
    /// @return random uniform double in range [min,max), i.e. including min but not max
    double random_double(double min, double max) {
        return min + (max-min)*random_double();
    }

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);

    double normal_double() {
        return distribution(generator);
    }

    // Vec3 Class

    Vec3::Vec3() : x(0), y(0), z(0) {}
    Vec3::Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3 Vec3::operator-() const { return Vec3(-x, -y, -z); }

    Vec3& Vec3::operator+=(const Vec3& v) {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }

    Vec3& Vec3::operator-=(const Vec3& v) {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }

    Vec3& Vec3::operator*=(const double t) {
        x *= t;
        y *= t;
        z *= t;
        return *this;
    }

    Vec3& Vec3::operator/=(const double t) {
        return *this *= 1/t;
    }

    double Vec3::abs() const {
        return std::sqrt(abs2());
    }

    double Vec3::abs2() const {
        return x*x + y*y + z*z;
    }

    // check performance of this against just using division...
    Vec3 Vec3::unit() const {
        return *this / abs();
    }

    Vec3 operator+(const Vec3& u, const Vec3& v){
        return Vec3(u.x + v.x, u.y + v.y, u.z + v.z);
    }

    Vec3 operator-(const Vec3& u, const Vec3& v) {
        return Vec3(u.x - v.x, u.y - v.y, u.z - v.z);
    }

    double dot(const Vec3& u, const Vec3& v) {
        return (u.x * v.x + u.y * v.y + u.z * v.z);
    }

    // removed this because of dodgy precedence rules (^ has lower precedence than addition and multiplication and more).
    // better just to have an explicit function, which has obvious precedence.
    // Vec3 operator^(const Vec3& u, const Vec3& v) {
    //     return Vec3(u.y * v.z - u.z * v.y,
    //                 u.z * v.x - u.x * v.z,
    //                 u.x * v.y - u.y * v.x);
    // }

    Vec3 cross(const Vec3& u, const Vec3& v) {
        return Vec3(u.y * v.z - u.z * v.y,
                    u.z * v.x - u.x * v.z,
                    u.x * v.y - u.y * v.x);
    }

    Vec3 pow(const Vec3& u, const double t) {
        return Vec3(std::pow(u.x, t), std::pow(u.y, t), std::pow(u.z, t));
    }

    Vec3 operator*(const Vec3& v, const double t) {
        return Vec3(t * v.x, t * v.y, t * v.z);
    }

    Vec3 operator*(const double t, const Vec3& v) {
        return v * t;
    }

    Vec3 operator*(const Vec3& u, const Vec3& v) { // item-wise multiplication
        return Vec3(u.x * v.x, u.y * v.y, u.z * v.z);
    }

    Vec3 operator/(const Vec3& v, const double t) {
        return v * (1/t);
    }

    std::ostream& operator<<(std::ostream& out, const Vec3& v) {
        return out << v.x << " " << v.y << " " << v.z << "\n";
    }

}