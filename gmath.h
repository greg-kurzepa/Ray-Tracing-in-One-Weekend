#pragma once

#include <cmath>
#include <iostream>
#include <cstdlib>

namespace gmath {

    extern const double pi;

    double random_double();
    double random_double(double min, double max);

    class Vec3 {
        public:
            double x;
            double y;
            double z;

            Vec3();
            Vec3(double x, double y, double z);

            Vec3 operator-() const;

            Vec3& operator+=(const Vec3& v);
            Vec3& operator*=(const double t);
            Vec3& operator/=(const double t);

            double abs() const;
            double abs2() const;
            Vec3 unit() const;
    };

    Vec3 operator+(const Vec3& u, const Vec3& v);
    Vec3 operator-(const Vec3& u, const Vec3& v);

    double dot(const Vec3& u, const Vec3& v); // dot product
    // Vec3 operator^(const Vec3& u, const Vec3& v); // cross product
    Vec3 cross(const Vec3& u, const Vec3& v);

    Vec3 operator*(const Vec3& v, const double t);
    Vec3 operator*(const double t, const Vec3& v);
    Vec3 operator/(const Vec3& v, const double t);

    std::ostream& operator<<(std::ostream& out, const Vec3& v);

    using Colour = Vec3;

}