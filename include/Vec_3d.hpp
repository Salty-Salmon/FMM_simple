#pragma once

#include <cmath>
#include <iostream>

#include <random>
#include <chrono>

class Vec_3d{
private:

public:
    double x, y, z;
//    double &operator[](size_t);
//
//    Vec_3d(double, double, double);
//    Vec_3d();
//
//    Vec_3d(Vec_3d const &) = default;
//    ~Vec_3d() = default;
//    Vec_3d& operator=(Vec_3d const &) = default;
//
//    Vec_3d operator+(Vec_3d const &) const;
//    Vec_3d operator-() const;
//    Vec_3d operator-(Vec_3d const &) const;
//    Vec_3d& operator+=(Vec_3d const &);
//    Vec_3d& operator-=(Vec_3d const &);
//
//    Vec_3d operator*(double const &) const;
//    friend Vec_3d operator*(double const &, Vec_3d const &);
//    Vec_3d operator/(double const &) const;
//    Vec_3d& operator*=(double const &);
//    Vec_3d& operator/=(double const &);
//
//    double operator*(Vec_3d const &) const;
//    double sqr() const;
//    double len() const;
//
//    friend std::ostream& operator<<(std::ostream& os, const Vec_3d& rha);

    double& operator[](size_t ind){
        if (ind == 0) {return x;}
        if (ind == 1) {return y;}
        if (ind == 2) {return z;}
        return x;
    }
    Vec_3d(double x, double y, double z):x(x), y(y), z(z){}
    Vec_3d():Vec_3d(0, 0, 0){}
    Vec_3d operator+(Vec_3d const &rha) const{
        return Vec_3d(x + rha.x, y + rha.y, z + rha.z);
    }
    Vec_3d operator-() const{
        return Vec_3d(-x, -y, -z);
    }
    Vec_3d operator-(Vec_3d const &rha) const{
        return *this + (-rha);
    }
    Vec_3d& operator+=(Vec_3d const &rha){
        *this = *this + rha;
        return *this;
    }
    Vec_3d& operator-=(Vec_3d const &rha){
        *this = *this - rha;
        return *this;
    }

    Vec_3d operator*(double const &k) const{
        return Vec_3d(k * this->x, k * this->y, k * this->z);
    }
    friend Vec_3d operator*(double const &k, Vec_3d const &rha){
        return rha * k;
    }
    Vec_3d operator/(double const &k) const{
        return *this * (1/k);
    }
    Vec_3d& operator*=(double const &k){
        *this = *this * k;
        return *this;
    }
    Vec_3d& operator/=(double const &k){
        *this = *this / k;
        return *this;
    }

    double operator*(Vec_3d const &rha) const{
        return x*rha.x + y*rha.y + z*rha.z;
    }
    double sqr() const{
        return (*this) * (*this);
    }
    double len() const{
        return std::sqrt(this->sqr());
    }

    friend std::ostream& operator<<(std::ostream &os, const Vec_3d &rha) {
        os << "(" << rha.x << ", " << rha.y << ", " << rha.z << ")";
        return os;
    }
};

inline int sign(int x){
    return (x > 0) - (x < 0);
}

inline double sqr(double x){
    return x*x;
}

inline double cube(double x){
    return x*x*x;
}

Vec_3d rotate_a_to_b(Vec_3d a, Vec_3d b, Vec_3d p);

double rand_uns(double min, double max);
Vec_3d rand_unit_vec();
Vec_3d rand_unit_ball_vec();
Vec_3d rand_unit_segment(Vec_3d axis, double theta_max);
