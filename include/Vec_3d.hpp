#pragma once

#include <cmath>
#include <iostream>
#include <iomanip>

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
        static double trash = 0;
        trash = 0;
        return trash;
    }
    double operator[](size_t ind) const{
        if (ind == 0) {return x;}
        if (ind == 1) {return y;}
        if (ind == 2) {return z;}
        return 0;
    }

    Vec_3d(double x, double y, double z):x(x), y(y), z(z){}
    Vec_3d():Vec_3d(0, 0, 0){}

    Vec_3d operator+() const{
        return Vec_3d(x, y, z);
    }
    Vec_3d operator+(Vec_3d const &rha) const{
        return Vec_3d(x + rha.x, y + rha.y, z + rha.z);
    }
    Vec_3d& operator+=(Vec_3d const &rha){
        *this = *this + rha;
        return *this;
    }

    Vec_3d operator-() const{
        return Vec_3d(-x, -y, -z);
    }
    Vec_3d operator-(Vec_3d const &rha) const{
        return *this + (-rha);
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
    Vec_3d cross(Vec_3d const &rha) const{
        return Vec_3d(this->y * rha.z - this->z * rha.y,
                      this->z * rha.x - this->x * rha.z,
                      this->x * rha.y - this->y * rha.x);
    }
    static Vec_3d cross(Vec_3d const &lha, Vec_3d const &rha){
        return lha.cross(rha);
    }
    double sqr() const{
        return (*this) * (*this);
    }
    double len() const{
        return std::sqrt(this->sqr());
    }

    static size_t out_precision;
    friend std::ostream& operator<<(std::ostream &os, const Vec_3d &rha) {
        os << std::fixed << std::setprecision(out_precision);
        os << "(";
        for (size_t i=0; i<3; ++i){
            os << std::setw(3+out_precision) << rha[i];
            if (i!=2) { os << ", "; }
        }
        os << ")";
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


inline unsigned long long nanoseconds(){ //takes ~(0.9 to 1.5) mcs itself
    static auto t_0 = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - t_0).count();
}

Vec_3d rotate_a_to_b(Vec_3d a, Vec_3d b, Vec_3d p);

double rand_uns(double min, double max);
Vec_3d rand_unit_vec();
Vec_3d rand_unit_ball_vec();
Vec_3d rand_unit_segment(Vec_3d axis, double theta_max);
