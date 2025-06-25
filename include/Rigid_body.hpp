#pragma once

#include <vector>
#include <iostream>
#include <iomanip>

#include "Particle.hpp"

class Rotation_matr{
private:
    double matr[3][3];

public:
    static Rotation_matr const identity;
    static double eps;

    Rotation_matr();
    Rotation_matr(double const (&matr_)[3][3]);
    Rotation_matr(Rotation_matr const &matr_);

    double* const operator[](size_t const &ind);

    Rotation_matr operator*(Rotation_matr const &rha) const;
    Rotation_matr& operator*=(Rotation_matr const &rha);
    Rotation_matr operator*(double const &rha) const;
    friend Rotation_matr operator*(double const &lha, Rotation_matr const &rha);
    Rotation_matr& operator*=(double const &rha);
    Rotation_matr operator/(double const &rha) const;
    Rotation_matr& operator/=(double const &rha);
    Vec_3d operator*(Vec_3d const &rha) const;
    friend Vec_3d operator*(Vec_3d const &lha, Rotation_matr const &rha);

    Rotation_matr operator+() const;
    Rotation_matr operator+(Rotation_matr const &rha) const;
    Rotation_matr& operator+=(Rotation_matr const &rha);
    Rotation_matr operator-() const;
    Rotation_matr operator-(Rotation_matr const &rha) const;
    Rotation_matr& operator-=(Rotation_matr const &rha);

    double norm() const;
    Rotation_matr transpose() const;
    static Rotation_matr diag(Vec_3d const &vec);
    Rotation_matr orthogonalize() const;
    std::pair<Rotation_matr, Vec_3d> diagonalize_q_form() const;

    static Rotation_matr delta_ij(size_t i, size_t j);
    static Rotation_matr rot_to_matr(Vec_3d const &lha);

    friend std::ostream& operator<<(std::ostream &os, const Rotation_matr &rha);
};



class Rigid_body{
private:
    static size_t rigid_body_id;
    struct Part{
        Particle * const pcl;
        Vec_3d const rel_pos;

        Part(Particle * const pcl_, Vec_3d const rel_pos_):
            pcl(pcl_),
            rel_pos(rel_pos_)
        {

        };
    };
public:
    Vec_3d pos;
    Vec_3d vel;
    Vec_3d acc;
    Vec_3d force;

    Rotation_matr rotation;
    Vec_3d angular_mom;
    Vec_3d angular_vel;
    Vec_3d torque;

    const double mass;
    const Vec_3d main_inertia;

    const std::vector<Part> parts;

    const size_t id;

    const Vec_3d main_inertia_inverse;

private:
    double calc_mass_fill_pos_vel_angmom(std::vector<Particle *> const &gas);
    Vec_3d calc_main_inertia_fill_rotation(std::vector<Particle *> const &gas);
    std::vector<Part> zip_parts(std::vector<Particle *> const &gas);
public:
    Rigid_body(std::vector<Particle *> const &gas);

    Rotation_matr calc_rotation_der();
    void fill_parts_pos_vel();
    void calc_force_acc();
    void calc_torque();
};
