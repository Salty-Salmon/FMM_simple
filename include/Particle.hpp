#pragma once

#include "Vec_3d.hpp"

class Particle{
private:

public:
    Vec_3d pos;
    Vec_3d vel;
    Vec_3d acc;
    Vec_3d force;

    double const mass;
    double const charge;
    double potential;

    size_t rigid_body_id;
    bool processed;
    Particle(Vec_3d pos, Vec_3d vel, double mass, double charge):
        pos(pos), vel(vel), acc(),
        mass(mass), charge(charge), potential(0),
        rigid_body_id(0),
        processed(false)
    {

    };

private:
    enum class Print_mode {print_verbose, print_minimal};
    static Print_mode print_mode;
public:
    static void set_print_verbose() {print_mode = Print_mode::print_verbose;};
    static void set_print_minimal() {print_mode = Print_mode::print_minimal;};
    friend std::ostream& operator <<(std::ostream& os, const Particle & rha){
        os << std::fixed << std::setprecision(Vec_3d::out_precision);
        if (print_mode == Print_mode::print_verbose){
            return os << "pos:"      << rha.pos       << "   "
                      << "vel:"      << rha.vel       << "   "
                      << "acc:"      << rha.acc       << "   "
                      << "mass:"     << rha.mass      << "   "
                      << "charge:"   << rha.charge    << "   "
                      << "pot_ener:" << rha.potential << "   "
                      << "body_id:"  << rha.rigid_body_id << "   "
                      << "\n";
        }else{
            size_t w = 3 + Vec_3d::out_precision;
            return os << std::setw(w) << rha.pos.x << " " << std::setw(w) << rha.pos.y << " " << std::setw(w) << rha.pos.z << " "
                      << std::setw(w) << rha.vel.x << " " << std::setw(w) << rha.vel.y << " " << std::setw(w) << rha.vel.z << " "
                      << std::setw(w) << rha.mass  << " " << std::setw(w) << rha.charge << " "
                      << std::setw(w) << rha.potential << " " << rha.rigid_body_id << "\n";
        }
    };
};
