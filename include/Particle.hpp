#pragma once

#include "Vec_3d.hpp"

class Particle{
private:

public:
    Vec_3d pos;
    Vec_3d vel;
    Vec_3d acc;
    Vec_3d force;

    double mass;
    double charge;
    double potential;

    bool processed;
    Particle(Vec_3d pos, Vec_3d vel, double mass, double charge):
        pos(pos), vel(vel), acc(),
        mass(mass), charge(charge), potential(0),
        processed(0)
    {

    };

    friend std::ostream& operator <<(std::ostream& os, const Particle & rha){
        return os << "pos:" << rha.pos << "   vel:" << rha.vel << "   acc:" << rha.acc << "   mass:" << rha.mass << "   charge:" << rha.charge << "\n";
    };
    std::ostream& print_minimal(std::ostream& os){
        return os << pos.x << " " << pos.y << " " << pos.z << " " << vel.x << " " << vel.y << " " << vel.z << " " << mass << " " << charge << " " << potential << "\n";
    }
};
