#pragma once

#include <cmath>

#include "Interaction_base.hpp"

class Interaction_trap: public Interaction_base{
private:
public:
    Vec_3d pos;
    double r_trap;         ///r = (pcl->pos-pos)
    double factor;         ///P(r) = (r.len() > r_trap) * factor * (r.len() - r_trap)^degree
    unsigned int degree;   ///F(r) = (r.len() > r_trap) * factor * degree * (r.len() - r_trap)^(degree-1)

    Interaction_trap(Vec_3d pos, double r_trap, double factor, unsigned int degree);
    double calc_energy (std::vector<Particle *> &gas);
    void   calc_force  (std::vector<Particle *> &gas);
    double calc        (std::vector<Particle *> &gas);
};

class Interaction_uniform_field: public Interaction_base{
private:
public:
    Vec_3d strength;

    Interaction_uniform_field(Vec_3d strength);
    double calc_energy (std::vector<Particle *> &gas);
    void   calc_force  (std::vector<Particle *> &gas);
    double calc        (std::vector<Particle *> &gas);
};
