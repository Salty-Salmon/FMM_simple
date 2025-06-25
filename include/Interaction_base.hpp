#pragma once

#include <vector>
#include <cmath>

#include "Particle.hpp"

class Interaction_base{
private:

public:
    virtual double calc_energy (std::vector<Particle *> const &gas) = 0;
    virtual void   calc_force  (std::vector<Particle *> const &gas) = 0;
    virtual double calc        (std::vector<Particle *> const &gas) = 0;

    virtual ~Interaction_base() = default;
};

class Interaction_empty: public Interaction_base{
private:

public:
    double calc_energy (std::vector<Particle *> const &gas){
        double energy = 0;
        return energy;
    };
    void   calc_force (std::vector<Particle *> const &gas){

    };
    double calc (std::vector<Particle *> const &gas){
        double energy = 0;
        return energy;
    };
};

class Interaction_trap: public Interaction_base{
private:
public:
    Vec_3d pos;
    double r_trap;         ///r = (pcl->pos-pos)
    double factor;         ///P(r) = (r.len() > r_trap) * factor * (r.len() - r_trap)^degree
    unsigned int degree;   ///F(r) = (r.len() > r_trap) * factor * degree * (r.len() - r_trap)^(degree-1)

    Interaction_trap(Vec_3d pos, double r_trap, double factor, unsigned int degree);
    double calc_energy (std::vector<Particle *> const &gas);
    void   calc_force  (std::vector<Particle *> const &gas);
    double calc        (std::vector<Particle *> const &gas);
};

class Interaction_uniform_field: public Interaction_base{
private:
public:
    Vec_3d strength;

    Interaction_uniform_field(Vec_3d strength);
    double calc_energy (std::vector<Particle *> const &gas);
    void   calc_force  (std::vector<Particle *> const &gas);
    double calc        (std::vector<Particle *> const &gas);
};
