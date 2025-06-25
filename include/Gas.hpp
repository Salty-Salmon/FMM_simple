#pragma once

#include <iostream>
#include <fstream>
#include <vector>

#include "Particle.hpp"
#include "Rigid_body.hpp"
#include "Interaction_base.hpp"

std::pair<Vec_3d, Vec_3d> get_bounding_box  (std::vector<Particle *> const &gas);
std::pair<Vec_3d, double> get_bounding_cube (std::vector<Particle *> const &gas);

void reset_force     (std::vector<Particle *> const &gas);
void reset_potential (std::vector<Particle *> const &gas);
void force_to_acc    (std::vector<Particle *> const &gas);

double calc_kinetic_energy   (std::vector<Particle *> const &gas);
Vec_3d calc_momentum         (std::vector<Particle *> const &gas);
Vec_3d calc_angular_mom      (std::vector<Particle *> const &gas);
double calc_potential_energy (std::vector<Particle *> const &gas,
                              std::vector<Interaction_base *> &inter_arr);

double init_redundant(std::vector<Particle *> const &gas,
                      std::vector<Rigid_body *> &bodies,
                      std::vector<Interaction_base *> &inter_arr);
double euler(std::vector<Particle *> const &gas,
             std::vector<Rigid_body *> &bodies,
             std::vector<Interaction_base *> &inter_arr,
             double dt, bool calc_ener);
double verlet (std::vector<Particle *> const &gas,
               std::vector<Rigid_body *> &bodies,
               std::vector<Interaction_base *> &inter_arr,
               double dt, bool calc_ener);



