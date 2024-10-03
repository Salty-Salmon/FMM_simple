#pragma once

#include <cmath>
#include <iomanip>
#include <vector>

#include "Interaction_base.hpp"
#include "Particle.hpp"
#include "Real_harmonics.hpp"
#include "Gas.hpp"

class FMM_octtree_node{
private:

public:
    Vec_3d pos;
    double r;

    std::vector<Particle *> gas;
    std::vector<Particle *> gas_neighbour;
    std::vector<Particle *> gas_interaction;

    FMM_octtree_node(Vec_3d pos, double r, std::vector<Particle *> &gas);
};

class Interaction_Coulomb_FMM_nlogn: public Interaction_base{
private:
    Multipole_calculator mc;
    Rel_dist_4force rel_dist_4force;
    const double mp_rel_dist;
    const double mp_dir_time_ratio;

    double gas_energy;
    bool must_calc_energy;
    bool must_calc_force;

    void calc_direct    (std::vector<Particle *> const &inner, std::vector<Particle *> const &outer);
    void calc_direct    (std::vector<Particle *> const &gas);
    void calc_multipole (std::vector<Particle *> const &inner, std::vector<Particle *> const &outer, Vec_3d origin);
    void calc_node_handler (FMM_octtree_node *curr);

    void subdivide_gas (FMM_octtree_node *curr, std::vector<Particle *> ans[2][2][2]);
    void complete_gas  (FMM_octtree_node *parent, FMM_octtree_node *curr);

    void calc_handler_recursion (FMM_octtree_node *curr);
    void calc_handler (std::vector<Particle *> &gas);

    void print_tree_recursion(FMM_octtree_node *curr, std::ostream& os);
public:
    double coulomb_constant;

    Interaction_Coulomb_FMM_nlogn (double coulomb_constant_, double max_err_force_, double mp_rel_dist_, double mp_dir_time_ratio_);
    double calc_energy (std::vector<Particle *> &gas);
    void   calc_force  (std::vector<Particle *> &gas);
    double calc        (std::vector<Particle *> &gas);

    void print_tree (std::vector<Particle *> &gas, std::ostream& os);
};

