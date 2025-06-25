#pragma once

#include <iostream>
#include <vector>

#include "Gas.hpp"
#include "Particle.hpp"
#include "Rigid_body.hpp"
#include "Interaction_base.hpp"


class Runge_Kutta_4{
public:
    static constexpr size_t calc_force_per_step = 4;
    static constexpr size_t STAGE_N = 4;
private:
    static constexpr double butcher_table[STAGE_N][STAGE_N] = {
        {  0,   0,   0,   0},
        {0.5,   0,   0,   0},
        {  0, 0.5,   0,   0},
        {  0,   0, 1.0,   0}
    };
    static constexpr double butcher_summ_k[STAGE_N] =
    {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
    //static constexpr double butcher_time_k[STAGE_N] =
    //{  0, 0.5, 0.5, 1.0};
    struct Pcl_data{
        Particle * const pcl;
        Vec_3d pos_der[STAGE_N];
        Vec_3d vel_der[STAGE_N];
        Pcl_data (Particle *pcl_): pcl(pcl_){
            for (size_t i=0; i<STAGE_N; ++i){
                pos_der[i] = Vec_3d();
                vel_der[i] = Vec_3d();
            }
        }
    };
    struct Body_data{
        Rigid_body * const body;
        Vec_3d pos_der[STAGE_N];
        Vec_3d vel_der[STAGE_N];
        Vec_3d angular_mom_der[STAGE_N];
        Rotation_matr rotation_der[STAGE_N];

        Body_data (Rigid_body *body_): body(body_){
            for (size_t i=0; i<STAGE_N; ++i){
                pos_der[i] = Vec_3d();
                vel_der[i] = Vec_3d();
                angular_mom_der[i] = Vec_3d();
                rotation_der[i] = Rotation_matr();
            }
        }
    };

    double energy_kin;
    double energy_pot;
    double energy_tot;
    Vec_3d momentum;
    Vec_3d angular_mom;

    std::vector<Pcl_data *> init_pdata(std::vector<Particle *> const &gas_);
    std::vector<Body_data *> init_bdata(std::vector<Rigid_body *> const &bodies_);

    void calc_derivatives(size_t i, bool last_step);
    void add_derivatives(size_t i, double dt);
    void add_butcher_row(size_t i, double dt);
    void integrate_substep(double dt, bool last_step);
public:
    const std::vector<Particle *> gas_raw;
    const std::vector<Pcl_data *> gas;
    const std::vector<Rigid_body *> bodies_raw;
    const std::vector<Body_data *> bodies;
    const std::vector<Interaction_base *> interactions;
    double dt;
    size_t substep_amm;

    Runge_Kutta_4(std::vector<Particle *> const &gas_,
                  std::vector<Rigid_body *> const &bodies_,
                  std::vector<Interaction_base *> const &interactions,
                  double dt_, size_t substep_amm_);
    ~Runge_Kutta_4();

    void integrate_dt();

    double get_energy_kin() {return energy_kin;};
    double get_energy_pot() {return energy_pot;};
    double get_energy_tot() {return energy_tot;};
    Vec_3d get_momentum() {return momentum;};
    Vec_3d get_angular_mom() {return angular_mom;};
};
