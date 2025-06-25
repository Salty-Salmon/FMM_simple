#include "../include/Integrator.hpp"

std::vector<Runge_Kutta_4::Pcl_data *> Runge_Kutta_4::init_pdata(std::vector<Particle *> const &gas_){
    std::vector<Runge_Kutta_4::Pcl_data *> gas;
    size_t maintained_amm = 0;
    for (auto pcl : gas_){
        if (!pcl->rigid_body_id){
            ++maintained_amm;
        }
    }
    gas.reserve(maintained_amm);
    for (auto pcl : gas_){
        if (!pcl->rigid_body_id){
            gas.push_back(new Pcl_data(pcl));
        }
    }
    return gas;
}

std::vector<Runge_Kutta_4::Body_data *> Runge_Kutta_4::init_bdata(std::vector<Rigid_body *> const &bodies_){
    std::vector<Runge_Kutta_4::Body_data *> bodies;
    bodies.reserve(bodies_.size());
    for (auto body : bodies_){
        bodies.push_back(new Body_data(body));
    }
    return bodies;
}

Runge_Kutta_4::Runge_Kutta_4(std::vector<Particle*> const &gas_,
                             std::vector<Rigid_body *> const &bodies_,
                             std::vector<Interaction_base *> const &interactions_,
                             double dt_, size_t substep_amm_):
    gas_raw(gas_),
    gas(init_pdata(gas_)),
    bodies_raw(bodies_),
    bodies(init_bdata(bodies_)),
    interactions(interactions_),
    dt(dt_),
    substep_amm(substep_amm_)
{
    calc_derivatives(0, true);

    energy_kin = calc_kinetic_energy(gas_raw);
    energy_tot = energy_kin + energy_pot;
    momentum = calc_momentum(gas_raw);
    angular_mom = calc_angular_mom(gas_raw);
}

Runge_Kutta_4::~Runge_Kutta_4(){
    for (auto pcl : gas){
        delete pcl;
    }
}

void Runge_Kutta_4::calc_derivatives(size_t i, bool last_step){
    for (auto pcl : gas){
        pcl->pos_der[i] = pcl->pcl->vel;
    }
    for (auto body : bodies){
        body->pos_der[i] = body->body->vel;
        body->rotation_der[i] = body->body->calc_rotation_der();
    }

    reset_force(gas_raw);
    if (last_step) {
        energy_pot = 0;
        reset_potential(gas_raw);
    }
    for (auto interaction : interactions){
        if(last_step){
            energy_pot += interaction->calc(gas_raw);
        }else{
            interaction->calc_force(gas_raw);
        }
    }
    force_to_acc(gas_raw);
    for (auto body : bodies){
        body->body->calc_force_acc();
        body->body->calc_torque();
    }

    for (auto pcl : gas){
        pcl->vel_der[i] = pcl->pcl->acc;
    }
    for (auto body : bodies){
        body->vel_der[i] = body->body->acc;
        body->angular_mom_der[i] = body->body->torque;
    }
}
void Runge_Kutta_4::add_derivatives(size_t i, double dt){
    for (auto pcl : gas){
        pcl->pcl->pos += pcl->pos_der[i] * dt;
        pcl->pcl->vel += pcl->vel_der[i] * dt;
    }
    for (auto body : bodies){
        body->body->pos += body->pos_der[i] * dt;
        body->body->vel += body->vel_der[i] * dt;
        body->body->angular_mom += body->angular_mom_der[i] * dt;
        body->body->rotation += body->rotation_der[i] * dt;

        body->body->fill_parts_pos_vel();
    }
}
void Runge_Kutta_4::add_butcher_row(size_t i, double dt){
    for(size_t j=0; j<i; ++j){
        if (butcher_table[i][j] != 0){
            add_derivatives(j, dt * butcher_table[i][j]);
        }
    }
}
void Runge_Kutta_4::integrate_substep(double dt, bool last_step){
    for(size_t i=1; i<STAGE_N; ++i){
        add_butcher_row(i, +dt);
        calc_derivatives(i, false);
        add_butcher_row(i, -dt);
    }
    for(size_t i=0; i<STAGE_N; ++i){
        add_derivatives(i, dt * butcher_summ_k[i]);
    }
    for (auto body : bodies){
        body->body->rotation = body->body->rotation.orthogonalize();
        body->body->fill_parts_pos_vel();
    }
    calc_derivatives(0, last_step);
}

void Runge_Kutta_4::integrate_dt(){
    for(size_t i=0; i<substep_amm-1; ++i){
        integrate_substep(dt/substep_amm, false);
    }
    integrate_substep(dt/substep_amm, true);

    energy_kin = calc_kinetic_energy(gas_raw);
    energy_tot = energy_kin + energy_pot;
    momentum = calc_momentum(gas_raw);
    angular_mom = calc_angular_mom(gas_raw);
}
