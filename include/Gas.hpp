#pragma once

#include <iostream>
#include <fstream>
#include <vector>

#include "Particle.hpp"
#include "Interaction_base.hpp"
#include "Stopwatch.hpp"

std::pair<Vec_3d, Vec_3d> get_bounding_box  (std::vector<Particle *> &gas);
std::pair<Vec_3d, double> get_bounding_cube (std::vector<Particle *> &gas);
Vec_3d get_vel_cm (std::vector<Particle *> &gas);

void reset_force     (std::vector<Particle *> &gas);
void reset_potential (std::vector<Particle *> &gas);
void force_to_acc    (std::vector<Particle *> &gas);

void read_gas    (std::vector<Particle *> &gas, const std::string &name);
void print_gas   (std::vector<Particle *> &gas, std::ostream& os);
void print_frame (std::vector<Particle *> &gas, int i);

double calc_kinetic_energy   (std::vector<Particle *> &gas);
double calc_potential_energy (std::vector<Particle *> &gas, std::vector<Interaction_base *> &inter_arr);
void verlet (std::vector<Particle *> &gas, std::vector<Interaction_base *> &inter_arr, double dt);

class Verlet_adaptive{
private:
    struct Particle_img{
        Particle * const pcl;
        Vec_3d pos_backup;
        Vec_3d vel_backup;
        Vec_3d acc_backup;

        Particle_img(Particle * pcl_): pcl(pcl_){

        };
    };
public:
    std::vector<Particle *> gas;
    std::vector<Particle_img *> gas_img;

    std::vector<Interaction_base *> inter_arr;
    double dt;
    double ener_err_max;

    double ener_curr;
    double subdivision;
    size_t subframe_amm;
    size_t itr_counter;

    Verlet_adaptive(std::vector<Particle *> &gas_, std::vector<Interaction_base *> &inter_arr_, double dt_, double ener_err_max_):
        gas(gas_), inter_arr(inter_arr_),
        dt(dt_), ener_err_max(ener_err_max_),
        ener_curr(0), subdivision(1.0), subframe_amm(1), itr_counter(1)
    {
        gas_img.reserve(gas_.size());
        for(auto pcl : gas_){
            gas_img.push_back(new Particle_img(pcl));
        }
        ener_curr = verlet(0.0);
    }

    ~Verlet_adaptive(){
        for(auto pcl_img : gas_img){
            delete pcl_img;
        }
    }
    Verlet_adaptive (Verlet_adaptive const &) = delete;
    Verlet_adaptive& operator=(Verlet_adaptive const &) = delete;

    double verlet(double dt){
        for(auto pcl : gas){
            pcl->vel += pcl->acc * (dt/2);
            pcl->pos += pcl->vel * dt;
        }

        double ener = 0;
        reset_force(gas);
        reset_potential(gas);
        for (auto inter : inter_arr){
            ener += inter->calc(gas);
        }
        force_to_acc(gas);

        for(auto pcl : gas){
            pcl->vel += pcl->acc * (dt/2);
        }
        ener += calc_kinetic_energy(gas);
        return ener;
    }

    double operator()(){
        for(auto pcl_img : gas_img){
            pcl_img->pos_backup = pcl_img->pcl->pos;
            pcl_img->vel_backup = pcl_img->pcl->vel;
            pcl_img->acc_backup = pcl_img->pcl->acc;
        }

        itr_counter = 0;
        double ener_next = std::nan("0");
        while( !(std::abs(ener_next - ener_curr) < ener_err_max) ){
            for(auto pcl_img : gas_img){
                pcl_img->pcl->pos = pcl_img->pos_backup;
                pcl_img->pcl->vel = pcl_img->vel_backup;
                pcl_img->pcl->acc = pcl_img->acc_backup;
            }
            ener_next = ener_curr;
            for(size_t i=0; i<subdivision && (std::abs(ener_next - ener_curr) < ener_err_max); ++i){
                ener_next = verlet( dt/std::ceil(subdivision) );
                ++itr_counter;
            }
            subdivision *= 2.0;
        }
        subdivision /= 2.0;
        subframe_amm = std::ceil(subdivision);

        subdivision = (subdivision - 1.0)/2.0 + 1.0;
        ener_curr = ener_next;

        return ener_curr;
    }
};
