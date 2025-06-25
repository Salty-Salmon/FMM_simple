#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>
#include <set>
#include <complex>
#include <climits>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "include/Gas.hpp"
#include "include/Gas_IO.hpp"
#include "include/Interaction_Coulomb_FMM.hpp"
#include "include/Interaction_Coulomb_pairwise.hpp"
#include "include/Interaction_short.hpp"
#include "include/Integrator.hpp"

std::vector<Particle *> give_some_gas (){
    std::vector<Particle *> gas;
    int size_x =  2;
    int size_y = 41;
    int size_z = 41;
    double dist = 1.1;
    for (int i_x=0; i_x<size_x; ++i_x){
        for (int i_y=0; i_y<size_y; ++i_y){
            for (int i_z=0; i_z<size_z; ++i_z){
                double x = 00 + i_x + (i_y+i_z)%2 * 0.5;
                double y = 00 + std::sqrt(3)/2 * (i_y - (size_y)*0.5 + (i_z%2)/3.0);
                double z = 00 + std::sqrt(6)/3 * (i_z - (size_z)*0.5);

                Vec_3d pos = Vec_3d(x, y, z) * dist;
                Vec_3d vel = Vec_3d(0.3, 0, 0).cross(pos);//(0, 0, 0);
                double mass = 1.0;
                double charge = 1.0;

                double R = 10;
                if(y*y + z*z <= R*R){
                    gas.push_back(new Particle(pos, vel, mass, charge));
                }
            }
        }
    }
    gas.push_back(new Particle(Vec_3d(20, 0, 8), Vec_3d(-6, 0, 0), 20, 1));
    return gas;
}

std::vector<Rigid_body *> give_some_bodies (std::vector<Particle *> &gas){
    std::vector<Rigid_body *> bodies;

    std::vector<Particle *> rigid_gas;
    for (size_t i=0; i<gas.size()-1; ++i){
        rigid_gas.push_back(gas[i]);
    }
    bodies.push_back(new Rigid_body(rigid_gas));

    return bodies;
}

Vec_3d substract_CoM_vel (std::vector<Particle *> const &gas){
    Vec_3d momentum;
    double mass = 0;
    for (auto pcl : gas){
        momentum += pcl->vel * pcl->mass;
        mass += pcl->mass;
    }
    Vec_3d CoM_vel = momentum/mass;
    for (auto pcl : gas){
        pcl->vel -= CoM_vel;
    }
    return CoM_vel;
}

std::vector<Interaction_base *> init_interactions(){
    std::vector<Interaction_base *> interactions;

    //interactions.push_back(new Interaction_Coulomb_pairwise (-1.00) );
    //interactions.push_back(new Interaction_Coulomb_FMM_nlogn (-1.00, 1E-4, 1.5, 100.0) );
    interactions.push_back(new Interaction_6_12_smoothed (1.0, 1.0, 2.3, 1.7, 2) );

    return interactions;
}

int main()
{
    std::vector<Particle *> gas = give_some_gas();
    substract_CoM_vel(gas);
    std::vector<Rigid_body *> bodies = give_some_bodies(gas);
    //read_gas(gas, "data/saved_frames/planet8k.txt");

    std::vector<Interaction_base *> interactions = init_interactions();

    Runge_Kutta_4 integrator(gas, bodies, interactions, 0.1, 30);

    double ener_0 = integrator.get_energy_tot();

    Report report;

    print_frame(gas, 0);

    for (int i=1; i<=400; ++i){
        size_t t_0 = nanoseconds();
        integrator.integrate_dt();
        size_t t_1 = nanoseconds();

        Report::Entry entry;
        entry.frame_i = i;

        entry.energy_tot  = integrator.get_energy_tot();
        entry.energy_kin  = integrator.get_energy_kin();
        entry.energy_pot  = integrator.get_energy_pot();
        entry.momentum    = integrator.get_momentum();
        entry.angular_mom = integrator.get_angular_mom();
        entry.delta_energy = entry.energy_tot - ener_0;
        entry.time_per_frame = t_1 - t_0;
        entry.time_per_pcl = (t_1 - t_0) / (integrator.calc_force_per_step *
                                            integrator.substep_amm * gas.size());
        std::cout << entry;
        report.add(entry);


        print_frame(gas, i);
    }

    std::ofstream f_out("data/report_arrays.txt");
    f_out << report;
    f_out.close();

    delete_everything(gas, bodies, interactions);
    return 0;
}
