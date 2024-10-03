#include "../include/Interaction_Coulomb_pairwise.hpp"

Interaction_Coulomb_pairwise::Interaction_Coulomb_pairwise (double factor): factor(factor) {}

double Interaction_Coulomb_pairwise::calc_energy_pcl_pcl(Particle *pcl_1, Particle *pcl_2){
    return factor * (pcl_1->charge * pcl_2->charge) / (pcl_1->pos - pcl_2->pos).len();
}

Vec_3d Interaction_Coulomb_pairwise::calc_force_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
    Vec_3d r = pcl_1->pos - pcl_2->pos;
    double r_len = r.len();
    return (factor * (pcl_1->charge * pcl_2->charge) / (r_len * r_len * r_len)) * r;
}

double Interaction_Coulomb_pairwise::calc_energy (std::vector<Particle *> &gas){
    double energy = 0;
    for (auto it_1 = gas.begin(); it_1 != gas.end(); ++it_1){
        Particle *pcl_1 = *it_1;
        for (auto it_2 = std::next(it_1); it_2 != gas.end(); ++it_2){
            Particle *pcl_2 = *it_2;
            double energy_1_2 = calc_energy_pcl_pcl(pcl_1, pcl_2);
            pcl_1->potential += energy_1_2;
            pcl_2->potential += energy_1_2;
            energy += energy_1_2;
        }
    }
    return energy;
}

void Interaction_Coulomb_pairwise::calc_force (std::vector<Particle *> &gas){
    for (auto it_1 = gas.begin(); it_1 != gas.end(); ++it_1){
        Particle *pcl_1 = *it_1;
        for (auto it_2 = std::next(it_1); it_2 != gas.end(); ++it_2){
            Particle *pcl_2 = *it_2;
            Vec_3d force = calc_force_pcl_pcl(pcl_1, pcl_2);
            pcl_1->force += force;
            pcl_2->force -= force;
        }
    }
}

double Interaction_Coulomb_pairwise::calc (std::vector<Particle *> &gas){
    calc_force(gas);
    return calc_energy(gas);
}
