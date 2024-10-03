#include "../include/Gas.hpp"

std::pair<Vec_3d, Vec_3d> get_bounding_box(std::vector<Particle *> &gas){
    Vec_3d coor_min;
    coor_min.x = coor_min.y = coor_min.z = std::numeric_limits<double>::infinity();
    Vec_3d coor_max = -1.0 * coor_min;
    for (auto pcl:gas){
        if (coor_min.x > pcl->pos.x){ coor_min.x = pcl->pos.x; }
        if (coor_min.y > pcl->pos.y){ coor_min.y = pcl->pos.y; }
        if (coor_min.z > pcl->pos.z){ coor_min.z = pcl->pos.z; }

        if (coor_max.x < pcl->pos.x){ coor_max.x = pcl->pos.x; }
        if (coor_max.y < pcl->pos.y){ coor_max.y = pcl->pos.y; }
        if (coor_max.z < pcl->pos.z){ coor_max.z = pcl->pos.z; }
    }
    Vec_3d pos = (coor_max + coor_min) / 2;
    Vec_3d dim = (coor_max - coor_min) / 2;
    return std::make_pair(pos, dim);
}

std::pair<Vec_3d, double> get_bounding_cube(std::vector<Particle *> &gas){
    std::pair<Vec_3d, Vec_3d> box = get_bounding_box(gas);
    double r = box.second.x;
    if (r < box.second.y){ r = box.second.y; }
    if (r < box.second.z){ r = box.second.z; }
    return std::make_pair(box.first, r);
}

Vec_3d get_vel_cm (std::vector<Particle *> &gas){
    Vec_3d momentum(0, 0, 0);
    double gas_mass = 0;
    for (auto pcl : gas){
        momentum += pcl->vel * pcl->mass;
        gas_mass += pcl->mass;
    }
    return momentum/gas_mass;
}

void reset_force (std::vector<Particle *> &gas){
    for(auto pcl : gas){
        pcl->force = Vec_3d(0, 0, 0);
    }
}

void reset_potential (std::vector<Particle *> &gas){
    for(auto pcl : gas){
        pcl->potential = 0;
    }
}

void force_to_acc (std::vector<Particle *> &gas){
    for(auto pcl : gas){
        pcl->acc = pcl->force / pcl->mass;
    }
}


void read_gas(std::vector<Particle *> &gas, const std::string &name){
    std::ifstream f_in;
    f_in.open (name);
    {
        std::string skip;
        std::getline(f_in, skip);
    }
    Vec_3d pos;
    Vec_3d vel;
    double mass;
    double charge;
    double potential;
    while (f_in >> pos.x){
        f_in >> pos.y >> pos.z >> vel.x >> vel.y >> vel.z >> mass >> charge >> potential;
        Particle *pcl = new Particle(pos, vel, mass, charge);
        pcl->potential = potential;
        gas.push_back(pcl);
    }
    f_in.close();
}

void print_gas(std::vector<Particle *> &gas, std::ostream& os){
    os << "x y z vx vy vz mass charge potential\n";
    for(auto pcl : gas){
        pcl->print_minimal(os);
    }
}

void print_frame(std::vector<Particle *> &gas, int i){
    std::ofstream f_out;
    f_out.open ("data/frames/frame" + std::to_string(i) + ".txt");
    print_gas(gas, f_out);
    f_out.close();
}


double calc_kinetic_energy(std::vector<Particle *> &gas){
    double energy = 0;
    for(auto pcl : gas){
        energy += pcl->vel.sqr() * pcl->mass;
    }
    return energy / 2;
}

double calc_potential_energy(std::vector<Particle *> &gas, std::vector<Interaction_base *> &inter_arr){
    double energy = 0;
    for (auto inter : inter_arr){
        energy += inter->calc_energy(gas);
    }
    return energy;
}

void verlet(std::vector<Particle *> &gas, std::vector<Interaction_base *> &inter_arr, double dt){
    for(auto pcl : gas){
        pcl->pos += pcl->vel * dt + pcl->acc * (dt*dt/2);
        pcl->vel += pcl->acc * (dt/2);
    }

    reset_force(gas);
    reset_potential(gas);
    for (auto inter : inter_arr){
        inter->calc_force(gas);
    }
    force_to_acc(gas);

    for(auto pcl : gas){
        pcl->vel += pcl->acc * (dt/2);
    }
}
