#include "../include/Gas_IO.hpp"

void read_gas(std::vector<Particle *> &gas, std::string const &name){
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

void print_gas(std::vector<Particle *> const &gas, std::ostream& os){
    os << "x y z vx vy vz mass charge potential body_id\n";
    Particle::set_print_minimal();
    for(auto pcl : gas){
        os << *pcl;
    }
}

void print_frame(std::vector<Particle *> const &gas, int i){
    std::ofstream f_out;
    f_out.open ("data/frames/frame" + std::to_string(i) + ".txt");
    print_gas(gas, f_out);
    f_out.close();
}

void delete_everything(std::vector<Particle *> &gas,
                       std::vector<Rigid_body *> &bodies,
                       std::vector<Interaction_base *> &inter_arr)
{
    for (auto pcl : gas){ delete pcl; }
    gas.clear();
    for (auto body : bodies){ delete body; }
    bodies.clear();
    for (auto interaction : inter_arr){ delete interaction; }
    inter_arr.clear();
}
