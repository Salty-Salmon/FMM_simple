#include "../include/Gas.hpp"

std::pair<Vec_3d, Vec_3d> get_bounding_box(std::vector<Particle *> const &gas){
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

std::pair<Vec_3d, double> get_bounding_cube(std::vector<Particle *> const &gas){
    std::pair<Vec_3d, Vec_3d> box = get_bounding_box(gas);
    double r = box.second.x;
    if (r < box.second.y){ r = box.second.y; }
    if (r < box.second.z){ r = box.second.z; }
    return std::make_pair(box.first, r);
}

void reset_force (std::vector<Particle *> const &gas){
    for(auto pcl : gas){
        pcl->force = Vec_3d(0, 0, 0);
    }
}

void reset_potential (std::vector<Particle *> const &gas){
    for(auto pcl : gas){
        pcl->potential = 0;
    }
}

void force_to_acc (std::vector<Particle *> const &gas){
    for(auto pcl : gas){
        pcl->acc = pcl->force / pcl->mass;
    }
}

double calc_kinetic_energy(std::vector<Particle *> const &gas){
    double energy = 0;
    for(auto pcl : gas){
        energy += pcl->vel.sqr() * pcl->mass;
    }
    return energy / 2;
}
Vec_3d calc_momentum (std::vector<Particle *> const &gas){
    Vec_3d momentum, angular_mom;
    for(auto pcl : gas){
        momentum += pcl->mass * pcl->vel;
    }
    return momentum;
}
Vec_3d calc_angular_mom (std::vector<Particle *> const &gas){
    Vec_3d angular_mom;
    for(auto pcl : gas){
        angular_mom += pcl->mass * pcl->pos.cross(pcl->vel);
    }
    return angular_mom;
}

double calc_potential_energy(std::vector<Particle *> const &gas,
                             std::vector<Interaction_base *> &inter_arr)
{
    double energy = 0;
    for (auto inter : inter_arr){
        energy += inter->calc_energy(gas);
    }
    return energy;
}

double init_redundant(std::vector<Particle *> const &gas,
                      std::vector<Rigid_body *> &bodies,
                      std::vector<Interaction_base *> &inter_arr)
{
    double ener = 0;
    reset_force(gas);
    reset_potential(gas);
    for (auto inter : inter_arr){
        ener += inter->calc(gas);
    }
    force_to_acc(gas);
    for (auto body : bodies){
        body->calc_force_acc();
        body->calc_torque();
    }

    ener += calc_kinetic_energy(gas);

    return ener;
}

double euler(std::vector<Particle *> const &gas,
             std::vector<Rigid_body *> &bodies,
             std::vector<Interaction_base *> &inter_arr,
             double dt, bool calc_ener)
{
    for (auto pcl : gas){
        if (!pcl->rigid_body_id){
            pcl->pos += pcl->vel * dt;
            pcl->vel += pcl->acc * dt;
        }
    }
    for (auto body : bodies){
        body->pos += body->vel * dt;
        body->vel += body->acc * dt;

        body->angular_mom += body->torque * dt;
        body->rotation += body->calc_rotation_der() * dt;
        body->rotation = body->rotation.orthogonalize();

        body->fill_parts_pos_vel();
    }

    double ener = 0;
    reset_force(gas);
    reset_potential(gas);
    for (auto inter : inter_arr){
        if(calc_ener){
            ener += inter->calc(gas);
        }else{
            inter->calc_force(gas);
        }
    }
    force_to_acc(gas);
    for (auto body : bodies){
        body->calc_force_acc();
        body->calc_torque();
    }

    if (calc_ener){
        ener += calc_kinetic_energy(gas);
    }
    return ener;
}

//double verlet(std::vector<Particle *> const &gas,
//              std::vector<Rigid_body *> &bodies,
//              std::vector<Interaction_base *> &inter_arr,
//              double dt, bool calc_ener)
//{
//    for(auto pcl : gas){
//        if (!pcl->rigid_body_id){
//            pcl->pos += pcl->vel * dt + pcl->acc * (dt*dt/2);
//            pcl->vel += pcl->acc * (dt/2);
//        }
//    }
//
//    double ener = 0;
//    reset_force(gas);
//    reset_potential(gas);
//    for (auto inter : inter_arr){
//        if(calc_ener){
//            ener += inter->calc(gas);
//        }else{
//            inter->calc_force(gas);
//        }
//
//    }
//    force_to_acc(gas);
//
//    for(auto pcl : gas){
//        if (!pcl->rigid_body_id){
//            pcl->vel += pcl->acc * (dt/2);
//        }
//    }
//    if (calc_ener){
//        ener += calc_kinetic_energy(gas);
//    }
//    return ener;
//}

