#include "../include/Interaction_Coulomb_FMM.hpp"


///FMM_octtree_node

FMM_octtree_node::FMM_octtree_node(Vec_3d pos, double r, std::vector<Particle *> &gas):
    pos(pos), r(r),
    gas(gas),
    gas_neighbour(0),
    gas_interaction(0)
{

}

///FMM implementation

Interaction_Coulomb_FMM_nlogn::Interaction_Coulomb_FMM_nlogn (double coulomb_constant_, double max_err_force_, double mp_rel_dist_, double mp_dir_time_ratio_):
    mc(0),
    rel_dist_4force(max_err_force_),
    mp_rel_dist(mp_rel_dist_),
    mp_dir_time_ratio(mp_dir_time_ratio_),
    gas_energy(0),
    must_calc_energy(false),
    must_calc_force(false),
    coulomb_constant(coulomb_constant_)
{
    mc.init(rel_dist_4force.get_degree(mp_rel_dist));
}

void Interaction_Coulomb_FMM_nlogn::calc_direct(std::vector<Particle *> const &inner, std::vector<Particle *> const &outer){
    for(auto pcl_1 : inner){
        for(auto pcl_2 : outer){
            double coulomb_factor = coulomb_constant * pcl_1->charge * pcl_2->charge;
            double dist = (pcl_1->pos - pcl_2->pos).len();
            if (must_calc_energy){
                double energy = coulomb_factor / dist;
                gas_energy += 1.0 * energy;
                pcl_1->potential += energy;
                pcl_2->potential += energy;
            }
            if (must_calc_force){
                Vec_3d force = (coulomb_factor / cube(dist)) * (pcl_1->pos - pcl_2->pos);
                pcl_1->force += force;
                pcl_2->force -= force;
            }
        }
    }
}

void Interaction_Coulomb_FMM_nlogn::calc_direct(std::vector<Particle *> const &gas){
    for(size_t i=0; i<gas.size(); ++i){
        Particle *pcl_1 = gas[i];
        for(size_t j=i+1; j<gas.size(); ++j){
            Particle *pcl_2 = gas[j];

            double coulomb_factor = coulomb_constant * pcl_1->charge * pcl_2->charge;
            double dist = (pcl_1->pos - pcl_2->pos).len();
            if (must_calc_energy){
                double energy = coulomb_factor / dist;
                gas_energy += 1.0 * energy;
                pcl_1->potential += energy;
                pcl_2->potential += energy;
            }
            if (must_calc_force){
                Vec_3d force = (coulomb_factor / cube(dist)) * (pcl_1->pos - pcl_2->pos);
                pcl_1->force += force;
                pcl_2->force -= force;
            }
        }
    }
}

void Interaction_Coulomb_FMM_nlogn::calc_multipole(std::vector<Particle *> const &inner, std::vector<Particle *> const &outer, Vec_3d origin){
    double max_r_inner = 0.0;
    for (auto pcl : inner){
        Vec_3d pos = pcl->pos - origin;
        if (max_r_inner < pos.sqr()) {
            max_r_inner = pos.sqr();
        }
    }
    max_r_inner = std::sqrt(max_r_inner);
    double min_r_outer = std::numeric_limits<double>::infinity();
    for (auto pcl : outer){
        Vec_3d pos = pcl->pos - origin;
        if (min_r_outer > pos.sqr()) {
            min_r_outer = pos.sqr();
        }
    }
    min_r_outer = std::sqrt(min_r_outer);
    double r_scale = min_r_outer;

    Multipole_table moment_inner(0);
    Multipole_table moment_inner_der_x(0);
    Multipole_table moment_inner_der_y(0);
    Multipole_table moment_inner_der_z(0);

    for (auto pcl : inner){
        Vec_3d pos = pcl->pos - origin;
        size_t p = rel_dist_4force.get_degree(min_r_outer/pos.len());
        moment_inner += mc.calc_regular(pos/r_scale, pcl->charge, p);
    }
    if (must_calc_force){
        moment_inner_der_x = mc.calc_regular_der_x(moment_inner);
        moment_inner_der_y = mc.calc_regular_der_y(moment_inner);
        moment_inner_der_z = mc.calc_regular_der_z(moment_inner);
    }
    Multipole_table moment_outer(0);
    Multipole_table moment_outer_der_x(0);
    Multipole_table moment_outer_der_y(0);
    Multipole_table moment_outer_der_z(0);
    for (auto pcl : outer){
        Vec_3d pos = pcl->pos - origin;
        size_t p = rel_dist_4force.get_degree(pos.len()/max_r_inner);
        Multipole_table moment_targ = mc.calc_irregular(pos/r_scale, pcl->charge, p);
        moment_outer += moment_targ;

        if(must_calc_energy){
            double energy = moment_inner * moment_targ;
            gas_energy += 0.5 * energy * (coulomb_constant / r_scale);
            pcl->potential += energy * (coulomb_constant / r_scale);
        }
        if(must_calc_force){
            Vec_3d force(0, 0, 0);
            force.x = moment_inner_der_x * moment_targ;
            force.y = moment_inner_der_y * moment_targ;
            force.z = moment_inner_der_z * moment_targ;
            pcl->force += force * (coulomb_constant / sqr(r_scale));
        }
    }
    if (must_calc_force){
        moment_outer_der_x = mc.calc_irregular_der_x(moment_outer);
        moment_outer_der_y = mc.calc_irregular_der_y(moment_outer);
        moment_outer_der_z = mc.calc_irregular_der_z(moment_outer);
    }
    for (auto pcl : inner){
        Vec_3d pos = pcl->pos - origin;
        size_t p = rel_dist_4force.get_degree(min_r_outer/pos.len());
        Multipole_table moment_targ = mc.calc_regular(pos/r_scale, pcl->charge, p);

        if(must_calc_energy){
            double energy = moment_outer * moment_targ;
            gas_energy += 0.5 * energy * (coulomb_constant / r_scale);
            pcl->potential += energy * (coulomb_constant / r_scale);
        }
        if(must_calc_force){
            Vec_3d force(0, 0, 0);
            force.x = moment_outer_der_x * moment_targ;
            force.y = moment_outer_der_y * moment_targ;
            force.z = moment_outer_der_z * moment_targ;
            pcl->force += force * (coulomb_constant / sqr(r_scale));
        }
    }
}

void Interaction_Coulomb_FMM_nlogn::calc_node_handler(FMM_octtree_node *curr){
    size_t size_inner = curr->gas.size();
    size_t size_outer = curr->gas_interaction.size();
    if (size_inner * size_outer >= (size_inner + size_outer) * mp_dir_time_ratio){
        calc_multipole(curr->gas, curr->gas_interaction, curr->pos);
    }else{
        calc_direct(curr->gas, curr->gas_interaction);
    }
}

void Interaction_Coulomb_FMM_nlogn::subdivide_gas(FMM_octtree_node *curr, std::vector<Particle *> ans[2][2][2]){
    for (auto pcl : curr->gas){
        int i_x=0, i_y=0, i_z=0;
        if (pcl->pos.x > curr->pos.x) { i_x=1; }
        if (pcl->pos.y > curr->pos.y) { i_y=1; }
        if (pcl->pos.z > curr->pos.z) { i_z=1; }
        ans[i_x][i_y][i_z].push_back(pcl);
    }
}

void Interaction_Coulomb_FMM_nlogn::complete_gas(FMM_octtree_node *parent, FMM_octtree_node *curr){
    double dist_interaction = curr->r * std::sqrt(3) * mp_rel_dist;
    for(auto pcl : parent->gas_neighbour){
        if (!pcl->processed){
            if( (pcl->pos - curr->pos).sqr() > sqr(dist_interaction)){
                curr->gas_interaction.push_back(pcl);
            }else{
                curr->gas_neighbour.push_back(pcl);
            }
        }
    }
}

void Interaction_Coulomb_FMM_nlogn::calc_handler_recursion(FMM_octtree_node *curr){
    calc_node_handler(curr);

    if(curr->gas_neighbour.size() == 1){
        curr->gas[0]->processed = true;
        return;
    }
    std::vector<Particle *> gas_son[2][2][2];
    subdivide_gas(curr, gas_son);
    for (int i_z=0; i_z<2; ++i_z){
        for (int i_y=0; i_y<2; ++i_y){
            for (int i_x=0; i_x<2; ++i_x){
                if (gas_son[i_x][i_y][i_z].size() != 0){
                    double new_r = curr->r / 2;
                    Vec_3d new_pos = curr->pos + new_r * Vec_3d((2*i_x-1), (2*i_y-1), (2*i_z-1));

                    FMM_octtree_node *son = new FMM_octtree_node(new_pos, new_r, gas_son[i_x][i_y][i_z]);
                    complete_gas(curr, son);

                    calc_handler_recursion(son);
                    delete son;
                }
            }
        }
    }
}

void Interaction_Coulomb_FMM_nlogn::calc_handler (std::vector<Particle *> &gas){
    gas_energy = 0;
    for(auto pcl : gas){
        pcl->processed = false;
    }

    std::pair<Vec_3d, double> cube = get_bounding_cube(gas);
    Vec_3d pos = cube.first;
    double r   = cube.second;
    FMM_octtree_node *root = new FMM_octtree_node(pos, r, gas);
    root->gas_neighbour = gas;

    calc_handler_recursion(root);
    delete root;
};

double Interaction_Coulomb_FMM_nlogn::calc_energy(std::vector<Particle *> &gas){
    must_calc_energy = true;
    must_calc_force  = false;
    calc_handler(gas);
    return gas_energy;
}

void Interaction_Coulomb_FMM_nlogn::calc_force (std::vector<Particle *> &gas){
    must_calc_energy = false;
    must_calc_force  = true;
    calc_handler(gas);
}

double Interaction_Coulomb_FMM_nlogn::calc (std::vector<Particle *> &gas){
    must_calc_energy = true;
    must_calc_force  = true;
    calc_handler(gas);
    return gas_energy;
}



void Interaction_Coulomb_FMM_nlogn::print_tree_recursion(FMM_octtree_node *curr, std::ostream& os){
    os << curr->pos.x << " ";
    os << curr->pos.y << " ";
    os << curr->pos.z << " ";
    os << curr->r << " ";
    os << curr->gas.size() << " ";
    os << curr->gas_neighbour.size() << " ";
    os << curr->gas_interaction.size() << " ";
    os << "\n";

    if(curr->gas_neighbour.size() == 1){
        return;
    }
    std::vector<Particle *> gas_son[2][2][2];
    subdivide_gas(curr, gas_son);
    for (int i_z=0; i_z<2; ++i_z){
        for (int i_y=0; i_y<2; ++i_y){
            for (int i_x=0; i_x<2; ++i_x){
                if (gas_son[i_x][i_y][i_z].size() != 0){
                    double new_r = curr->r / 2;
                    Vec_3d new_pos = curr->pos + new_r * Vec_3d((2*i_x-1), (2*i_y-1), (2*i_z-1));
//                    std::pair<Vec_3d, double> cube = get_bounding_cube(gas_son[i_x][i_y][i_z]);
//                    Vec_3d new_pos = cube.first;
//                    double new_r   = cube.second;

                    FMM_octtree_node *son = new FMM_octtree_node(new_pos, new_r, gas_son[i_x][i_y][i_z]);
                    complete_gas(curr, son);

                    print_tree_recursion(son, os);
                    delete son;
                }
            }
        }
    }
}

void Interaction_Coulomb_FMM_nlogn::print_tree (std::vector<Particle *> &gas, std::ostream& os){
    for(auto pcl : gas){
        pcl->processed = false;
    }
    os << "x y z r gas_size gas_neigh_size gas_inter_size\n";

    std::pair<Vec_3d, double> cube = get_bounding_cube(gas);
    Vec_3d pos = cube.first;
    double r   = cube.second;
    FMM_octtree_node *root = new FMM_octtree_node(pos, r, gas);
    root->gas_neighbour = gas;

    print_tree_recursion(root, os);
    delete root;
};
