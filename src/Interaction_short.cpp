#include "../include/Interaction_short.hpp"

///Collision_octtree::Node

size_t Collision_octtree::Node::id_total = 0;

size_t Collision_octtree::Node::get_id(const Node *curr){
    if (curr){
        return curr->id;
    }
    return 0;
}

Collision_octtree::Node::Node(Vec_3d pos, double r, Node *parent):
    pos(pos), r(r), parent(parent), son{{},{}}, neighbour{{},{},{}}, id(++id_total)
{
    neighbour[1][1][1] = this;
};

std::vector<Collision_octtree::Node *> Collision_octtree::Node::get_neighbour_list(){
    std::vector<Node *> neighbour_list;
    neighbour_list.reserve(27);
    for (int n_z=0; n_z<3; ++n_z){
        for (int n_y=0; n_y<3; ++n_y){
            for (int n_x=0; n_x<3; ++n_x){
                Node *neighbour = this->neighbour[n_x][n_y][n_z];
                if (neighbour != NULL){
                    neighbour_list.push_back(neighbour);
                }
            }
        }
    }
    return neighbour_list;
}

///Collision_octtree
bool Collision_octtree::print_pointers = false;

void Collision_octtree::add_pcl(Collision_octtree::Node *curr, Particle *pcl){
    if (curr->r <= pcl_diam){
        curr->gas.push_back(pcl);
        return;
    }
    int i_x=0, i_y=0, i_z=0;
    if (pcl->pos.x > curr->pos.x) { i_x=1; }
    if (pcl->pos.y > curr->pos.y) { i_y=1; }
    if (pcl->pos.z > curr->pos.z) { i_z=1; }

    if (curr->son[i_x][i_y][i_z] == NULL){
        double new_r = curr->r / 2;
        Vec_3d new_pos = curr->pos + new_r * Vec_3d((2*i_x-1), (2*i_y-1), (2*i_z-1));

        curr->son[i_x][i_y][i_z] = new Node(new_pos, new_r, curr);

        for (int n_z=0; n_z<3; ++n_z){
            for (int n_y=0; n_y<3; ++n_y){
                for (int n_x=0; n_x<3; ++n_x){
                    int n_x_curr = (i_x + n_x + 1)/2;
                    int n_y_curr = (i_y + n_y + 1)/2;
                    int n_z_curr = (i_z + n_z + 1)/2;
                    Node *curr_neighbour = curr->neighbour[n_x_curr][n_y_curr][n_z_curr];
                    if(curr_neighbour != NULL){
                        int i_x_neig = (i_x + n_x + 1)%2;
                        int i_y_neig = (i_y + n_y + 1)%2;
                        int i_z_neig = (i_z + n_z + 1)%2;
                        Node *son_neighbour = curr_neighbour->son[i_x_neig][i_y_neig][i_z_neig];
                        if (son_neighbour != NULL){
                            curr->son[i_x][i_y][i_z]->neighbour[n_x][n_y][n_z] = son_neighbour;
                            son_neighbour->neighbour[2-n_x][2-n_y][2-n_z] = curr->son[i_x][i_y][i_z];
                        }
                    }
                }
            }
        }
    }
    add_pcl(curr->son[i_x][i_y][i_z], pcl);
};

void Collision_octtree::dtor_recursion(Node *curr){
    if (curr == NULL){
        return;
    }
    for (int i_z=0; i_z<2; ++i_z){
        for (int i_y=0; i_y<2; ++i_y){
            for (int i_x=0; i_x<2; ++i_x){
                dtor_recursion(curr->son[i_x][i_y][i_z]);
            }
        }
    }
    delete curr;
}

bool Collision_octtree::check_collision(Particle *pcl_1, Particle *pcl_2){
    return (pcl_1->pos - pcl_2->pos).sqr() < pcl_diam * pcl_diam;
}

void Collision_octtree::get_collisions_recursion(Node *curr, std::vector<std::pair<Particle *, Particle *> > &collisions){
    if (curr == NULL){
        return;
    }
    if (curr->gas.size()){
        std::vector<Collision_octtree::Node *> neighbour_list = curr->get_neighbour_list();
        while (curr->gas.size()){
            Particle *pcl_1 = curr->gas[curr->gas.size()-1];
            curr->gas.pop_back();
            for (auto neighbour : neighbour_list){
                for (auto pcl_2 : neighbour->gas){
                    if ( check_collision(pcl_1, pcl_2) ){
                        collisions.push_back(std::make_pair(pcl_1, pcl_2));
                    }
                }
            }
        }
    }
    for (int i_z=0; i_z<2; ++i_z){
        for (int i_y=0; i_y<2; ++i_y){
            for (int i_x=0; i_x<2; ++i_x){
                get_collisions_recursion(curr->son[i_x][i_y][i_z], collisions);
            }
        }
    }
}

void Collision_octtree::handle_collisions_recursion(Node *curr, Functor_base &functor){
    if (curr == NULL){
        return;
    }
    if (curr->gas.size()){
        std::vector<Collision_octtree::Node *> neighbour_list = curr->get_neighbour_list();
        while (curr->gas.size()){
            Particle *pcl_1 = curr->gas[curr->gas.size()-1];
            curr->gas.pop_back();
            for (auto neighbour : neighbour_list){
                for (auto pcl_2 : neighbour->gas){
                    if ( check_collision(pcl_1, pcl_2) ){
                        functor.handle_collision(pcl_1, pcl_2);
                    }
                }
            }
        }
    }
    for (int i_z=0; i_z<2; ++i_z){
        for (int i_y=0; i_y<2; ++i_y){
            for (int i_x=0; i_x<2; ++i_x){
                handle_collisions_recursion(curr->son[i_x][i_y][i_z], functor);
            }
        }
    }
}

void Collision_octtree::print_recursion (Node *curr, std::ostream& os){
    if (curr == NULL){
        return;
    }
    os << curr->pos.x << " " << curr->pos.y << " " << curr->pos.z << " " << curr->r << " " << curr->gas.size() << "\n";
    if (print_pointers){
        size_t (*get_id)(const Node *) = Collision_octtree::Node::get_id;
        os << get_id(curr) << "\n";
        for (int i_z=0; i_z<2; ++i_z){
            for (int i_y=0; i_y<2; ++i_y){
                os << "|";
                for (int i_x=0; i_x<2; ++i_x){
                    os << std::setw(2) << get_id(curr->son[i_x][i_y][i_z]) << " ";
                }
                os << "|   ";
            }
            os << "\n";
        }
        os << "\n";
        for (int n_z=0; n_z<3; ++n_z){
            for (int n_y=0; n_y<3; ++n_y){
                os << "|";
                for (int n_x=0; n_x<3; ++n_x){
                    os << std::setw(2) << get_id(curr->neighbour[n_x][n_y][n_z]) << " ";
                }
                os << "|   ";
            }
            os << "\n";
        }
        os << "\n\n";
    }
    for (int i_z=0; i_z<2; ++i_z){
        for (int i_y=0; i_y<2; ++i_y){
            for (int i_x=0; i_x<2; ++i_x){
                print_recursion(curr->son[i_x][i_y][i_z], os);
            }
        }
    }
};

Collision_octtree::Collision_octtree (std::vector<Particle *> &gas, double cell_size, double pcl_diam):
    cell_size(cell_size), pcl_diam(pcl_diam)
{
    std::pair<Vec_3d, double> cube = get_bounding_cube(gas);
    Vec_3d pos = cube.first;
    double r   = cube.second;

    r = pcl_diam * std::pow(2, std::ceil(std::log2(r/pcl_diam)) + 1E-9);
    root = new Node(pos, r, NULL);
    for (auto pcl : gas){
        add_pcl(root, pcl);
    }
}

Collision_octtree::~Collision_octtree (){
    dtor_recursion(root);
}

std::vector<std::pair<Particle *, Particle *> > Collision_octtree::get_collisions(){
    std::vector<std::pair<Particle *, Particle *> > collisions;
    get_collisions_recursion(root, collisions);
    return collisions;
}

void Collision_octtree::handle_collisions(Functor_base &functor){
    handle_collisions_recursion(root, functor);
}

std::ostream& operator<<(std::ostream& os, Collision_octtree & rha){
    os << "x y z r gas_size\n";
    rha.print_recursion(rha.root, os);
    return os;
}

void Collision_octtree::print_tree_paraview(int i){
    std::ofstream f_out;
    f_out.open ("data/frames/tree" + std::to_string(i) + ".txt");
    print_pointers = false;
    f_out << *this;
    f_out.close();
}

///smooth_func

smooth_func::smooth_func(double r_0, double r_1, unsigned int degree):
    r_0(r_0), r_1(r_1), degree(degree) {}

double smooth_func::slow_pow(double x, unsigned int n){
    double pow = 1.0;
    for (unsigned int i=0; i<n; ++i){
        pow *= x;
    }
    return pow;
}

double smooth_func::operator()(double x){
    double x_normed = (x-r_0)/(r_1-r_0);
    if(x_normed < 0){
        return 0;
    }
    if(x_normed > 1){
        return 1;
    }
    double x_pow_0 = slow_pow(x_normed,     degree);
    double x_pow_1 = slow_pow(1 - x_normed, degree);
    return x_pow_0 / (x_pow_0 + x_pow_1);
}

double smooth_func::derivative(double x){
    double x_normed = (x-r_0)/(r_1-r_0);
    if(x_normed < 0 || x_normed > 1){
        return 0;
    }
    double x_pow_0 = slow_pow(x_normed,     degree-1);
    double x_pow_1 = slow_pow(1 - x_normed, degree-1);
    return (degree*x_pow_0*x_pow_1) / (slow_pow(x_normed*x_pow_0 + (1-x_normed)*x_pow_1, 2) * (r_1 - r_0));
}


///Interaction_6_12_smoothed

Interaction_6_12_smoothed::Interaction_6_12_smoothed(double eps, double sigma, double r_max, double r_unmod, unsigned int degree):
    eps(eps), sigma(sigma), smooth_f(r_max, r_unmod, degree) { }

double Interaction_6_12_smoothed::calc_energy_pcl_pcl_unmod (Particle *pcl_1, Particle *pcl_2){
    double dist_sqr = (pcl_2->pos - pcl_1->pos).sqr();

    double sigma_div_dist_2  = (sigma*sigma) / dist_sqr;
    double sigma_div_dist_6  = sigma_div_dist_2 * sigma_div_dist_2 * sigma_div_dist_2;
    double sigma_div_dist_12 = sigma_div_dist_6 * sigma_div_dist_6;

    return 4 * eps * (sigma_div_dist_12 - sigma_div_dist_6);
}

Vec_3d Interaction_6_12_smoothed::calc_force_pcl_pcl_unmod (Particle *pcl_1, Particle *pcl_2){
    double dist_sqr = (pcl_2->pos - pcl_1->pos).sqr();
    double sigma_sqr = sigma*sigma;

    double sigma_div_dist_2  = sigma_sqr/dist_sqr;
    double sigma_div_dist_4  = sigma_div_dist_2 * sigma_div_dist_2;
    double sigma_div_dist_8  = sigma_div_dist_4 * sigma_div_dist_4;
    double sigma_div_dist_14 = sigma_div_dist_2 * sigma_div_dist_4 * sigma_div_dist_8;

    return (-24*eps/sigma_sqr) * (2*sigma_div_dist_14 - sigma_div_dist_8) * (pcl_2->pos - pcl_1->pos);
}

double Interaction_6_12_smoothed::calc_energy_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
    return calc_energy_pcl_pcl_unmod(pcl_1, pcl_2) * smooth_f((pcl_2->pos - pcl_1->pos).len());
}

Vec_3d Interaction_6_12_smoothed::calc_force_pcl_pcl (Particle *pcl_1, Particle *pcl_2){
    Vec_3d delta_pos = pcl_2->pos - pcl_1->pos;
    double dist = (delta_pos).len();
    return calc_force_pcl_pcl_unmod  (pcl_1, pcl_2) * smooth_f(dist) +
          (calc_energy_pcl_pcl_unmod (pcl_1, pcl_2) * smooth_f.derivative(dist)/dist) * delta_pos;

}

double Interaction_6_12_smoothed::calc_energy (std::vector<Particle *> &gas){
    double energy = 0;
    Collision_octtree tree(gas, smooth_f.r_0, smooth_f.r_0);
    std::vector<std::pair<Particle *, Particle *> > collisions = tree.get_collisions();
    for (auto pair : collisions) {
        energy += calc_energy_pcl_pcl(pair.first, pair.second);
    }
    return energy;
}

void Interaction_6_12_smoothed::calc_force (std::vector<Particle *> &gas){
    Collision_octtree tree(gas, smooth_f.r_0, smooth_f.r_0);
    Functor_calc_force functor(this);
    tree.handle_collisions(functor);
}

double Interaction_6_12_smoothed::calc (std::vector<Particle *> &gas){
    double energy = 0;
    Collision_octtree tree(gas, smooth_f.r_0, smooth_f.r_0);
    std::vector<std::pair<Particle *, Particle *> > collisions = tree.get_collisions();
    for (auto pair : collisions) {
        energy      += calc_energy_pcl_pcl (pair.first, pair.second);
        Vec_3d force = calc_force_pcl_pcl (pair.first, pair.second);
        pair.first->force  += force;
        pair.second->force -= force;
    }
    return energy;
}
