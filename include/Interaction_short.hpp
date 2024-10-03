#pragma once

#include <cmath>

#include "Interaction_base.hpp"
#include "Gas.hpp"

class Collision_octtree{
public:
    class Functor_base{
    private:

    public:
        virtual void handle_collision (Particle *pcl_1, Particle *pcl_2) = 0;
        virtual ~Functor_base() = default;
    };
private:
    class Node{
    private:
        static size_t id_total;
    public:
        Vec_3d pos;
        double r;

        std::vector<Particle *> gas;

        Node *parent;
        Node *son[2][2][2];
        Node *neighbour[3][3][3];

        const size_t id;

        static size_t get_id(const Node *curr);
        Node(Vec_3d pos, double r, Node *parent);
        std::vector<Node *> get_neighbour_list();
    };
    void add_pcl(Node *curr, Particle *pcl);
    void dtor_recursion (Node *curr);
    void get_collisions_recursion(Node *curr, std::vector<std::pair<Particle *, Particle *> > &collisions);
    void handle_collisions_recursion(Node *curr, Functor_base &functor);

    void print_recursion (Node *curr, std::ostream& os);
public:
    const double cell_size;
    const double pcl_diam;
    Node *root;

    Collision_octtree (std::vector<Particle *> &gas, double cell_size, double pcl_diam);
    ~Collision_octtree ();
    Collision_octtree (Collision_octtree const &) = delete;
    Collision_octtree& operator=(Collision_octtree const &) = delete;

    bool check_collision (Particle *pcl_1, Particle *pcl_2);
    std::vector<std::pair<Particle *, Particle *> > get_collisions();
    void handle_collisions(Functor_base &functor);

    static bool print_pointers;
    friend std::ostream& operator<<(std::ostream& os, Collision_octtree & rha);
    void print_tree_paraview (int i);
};


struct smooth_func{
private:
    static double slow_pow(double x, unsigned int n);
public:
    double r_0;
    double r_1;
    unsigned int degree;

    smooth_func(double r_0, double r_1, unsigned int degree);
    double operator()(double x);
    double derivative(double x);
};

struct Interaction_6_12_smoothed: public Interaction_base{
private:
    class Functor_calc_force: public Collision_octtree::Functor_base{
    private:
        Interaction_6_12_smoothed *overlord;
    public:
        Functor_calc_force(Interaction_6_12_smoothed *overlord_):overlord(overlord_){};

        void handle_collision (Particle *pcl_1, Particle *pcl_2){
            Vec_3d force = overlord->calc_force_pcl_pcl (pcl_1, pcl_2);
            pcl_1->force += force;
            pcl_2->force -= force;
        };
    };
public:
    double eps;
    double sigma;
    smooth_func smooth_f;

    Interaction_6_12_smoothed(double eps, double sigma, double r_max, double r_unmod, unsigned int degree);
    double calc_energy_pcl_pcl_unmod (Particle *pcl_1, Particle *pcl_2);
    Vec_3d calc_force_pcl_pcl_unmod  (Particle *pcl_1, Particle *pcl_2);
    double calc_energy_pcl_pcl (Particle *pcl_1, Particle *pcl_2);
    Vec_3d calc_force_pcl_pcl  (Particle *pcl_1, Particle *pcl_2);
    double calc_energy (std::vector<Particle *> &gas);
    void   calc_force  (std::vector<Particle *> &gas);
    double calc        (std::vector<Particle *> &gas);
};
