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
#include "include/Interaction_Coulomb_FMM.hpp"
#include "include/Interaction_Coulomb_pairwise.hpp"
#include "include/Interaction_short.hpp"
#include "include/Interaction_external.hpp"
#include "include/Stopwatch.hpp"

template<typename T1, typename T2>
std::ostream &operator<<(std::ostream& os, std::pair<T1, T2> const &pair){
    os << "( " << pair.first << ", " << pair.second << ") ";
    return os;
}

template<typename T>
std::ostream &operator<<(std::ostream& os, std::vector<T> const &vec){
    os << "[";
    for (size_t i=0; i+1<vec.size(); ++i){
        os << vec[i] << ", ";
    }
    if(vec.size() != 0){
        os << vec[vec.size()-1];
    }
    os << "]";
    return os;
}

int main()
{
    std::vector<Particle *> gas;
    //read_gas(gas, "data/saved_frames/planet8k.txt");
    {
        int size_x = 2;
        int size_y = 2;
        int size_z = 5;
        double dist = 1.0 * 1.00;
        for (int i_y=0; i_y<size_y; ++i_y){
            for (int i_x=0; i_x<size_x; ++i_x){
                for (int i_z=0; i_z<size_z; ++i_z){
                    double x = i_x + (i_y+i_z)%2 * 0.5 + rand_uns(0, 1E-9);
                    double y = std::sqrt(3)/2 * (i_y - (size_y-1)*0.5 + (i_z%2)/3.0) + rand_uns(0, 1E-9);
                    double z = std::sqrt(6)/3 * (i_z - (size_z-1)*0.5) + rand_uns(0, 1E-9);

                    Vec_3d pos(dist * x, dist * y, dist * z);
                    Vec_3d vel(0, 0, 0);
                    double mass = 1.0;
                    double charge = 1.0;
                    gas.push_back(new Particle(pos, vel, mass, charge));
                }
            }
        }
    }
//    {
//        int size_x = 2;
//        int size_y = 2;
//        int size_z = 2;
//        double dist = 1.0;// * 1.12;
//        for (int i_y=0; i_y<size_y; ++i_y){
//            for (int i_x=0; i_x<size_x; ++i_x){
//                for (int i_z=0; i_z<size_z; ++i_z){
//                    double x = i_x + (i_y+i_z)%2 * 0.5;
//                    double y = std::sqrt(3)/2 * (i_y - (size_y-1)*0.5 + (i_z%2)/3.0);
//                    double z = std::sqrt(6)/3 * (i_z - (size_z-1)*0.5);
//
//                    Vec_3d pos(dist * x + 20, dist * y, dist * z);
//                    Vec_3d vel(0, 0, 0);
//                    double mass = 1.0;
//                    double charge = 1.0;
//                    gas.push_back(new Particle(pos, vel, mass, charge));
//                }
//            }
//        }
//    }
    {
        Vec_3d vel_cm = get_vel_cm(gas);
        for (auto pcl : gas){
            pcl->vel -= vel_cm;
        }
    }

    Interaction_Coulomb_FMM_nlogn inter_0 (-1.00, 1E-4, 1.5, 100.0);
    Interaction_Coulomb_pairwise  inter_1 (-1.00);
    Interaction_6_12_smoothed     inter_2 (1.0, 1.0, 3.0, 2.0, 2);

    std::vector<Interaction_base *> inter_arr;
    //inter_arr.push_back(&inter_0);
    inter_arr.push_back(&inter_1);
    inter_arr.push_back(&inter_2);

    double ener_potential = calc_potential_energy(gas, inter_arr);
    double ener_kinetic   = calc_kinetic_energy(gas);
    std::cout << "total potential energy:" << ener_potential << "\n";
    std::cout << "total energy:          " << ener_potential + ener_kinetic << "\n";
    double ener_err_max = 1E-4 * std::abs(ener_potential);
    double dt = 0.01;
    Verlet_adaptive integrator(gas, inter_arr, dt, ener_err_max);

    std::vector<double> energy_arr {ener_potential + ener_kinetic};
    std::vector<double> delta_energy_arr {0.0};
    std::vector<size_t> time_total_arr {0};
    std::vector<size_t> time_per_pcl_arr {0};
    std::vector<size_t> subframe_arr {0};
    std::vector<size_t> iteration_arr {0};

    print_frame(gas, 0);
    {
        std::ofstream f_out("data/frames/tree" + std::to_string(0) + ".txt");
        inter_0.print_tree(gas, f_out);
        f_out.close();
    }
    for (int i=1; i<=100000; ++i){
        size_t t_0 = Stopwatch::get_nanoseconds();
        //double energy = integrator();
        size_t subfr = 1;
        double energy = integrator.verlet(dt / subfr);
        for(size_t i=0; i<subfr-1; ++i){
            energy = integrator.verlet(dt / subfr);
        }
        integrator.subframe_amm = integrator.itr_counter = subfr;
//        for(auto pcl : gas){ pcl->vel -= pcl->vel * 0.1; }
        size_t t_1 = Stopwatch::get_nanoseconds();

        energy_arr.push_back(energy);
        delta_energy_arr.push_back(energy - ener_kinetic - ener_potential);
        time_total_arr.push_back( t_1 - t_0 );
        time_per_pcl_arr.push_back( (t_1 - t_0) / (integrator.itr_counter * gas.size()) );
        subframe_arr.push_back(integrator.subframe_amm);
        iteration_arr.push_back(integrator.itr_counter);

        std::cout << "___________\n";
        std::cout << "frame " << std::setw(4) << i << "|\n";
        std::cout << "__________|_______________________\n";
        std::cout.precision(6);
        std::cout << std::fixed;
        std::cout << "total energy" << " = " << energy_arr[i] << "\n";
        std::cout << "delta energy" << " = " << delta_energy_arr[i] << "\n";
        std::cout << "motion, sec:         " << time_total_arr[i] * 1.0E-9 << "\n";
        std::cout << "motion per pcl, mcs: " << time_per_pcl_arr[i] * 1.0E-3 << "\n";
        std::cout << "subframe_amm: " << subframe_arr[i] << "\n";
        std::cout << "itr_amm:      " << iteration_arr[i] << "\n";
        std::cout << "__________________________________\n\n";

        if (i%100000 == 0){
            i/=100000;
            print_frame(gas, i);

            {
                std::ofstream f_out("data/frames/tree" + std::to_string(i) + ".txt");
                inter_0.print_tree(gas, f_out);
                f_out.close();
            }
            i*=100000;
        }
    }

    {
        std::ofstream f_out("data/report_arrays.txt");
        f_out.precision(15);
        f_out << "total_energy" << " = " << energy_arr << "\n";
        f_out << "delta_energy" << " = " << delta_energy_arr << "\n";
        f_out << "frame_time" << " = " << time_total_arr << "\n";
        f_out << "motion per pcl" << " = " << time_per_pcl_arr << "\n";
        f_out << "subframe_amm" << " = " << subframe_arr << "\n";
        f_out << "iteration_amm" << " = " << iteration_arr << "\n";
        f_out.close();
    }

    for (auto pcl : gas){
        delete pcl;
    }
    gas.clear();

    return 0;
}
