#pragma once

#include <iostream>
#include <fstream>
#include <vector>

#include "Particle.hpp"
#include "Rigid_body.hpp"
#include "Interaction_base.hpp"

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

void read_gas    (std::vector<Particle *> &gas, std::string const &name);
void print_gas   (std::vector<Particle *> const &gas, std::ostream& os);
void print_frame (std::vector<Particle *> const &gas, int i);

void delete_everything(std::vector<Particle *> &gas,
                       std::vector<Rigid_body *> &bodies,
                       std::vector<Interaction_base *> &inter_arr);


class Report{
private:

public:
    struct Entry{
        size_t frame_i;

        double energy_tot;
        double energy_kin;
        double energy_pot;
        double delta_energy;
        Vec_3d momentum;
        Vec_3d angular_mom;
        size_t time_per_frame;
        size_t time_per_pcl;

        friend std::ostream& operator<<(std::ostream &os, const Entry &entry){
            os << "__________\n";
            os << "frame " << std::setw(4) << entry.frame_i << "|\n";
            os << "__________|_______________________\n";
            os << std::setprecision(5) << std::fixed << std::showpos;
            os << "total energy     " << " = " << entry.energy_tot << "\n";
            os << "kinetic energy   " << " = " << entry.energy_kin << "\n";
            os << "potential energy " << " = " << entry.energy_pot << "\n";
            os << "delta energy     " << " = " << entry.delta_energy << "\n";
            os << "momentum         " << " = " << entry.momentum << "\n";
            os << "angular momentum " << " = " << entry.angular_mom << "\n";
            os << "motion, sec:         " << entry.time_per_frame * 1.0E-9 << "\n";
            os << "motion per pcl, mcs: " << entry.time_per_pcl * 1.0E-3 << "\n";
            os << "__________________________________\n\n";
            return os;
        }
    };
    std::vector<Entry> entries;

    void add(Entry &entry){
        entries.push_back(entry);
    }

    template <typename T>
    std::vector<T> unroll_member(T Entry::* member) const{
        std::vector<T> arr;
        arr.reserve(entries.size());
        for (auto entry : entries){
            arr.push_back(entry.*member);
        }
        return arr;
    }

    friend std::ostream& operator<<(std::ostream &os, const Report &report){
        os << std::setprecision(15) << std::fixed;
        os << "total_energy"     << " = " << report.unroll_member(&Entry::energy_tot)         << "\n";
        os << "kinetic_energy"   << " = " << report.unroll_member(&Entry::energy_kin)     << "\n";
        os << "potential_energy" << " = " << report.unroll_member(&Entry::energy_pot)     << "\n";
        os << "delta_energy"     << " = " << report.unroll_member(&Entry::delta_energy)   << "\n";
        os << "frame_time"       << " = " << report.unroll_member(&Entry::time_per_frame) << "\n";
        os << "motion per pcl"   << " = " << report.unroll_member(&Entry::time_per_pcl)   << "\n";
        return os;
    }
};
