#pragma once

#include "Vec_3d.hpp"
#include <iostream>

class Multipole_table{
private:
    size_t size_order;
    double *arr;
    //size_t size_moment;
public:
    enum Harmonic_type {
        regular, irregular
    };

    Multipole_table(size_t size_order_):size_order(size_order_), arr(nullptr){
        if (size_order != 0){
            arr = new double[size_order*size_order];
        }
        for (size_t i=0; i<size_order*size_order; ++i){
            arr[i] = 0;
        }
    };

    ~Multipole_table(){
        if (arr != nullptr){
            delete[] arr;
        }
    };
    Multipole_table(Multipole_table const &src):Multipole_table(src.size()){
        for (size_t i=0; i<size_order*size_order; ++i){
            arr[i] = src.arr[i];
        }
    };
    Multipole_table& operator=(Multipole_table const &src){
        if (this == &src) return *this;
        Multipole_table tmp(src);
        std::swap(size_order, tmp.size_order);
        std::swap(arr,        tmp.arr);
        return *this;
    };

    size_t size() const{
        return size_order;
    };

    double *operator[](size_t l) const{
        return arr + l*(l+1);
    };

    double ind(size_t ind) const{
        return arr[ind];
    };
    double &ind(size_t ind){
        return arr[ind];
    };

    Multipole_table &operator+=(Multipole_table const &lha){
        if (size_order >= lha.size_order){
            for (size_t i=0; i<lha.size_order*lha.size_order; ++i){
                arr[i] += lha.arr[i];
            }
        } else {
            double *arr_new = new double[lha.size_order*lha.size_order];
            for (size_t i=0; i<size_order*size_order; ++i){
                arr_new[i] = arr[i] + lha.arr[i];
            }
            for (size_t i=size_order*size_order; i<lha.size_order*lha.size_order; ++i){
                arr_new[i] = lha.arr[i];
            }
            size_order = lha.size_order;
            delete[] arr;
            arr = arr_new;
        }
        return *this;
    };
    Multipole_table operator+(Multipole_table const &lha) const{
        if (size_order < lha.size_order){
            return lha + *this;
        }
        Multipole_table ans(*this);
        ans += lha;
        return ans;
    };

    double operator *(Multipole_table const &lha) const{
        double ans_0 = 0;
        double ans_1 = 0;
        double ans_2 = 0;
        double ans_3 = 0;
        size_t span = std::min(size_order*size_order, lha.size_order*lha.size_order);

        size_t i=0;
        for (; i+3<span; i+=4){
            ans_0 += arr[i+0] * lha.arr[i+0];
            ans_1 += arr[i+1] * lha.arr[i+1];
            ans_2 += arr[i+2] * lha.arr[i+2];
            ans_3 += arr[i+3] * lha.arr[i+3];
        }
        for (; i<span; ++i){
            ans_0 += arr[i+0] * lha.arr[i+0];
        }
        return ans_0 + ans_1 + ans_2 + ans_3;
    };


    Multipole_table &operator*=(double const k){
        for (size_t i=0; i<size_order*size_order; ++i){
            arr[i] *= k;
        }
        return *this;
    };
    Multipole_table operator*(double const k) const{
        Multipole_table ans(*this);
        ans *= k;
        return ans;
    };
    friend Multipole_table operator*(double const k, Multipole_table const &rha){
        return rha * k;
    }

    Multipole_table &pow_this(double k, Harmonic_type mode){
        double k_pow = 1.0;
        if(mode == irregular){
            k = 1.0/k;
            k_pow *= k;
        }
        for (int l=0; l<(signed)size_order; ++l){
            for (int m=-l; m<=l; ++m){
                *this[l][m] *= k_pow;
            }
            k_pow *= k;
        }
        return *this;
    };
    Multipole_table pow(double k, Harmonic_type mode) const{
        Multipole_table ans(*this);
        ans.pow_this(k, mode);
        return ans;
    };

    void print(Multipole_table const &mt){
    printf("\n");
    for (size_t l=0; l<mt.size(); ++l){
        printf("l=%2lu  ", l);
        for (size_t i=0; i<(mt.size()-l-1)*8; ++i){
            printf(" ");
        }
        for (int m=-l; m<=(signed)l; ++m){
            printf("%+6.3lf  ", mt[l][m]);

        }
        printf("\n");
    }
    printf("\n");
}
};

class Multipole_calculator{
private:
    size_t p_max;
    double *recc_ll_11;
    double *recc_ll_12;
    double *recc_lm_10;
    double *recc_lm_20;
    Multipole_table weight_reg_der_x_dec_m;
    Multipole_table weight_reg_der_x_inc_m;
    Multipole_table weight_reg_der_y_dec_m;
    Multipole_table weight_reg_der_y_inc_m;
    Multipole_table weight_reg_der_z;

    Multipole_table weight_irr_der_x_dec_m;
    Multipole_table weight_irr_der_x_inc_m;
    Multipole_table weight_irr_der_y_dec_m;
    Multipole_table weight_irr_der_y_inc_m;
    Multipole_table weight_irr_der_z;

public:
    void init(size_t p_max_);
    Multipole_calculator(size_t p_max_);
    ~Multipole_calculator();

    Multipole_table calc_regular   (Vec_3d r, double charge, size_t p);
    Multipole_table calc_spherical (Vec_3d r, double charge, size_t p);
    Multipole_table calc_irregular (Vec_3d r, double charge, size_t p);

    Multipole_table calc_regular_der_x (Multipole_table const & mt);
    Multipole_table calc_regular_der_y (Multipole_table const & mt);
    Multipole_table calc_regular_der_z (Multipole_table const & mt);

    Multipole_table calc_irregular_der_x (Multipole_table const & mt);
    Multipole_table calc_irregular_der_y (Multipole_table const & mt);
    Multipole_table calc_irregular_der_z (Multipole_table const & mt);
};

struct Rel_dist_4force{
    double const eps;
    std::vector<double> boundary_rel_dist;

    Rel_dist_4force(double eps_):eps(eps_), boundary_rel_dist(2){
        boundary_rel_dist[1] = boundary_rel_dist[0] = std::numeric_limits<double>::infinity();
    }

    double get_rel_dist(size_t p){
        double r = std::pow(p/eps, 1.0/(p-1));
        for(int i=0; i<10; ++i){
            r = r + r * ( ((r+1)*p - 1)/(p * (p-1) * r) ) * ( std::log(std::pow(1/r, p-1) * (p + (p-1)/r)) - std::log(eps) );
        }
        return r;
    }

    size_t get_degree (double rel_dist){
        size_t p=1;
        while (boundary_rel_dist[p] > rel_dist){
            ++p;
            if (p == boundary_rel_dist.size()){
                boundary_rel_dist.push_back(get_rel_dist(p));
            }
        }
        return p;
    }
};
