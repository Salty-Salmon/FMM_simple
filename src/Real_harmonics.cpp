#include "../include/Real_harmonics.hpp"

void Multipole_calculator::init(size_t p_max_){
    if(recc_ll_11!=nullptr){
        delete[] recc_ll_11;
        delete[] recc_ll_12;
        delete[] recc_lm_10;
        delete[] recc_lm_20;
    }
    p_max = std::max(p_max_, (size_t)2);

    recc_ll_11 = new double[p_max];
    recc_ll_12 = new double[p_max];
    recc_lm_10 = new double[p_max*(p_max+1)/2];
    recc_lm_20 = new double[p_max*(p_max+1)/2];

    recc_ll_11[0] = 0;
    for (size_t l=1; l<p_max; ++l){
        recc_ll_11[l] = std::sqrt( (1.0*(2*l-1)) / (2*l) );
    }

    recc_ll_12[0] = 0;
    recc_ll_12[1] = 0;
    for (size_t l=2; l<p_max; ++l){
        recc_ll_12[l] = std::sqrt( (1.0*(2*l-1)) / (2*l-2) );
    }

    recc_lm_10[0] = 0;
    recc_lm_20[0] = 0;
    for (size_t l=1; l<p_max; ++l){
        recc_lm_10[l*(l+1)/2 + l-1] = 0;
        recc_lm_10[l*(l+1)/2 + l] = 0;
        recc_lm_20[l*(l+1)/2 + l-1] = 0;
        recc_lm_20[l*(l+1)/2 + l] = 0;
    }

    for (size_t l=2; l<p_max; ++l){
        for (size_t m=0; m<=l-2; ++m){
            recc_lm_10[l*(l+1)/2 + m] = (2*l-1) / std::sqrt(l*l - m*m);
            recc_lm_20[l*(l+1)/2 + m] = -1.0 * std::sqrt( 1.0 * (sqr(l-1)-m*m)/(l*l-m*m) );
        }
    }
    /// derivatives of regular harmonics
    weight_reg_der_x_dec_m = Multipole_table(p_max);
    for (int l=0; l<(signed)p_max; ++l){
        for (int m=-l; m<=l; ++m){
            int sign = (m > 0) - (m < 0);
            weight_reg_der_x_dec_m[l][m] = (std::sqrt( (l+m-1) * (l+m) ) / 2) * sign * (1 + (std::sqrt(2) - 1) * (m == 1));
        }
    }
    weight_reg_der_x_inc_m = Multipole_table(p_max);
    for (int l=0; l<(signed)p_max; ++l){
        for (int m=-l; m<=l; ++m){
            int sign = ((m+1 > 0) - (m+1 < 0));
            weight_reg_der_x_inc_m[l][m] = (std::sqrt( (l-m-1) * (l-m) ) / 2) * (-sign) * (1 + (std::sqrt(2) - 1) * (m == 0));
        }
    }
    weight_reg_der_y_dec_m = Multipole_table(p_max);
    for (int l=0; l<(signed)p_max; ++l){
        for (int m=-l; m<=l; ++m){
            int sign = (m > 0) - (m < 0);
            weight_reg_der_y_dec_m[l][m] = (std::sqrt( (l+m-1) * (l+m) ) / 2) * (-sign) * (1 - (m == 1));
        }
    }
    weight_reg_der_y_inc_m = Multipole_table(p_max);
    for (int l=0; l<(signed)p_max; ++l){
        for (int m=-l; m<=l; ++m){
            int sign = ((m+1 > 0) - (m+1 < 0));
            weight_reg_der_y_inc_m[l][m] = (std::sqrt( (l-m-1) * (l-m) ) / 2) * (-sign + std::sqrt(2) * (m == -1)) * (1 + (std::sqrt(2) - 1) * (m == 0));
        }
    }
    weight_reg_der_z = Multipole_table(p_max);
    for (int l=0; l<(signed)p_max; ++l){
        for (int m=-l; m<=l; ++m){
            weight_reg_der_z[l][m] = std::sqrt(sqr(l)-sqr(m));
        }
    }
    /// derivatives of irregular harmonics
    weight_irr_der_x_dec_m = Multipole_table(p_max);
    for (int l=0; l<(signed)p_max; ++l){
        for (int m=-l; m<=l; ++m){
            int sign = (m > 0) - (m < 0);
            weight_irr_der_x_dec_m[l][m] = (std::sqrt( (l-m+1) * (l-m+2) ) / 2) * sign * (1 + (std::sqrt(2) - 1) * (m == 1));
        }
    }
    weight_irr_der_x_inc_m = Multipole_table(p_max);
    for (int l=0; l<(signed)p_max; ++l){
        for (int m=-l; m<=l; ++m){
            int sign = ((m+1 > 0) - (m+1 < 0));
            weight_irr_der_x_inc_m[l][m] = (std::sqrt( (l+m+1) * (l+m+2) ) / 2) * (-sign) * (1 + (std::sqrt(2) - 1) * (m == 0));
        }
    }
    weight_irr_der_y_dec_m = Multipole_table(p_max);
    for (int l=0; l<(signed)p_max; ++l){
        for (int m=-l; m<=l; ++m){
            int sign = (m > 0) - (m < 0);
            weight_irr_der_y_dec_m[l][m] = (std::sqrt( (l-m+1) * (l-m+2) ) / 2) * (-sign) * (1 - (m == 1));
        }
    }
    weight_irr_der_y_inc_m = Multipole_table(p_max);
    for (int l=0; l<(signed)p_max; ++l){
        for (int m=-l; m<=l; ++m){
            int sign = ((m+1 > 0) - (m+1 < 0));
            weight_irr_der_y_inc_m[l][m] = (std::sqrt( (l+m+1) * (l+m+2) ) / 2) * (-sign + std::sqrt(2) * (m == -1)) * (1 + (std::sqrt(2) - 1) * (m == 0));
        }
    }
    weight_irr_der_z = Multipole_table(p_max);
    for (int l=0; l<(signed)p_max; ++l){
        for (int m=-l; m<=l; ++m){
            weight_irr_der_z[l][m] = -1.0 * std::sqrt(sqr(l+1)-sqr(m));
        }
    }
};

Multipole_calculator::Multipole_calculator(size_t p_max_):
    p_max(0),
    recc_ll_11(nullptr), recc_lm_10(nullptr), recc_lm_20(nullptr),
    weight_reg_der_x_dec_m(0),
    weight_reg_der_x_inc_m(0),
    weight_reg_der_y_dec_m(0),
    weight_reg_der_y_inc_m(0),
    weight_reg_der_z(0),
    weight_irr_der_x_dec_m(0),
    weight_irr_der_x_inc_m(0),
    weight_irr_der_y_dec_m(0),
    weight_irr_der_y_inc_m(0),
    weight_irr_der_z(0)
{
    init(p_max_);
};

Multipole_calculator::~Multipole_calculator(){
    delete[] recc_ll_11;
    delete[] recc_ll_12;
    delete[] recc_lm_10;
    delete[] recc_lm_20;
};

Multipole_table Multipole_calculator::calc_regular (Vec_3d r, double charge, size_t p){
    if (p > p_max){
        p = p_max;
    }
    if (p < 2){
        if (p==0){
            return Multipole_table(0);
        }
        if (p==1){
            Multipole_table ans(1);
            ans[0][0] = charge;
            return ans;
        }
    }
    double *poly_z  = new double[p*(p+1)/2];
    double *re_x_iy = new double[p];
    double *im_x_iy = new double[p];

    double r_sqr = r.sqr();
    double x = r.x;
    double y = r.y;
    double z = r.z;

    poly_z[0] = charge;
    poly_z[1] = poly_z[0] * z;
    poly_z[1+1] = poly_z[0] * recc_ll_11[1];

    for (size_t l=2; l<p; ++l){
        size_t l_pos = l*(l+1)/2;
        size_t l_pos_1 = l_pos - l;
        size_t l_pos_2 = l_pos_1 - (l-1);

        for (size_t m=0; m<=l-2; ++m){
            poly_z[l_pos + m] = recc_lm_10[l_pos + m]*z*poly_z[l_pos_1 + m] + recc_lm_20[l_pos + m]*r_sqr*poly_z[l_pos_2 + m];
        }
        poly_z[l_pos + l-1] = recc_ll_12[l]*poly_z[l_pos_1 + (l-2)];
        poly_z[l_pos + l]   = recc_ll_11[l]*poly_z[l_pos_1 + (l-1)];
    }

    re_x_iy[0] = std::sqrt(2);
    im_x_iy[0] = 0;
    for (size_t m=1; m<p; ++m){
        re_x_iy[m] = x*re_x_iy[m-1] - y*im_x_iy[m-1];
        im_x_iy[m] = x*im_x_iy[m-1] + y*re_x_iy[m-1];
    }

    Multipole_table ans(p);
    for (size_t l=0; l<p; ++l){
        size_t l_pos = l*(l+1)/2;
        ans[l][0] = poly_z[l_pos];
        for (size_t m=1; m<=l; ++m){
            ans[l][+m] = poly_z[l_pos + m] * re_x_iy[m];
            ans[l][-m] = poly_z[l_pos + m] * im_x_iy[m];
        }
    }
    delete[] poly_z;
    delete[] re_x_iy;
    delete[] im_x_iy;

    return ans;
};

Multipole_table Multipole_calculator::calc_spherical (Vec_3d r, double charge, size_t p){
    return calc_regular(r/r.len(), charge, p);
};

Multipole_table Multipole_calculator::calc_irregular (Vec_3d r, double charge, size_t p){
    return calc_regular(r/r.sqr(), charge/r.len(), p);
};


Multipole_table Multipole_calculator::calc_regular_der_x (Multipole_table const & mt){
    size_t p = mt.size();
    Multipole_table ans(p);
    for (int l=1; l<(signed)p; ++l){
        for (int m=-l+2; m<=l; ++m){
            ans[l][m] += weight_reg_der_x_dec_m[l][m] * mt[l-1][m-1];
        }
    }
    for (int l=1; l<(signed)p; ++l){
        for (int m=-l; m<=l-2; ++m){
            ans[l][m] += weight_reg_der_x_inc_m[l][m] * mt[l-1][m+1];
        }
    }
    return ans;
}

Multipole_table Multipole_calculator::calc_regular_der_y (Multipole_table const & mt){
    size_t p = mt.size();
    Multipole_table ans(p);
    for (int l=1; l<(signed)p; ++l){
        for (int m=-l+2; m<=l; ++m){
            ans[l][m] += weight_reg_der_y_dec_m[l][m] * mt[l-1][-(m-1)];
        }
    }
    for (int l=1; l<(signed)p; ++l){
        for (int m=-l; m<=l-2; ++m){
            ans[l][m] += weight_reg_der_y_inc_m[l][m] * mt[l-1][-(m+1)];
        }
    }
    return ans;
}

Multipole_table Multipole_calculator::calc_regular_der_z (Multipole_table const & mt){
    size_t p = mt.size();
    Multipole_table ans(p);
    for (int l=1; l<(signed)p; ++l){
        ans[l][-l] = 0;
        for (int m=-l+1; m<=l-1; ++m){
            ans[l][m] = weight_reg_der_z[l][m] * mt[l-1][m];
        }
        ans[l][l] = 0;
    }
    return ans;
}


Multipole_table Multipole_calculator::calc_irregular_der_x (Multipole_table const & mt){
    size_t p = mt.size();
    Multipole_table ans(p-1);
    if (p==0){ return Multipole_table(0); }
    for (int l=0; l<(signed)p-1; ++l){
        for (int m=-l; m<=l; ++m){
            ans[l][m] += weight_irr_der_x_dec_m[l][m] * mt[l+1][m-1];
            ans[l][m] += weight_irr_der_x_inc_m[l][m] * mt[l+1][m+1];
        }
    }
    return ans;
}

Multipole_table Multipole_calculator::calc_irregular_der_y (Multipole_table const & mt){
    size_t p = mt.size();
    Multipole_table ans(p-1);
    if (p==0){ return Multipole_table(0); }
    for (int l=0; l<(signed)p-1; ++l){
        for (int m=-l; m<=l; ++m){
            ans[l][m] += weight_irr_der_y_dec_m[l][m] * mt[l+1][-(m-1)];
            ans[l][m] += weight_irr_der_y_inc_m[l][m] * mt[l+1][-(m+1)];
        }
    }
    return ans;
}

Multipole_table Multipole_calculator::calc_irregular_der_z (Multipole_table const & mt){
    size_t p = mt.size();
    if (p==0){ return Multipole_table(0); }
    Multipole_table ans(p-1);
    for (int l=0; l<(signed)p-1; ++l){
        for (int m=-l; m<=l; ++m){
            ans[l][m] = weight_irr_der_z[l][m] * mt[l+1][m];
        }
    }
    return ans;
}
