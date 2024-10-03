#include "../include/Vec_3d.hpp"

Vec_3d rotate_a_to_b(Vec_3d a, Vec_3d b, Vec_3d p){
    const double cos_min = 1E-9 - 1.0;

    a /= a.len();
    b /= b.len();
    double a_dot_b = a * b;
    if(a_dot_b < cos_min){
        return -p;
    }
    double inverse_cos_plus_1 = 1 / (1 + a_dot_b);

    double p_dot_a = p * a;
    double p_dot_b = p * b;
    double k_a = -1.0 * (p_dot_a + p_dot_b) * inverse_cos_plus_1;
    double k_b = (p_dot_a * (1 + 2 * a_dot_b) - p_dot_b) * inverse_cos_plus_1;

    return p + k_a*a + k_b*b;
}

double rand_uns(double min, double max) {
    static unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count();
    static std::default_random_engine e(seed);
    std::uniform_real_distribution<double> d(min, max);
    return d(e);
}

Vec_3d rand_unit_vec(){
    double phi = rand_uns(0, 4*std::acos(0));
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);

    double cos_theta = rand_uns(-1, 1);
    double sin_theta = std::sqrt(1 - sqr(cos_theta));

    Vec_3d ans(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
    return ans;
}

Vec_3d rand_unit_ball_vec(){
    double phi = rand_uns(0, 4*std::acos(0));
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);

    double cos_theta = rand_uns(-1, 1);
    double sin_theta = std::sqrt(1 - sqr(cos_theta));

    double r = std::pow(rand_uns(0, 1), 1.0 / 3.0);

    Vec_3d ans(r * sin_theta * cos_phi, r * sin_theta * sin_phi, r * cos_theta);
    return ans;
}

Vec_3d rand_unit_segment(Vec_3d axis, double theta_max){
    double phi = rand_uns(0, 4*std::acos(0));
    double cos_phi = std::cos(phi);
    double sin_phi = std::sin(phi);

    double cos_theta = rand_uns(std::cos(theta_max), 1);
    double sin_theta = std::sqrt(1 - sqr(cos_theta));

    Vec_3d deviation(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);
    Vec_3d ans = rotate_a_to_b(Vec_3d(0, 0, 1), axis, deviation);

    return ans;
}
