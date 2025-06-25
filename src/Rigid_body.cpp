#include "../include/Rigid_body.hpp"

const Rotation_matr Rotation_matr::identity = Rotation_matr ({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}});
double Rotation_matr::eps = 1E-10;

Rotation_matr::Rotation_matr(){
    for (size_t i=0; i<3; ++i){
        for (size_t j=0; j<3; ++j){
            matr[i][j] = 0;
        }
    }
}
Rotation_matr::Rotation_matr(double const (&matr_)[3][3]){
    for (size_t i=0; i<3; ++i){
        for (size_t j=0; j<3; ++j){
            matr[i][j] = matr_[i][j];
        }
    }
}
Rotation_matr::Rotation_matr(Rotation_matr const &matr_): Rotation_matr(matr_.matr) {

}

double* const Rotation_matr::operator[](size_t const &i){
    return matr[i];
}

Rotation_matr Rotation_matr::operator*(Rotation_matr const &rha) const{
    Rotation_matr ans;
    for (size_t i=0; i<3; ++i){
        for (size_t j=0; j<3; ++j){
            for (size_t k=0; k<3; ++k){
                ans.matr[i][j] += matr[i][k] * rha.matr[k][j];
            }
        }
    }
    return ans;
}
Rotation_matr& Rotation_matr::operator*=(Rotation_matr const &rha){
    *this = *this * rha;
    return *this;
}
Rotation_matr Rotation_matr::operator*(double const &rha) const{
    Rotation_matr ans;
    for (size_t i=0; i<3; ++i){
        for (size_t j=0; j<3; ++j){
            ans[i][j] = matr[i][j] * rha;
        }
    }
    return ans;
}
Rotation_matr operator*(double const &lha, Rotation_matr const &rha){
    return rha * lha;
}
Rotation_matr& Rotation_matr::operator*=(double const &rha){
    *this = *this * rha;
    return *this;
}
Rotation_matr Rotation_matr::operator/(double const &rha) const{
    return *this * (1/rha);
}
Rotation_matr& Rotation_matr::operator/=(double const &rha){
    *this = *this / rha;
    return *this;
}
Vec_3d Rotation_matr::operator*(Vec_3d const &rha) const{
    Vec_3d ans;
    for (size_t i=0; i<3; ++i){
        for (size_t j=0; j<3; ++j){
            ans[i] += matr[i][j] * rha[j];
        }
    }
    return ans;
}
Vec_3d operator*(Vec_3d const &lha, Rotation_matr const &rha){
    Vec_3d ans;
    for (size_t i=0; i<3; ++i){
        for (size_t j=0; j<3; ++j){
            ans[j] += lha[i] * rha.matr[i][j] ;
        }
    }
    return ans;
}

Rotation_matr Rotation_matr::operator+() const{
    Rotation_matr ans = *this;
    return ans;
}
Rotation_matr Rotation_matr::operator+(Rotation_matr const &rha) const{
    Rotation_matr ans;
    for (size_t i=0; i<3; ++i){
        for (size_t j=0; j<3; ++j){
            ans[i][j] = matr[i][j] + rha.matr[i][j];
        }
    }
    return ans;
}
Rotation_matr& Rotation_matr::operator+=(Rotation_matr const &rha){
    *this = *this + rha;
    return *this;
}
Rotation_matr Rotation_matr::operator-() const{
    Rotation_matr ans;
    for (size_t i=0; i<3; ++i){
        for (size_t j=0; j<3; ++j){
            ans[i][j] = -matr[i][j];
        }
    }
    return ans;
}
Rotation_matr Rotation_matr::operator-(Rotation_matr const &rha) const{
    Rotation_matr ans;
    return *this + (-rha);
}
Rotation_matr& Rotation_matr::operator-=(Rotation_matr const &rha){
    *this = *this - rha;
    return *this;
}

double Rotation_matr::norm() const{
    double sqr_sum = 0;
      for (size_t i=0; i<3; ++i){
        for (size_t j=0; j<3; ++j){
            sqr_sum += sqr(matr[i][j]);
        }
    }
    return std::sqrt(sqr_sum);
}
Rotation_matr Rotation_matr::transpose() const{
    Rotation_matr ans;
    for (size_t i=0; i<3; ++i){
        for (size_t j=0; j<3; ++j){
            ans[i][j] = matr[j][i];
        }
    }
    return ans;
}
Rotation_matr Rotation_matr::diag(Vec_3d const &vec){
    Rotation_matr ans;
    for (size_t i=0; i<3; ++i){
        ans[i][i] = vec[i];
    }
    return ans;
}
Rotation_matr Rotation_matr::orthogonalize() const{
    /**
    Bring to closest orthogonal in Frobenius norm.
    proven to converge in case ||Q_0||_F < 1.
    **/
    Rotation_matr A = *this;
    Rotation_matr Q = identity - A.transpose() * A;
    while (Q.norm() >= eps) {
        A = A * (identity + Q * (0.5 * identity + 0.375 * Q));
        Q = identity - A.transpose() * A;
    }
    return A;
}
std::pair<Rotation_matr, Vec_3d> Rotation_matr::diagonalize_q_form() const{
    std::pair<Rotation_matr, Vec_3d> ans (identity, Vec_3d());
    Rotation_matr q_form = *this;
    ///TODO add symmetricity check
    for (size_t i=0; i<3-1; ++i){     // for strictly upper
        for (size_t j=i+1; j<3; ++j){ // triangular elements
            /**
            diagonalize symm submatrix: |a b|
            using rotations             |b d|
            **/
            double a = q_form[i][i];
            double b = q_form[i][j];
            double d = q_form[j][j];

            double denom = std::sqrt(sqr(2*b) + sqr(a-d));
            if (denom >= eps){
                double alpha = std::atan2(a-d, 2*b);
                /// TODO: remake without trigonometric funcs
                double cos_alpha = std::cos(alpha);
                double sin_alpha = std::sin(alpha);

                Rotation_matr rotation_ij = Rotation_matr::identity;
                rotation_ij[i][i] =   rotation_ij[j][j] = cos_alpha;
                rotation_ij[j][i] = -(rotation_ij[i][j] = sin_alpha);

                ans.first *= rotation_ij;

                q_form = rotation_ij.transpose() * q_form * rotation_ij;
            }
        }
    }
    for (size_t i=0; i<3; ++i){
        ans.second[i] = q_form[i][i];
    }
    return ans;
}

Rotation_matr Rotation_matr::delta_ij(size_t i, size_t j){
    Rotation_matr ans;
    ans[i][j] = 1;
    return ans;
}
Rotation_matr Rotation_matr::rot_to_matr(Vec_3d const &lha){
    Rotation_matr ans ({
    {     0, -lha.z,  lha.y},
    { lha.z,      0, -lha.x},
    {-lha.y,  lha.x,      0} });
    return ans;
}

std::ostream& operator<<(std::ostream &os, const Rotation_matr &rha) {
    os << std::fixed << std::setprecision(Vec_3d::out_precision);
    size_t w = 3 + Vec_3d::out_precision;
    for (size_t i=0; i<3; ++i){
        os << "/|\\"[i];
        for (size_t j=0; j<3; ++j){
            os << std::setw(w) << rha.matr[i][j] << " ";
        }
        os << "\\|/"[i] << "\n";
    }
    return os;
}

///Rigid_body

size_t Rigid_body::rigid_body_id = 0;


double Rigid_body::calc_mass_fill_pos_vel_angmom(std::vector<Particle*>const& gas){
    double mass = 0;
    Vec_3d pos_sum;
    Vec_3d momentum;
    for (auto pcl : gas){
        mass += pcl->mass;
        pos_sum += pcl->mass * pcl->pos;
        momentum += pcl->mass * pcl->vel;
    }
    pos = pos_sum / mass;
    vel = momentum / mass;
    for (auto pcl : gas){
        angular_mom += pcl->pos.cross(pcl->mass * pcl->vel);
    }
    angular_mom -= pos.cross(momentum);
    return mass;
}

Vec_3d Rigid_body::calc_main_inertia_fill_rotation(std::vector<Particle *> const &gas){
    Rotation_matr inertia;
    for (auto pcl : gas){
        Vec_3d rel_pos = pcl->pos - pos;
        for (size_t i=0; i<3; ++i){
            for (size_t j=0; j<3; ++j){
                inertia[i][j] -= pcl->mass * rel_pos[i] * rel_pos[j];
            }
            inertia[i][i] += pcl->mass * rel_pos.sqr();
        }
    }
    std::pair<Rotation_matr, Vec_3d> rotation_and_main = inertia.diagonalize_q_form();

    rotation = rotation_and_main.first;
    Vec_3d main_inertia = rotation_and_main.second;

    if (main_inertia[0] <= Rotation_matr::eps &&
        main_inertia[1] <= Rotation_matr::eps){
        main_inertia = Vec_3d(1, 1, 1);
    }
    if (main_inertia[2] <= Rotation_matr::eps){
        main_inertia[2] = 1;
    }

    return main_inertia;
}
std::vector<Rigid_body::Part> Rigid_body::zip_parts(std::vector<Particle *> const &gas){
    std::vector<Rigid_body::Part> ans;
    ans.reserve(gas.size());
    for (auto pcl : gas){
        Vec_3d pcl_rel_pos = rotation.transpose() * (pcl->pos - pos);
        ans.push_back(Part(pcl, pcl_rel_pos));
        pcl->rigid_body_id = rigid_body_id;
    }
    return ans;
}

Rigid_body::Rigid_body(std::vector<Particle *> const &gas):
    mass(calc_mass_fill_pos_vel_angmom(gas)),
    main_inertia(calc_main_inertia_fill_rotation(gas)),
    parts(zip_parts(gas)),
    id(++rigid_body_id),
    main_inertia_inverse(Vec_3d(1.0/main_inertia[0],
                                1.0/main_inertia[1],
                                1.0/main_inertia[2]))
{
    calc_rotation_der();
}

Rotation_matr Rigid_body::calc_rotation_der(){
    Vec_3d angular_mom_rel = rotation.transpose() * angular_mom;
    Vec_3d angular_vel_rel = Rotation_matr::diag(main_inertia_inverse) * angular_mom_rel;
    angular_vel = rotation * angular_vel_rel;
    return Rotation_matr::rot_to_matr(angular_vel) * rotation;
}
void Rigid_body::fill_parts_pos_vel(){
    for (auto part : parts){
        part.pcl->pos = pos + rotation * part.rel_pos;
        part.pcl->vel = vel + angular_vel.cross(part.pcl->pos - pos);
    }
}
void Rigid_body::calc_force_acc(){
    force = Vec_3d();
    for (auto part : parts){
        force += part.pcl->force;
    }
    acc = force / mass;
}
void Rigid_body::calc_torque(){
    torque = Vec_3d();
    for (auto part : parts){
        torque += (part.pcl->pos - pos).cross(part.pcl->force);
    }
}

