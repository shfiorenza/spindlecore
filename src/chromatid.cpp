// #include "cglass/chromatid.hpp"
#include "cglass/chromosome.hpp"
#include <iostream>

void Chromatid::Init(Chromosome *parent, double *u, double *v, double *w) {
  parent_ = parent;
  printf("hello\n");
  // zero_temperature_ = false;

  // sparams_ = sparams;
  color_ = parent_->sparams_->color;
  draw_ = draw_type::_from_string(parent_->sparams_->draw_type.c_str());
  length_ = parent_->sparams_->length;
  diameter_ = parent_->sparams_->diameter;
  // SF TODO check these gamma/diffusion coefficients
  gamma_par_ = 1.0 / (diameter_);
  gamma_perp_ = 2 * gamma_par_;
  gamma_rot_ = 3.0 / CUBE(diameter_);
  double noise_tr = parent_->sparams_->translational_noise;
  double noise_rot = parent_->sparams_->rotational_noise;
  diffusion_par_ = noise_tr * sqrt(24.0 * diameter_ / delta_);
  diffusion_perp_ = noise_tr * sqrt(48.0 * diameter_ / delta_);
  diffusion_rot_ = noise_rot * sqrt(8.0 * CUBE(diameter_) / delta_);
  for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
    u_[i_dim] = u[i_dim];
    v_[i_dim] = v[i_dim];
    w_[i_dim] = w[i_dim];
  }
  kc.Init(this);
  //   double pos[3] = {0, 0, 0};
  //   double u[3] = {0, 0, 0};
  //   Logger::Trace("Inserting CHROMATID object %d randomly", GetOID());
  //   rng_.RandomCoordinate(space_, pos, diameter_);
  //   rng_.RandomUnitVector(n_dim_, u);
  //   InsertAt(pos, u);
  //   InsertRod("random");
}

void Chromatid::UpdatePosition() {
  double new_pos[3];
  int sign{second_sister_ ? -1 : 1};
  for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
    new_pos[i_dim] = position_[i_dim] + 0.5 * sign * v_[i_dim] * diameter_;
  }
  kc.SetPosition(new_pos);
}