#include "cglass/chromatid.hpp"
#include <iostream>

void Chromatid::Init(chromosome_parameters *sparams) {
  printf("hello\n");
  // zero_temperature_ = false;

  sparams_ = sparams;
  color_ = sparams->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = sparams->length;
  diameter_ = sparams->diameter;
  // SF TODO check these gamma/diffusion coefficients
  gamma_par_ = 1.0 / (diameter_);
  gamma_perp_ = 2 * gamma_par_;
  gamma_rot_ = 3.0 / CUBE(diameter_);
  double noise_tr = sparams->translational_noise;
  double noise_rot = sparams->rotational_noise;
  diffusion_par_ = noise_tr * sqrt(24.0 * diameter_ / delta_);
  diffusion_perp_ = noise_tr * sqrt(48.0 * diameter_ / delta_);
  diffusion_rot_ = noise_rot * sqrt(8.0 * CUBE(diameter_) / delta_);
  kc.Init();
  //   double pos[3] = {0, 0, 0};
  //   double u[3] = {0, 0, 0};
  //   Logger::Trace("Inserting CHROMATID object %d randomly", GetOID());
  //   rng_.RandomCoordinate(space_, pos, diameter_);
  //   rng_.RandomUnitVector(n_dim_, u);
  //   InsertAt(pos, u);
  //   InsertRod("random");
}