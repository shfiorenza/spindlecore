#ifndef _CGLASS_CHROMATID_H_
#define _CGLASS_CRROMATID_H_

#include "kinetochore.hpp"
#include "spherocylinder.hpp"
class Chromosome;

class Chromatid : public BrRod {
  friend Chromosome;
  friend Kinetochore;

private:
  bool second_sister_;

  Chromosome *parent_;
  Kinetochore kc;

  chromosome_parameters *sparams_;

protected:
  double r_[3];
  double u_[3];
  double u_eff_[3];
  double v_[3];
  double w_[3];

public:
  Chromatid(unsigned long seed) : BrRod(seed), kc(seed) {
    printf("  NEW chromatid\n");
    SetSID(species_id::chromosome);
  }
  void Init(Chromosome *chromo, double *u, double *v, double *w);

  void GetBodyFrame() {
    BrRod::GetBodyFrame();
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      v_[i_dim] = body_frame_[i_dim];
      w_[i_dim] = body_frame_[3 + i_dim];
    }
  }

  void UpdatePosition();
};

#endif