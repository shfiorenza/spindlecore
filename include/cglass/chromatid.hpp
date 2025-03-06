#ifndef _CGLASS_CHROMATID_H_
#define _CGLASS_CRROMATID_H_

#include "kinetochore.hpp"
#include "spherocylinder.hpp"
class Chromosome;

class Chromatid : public BrRod {

private:
  Kinetochore kc;

  chromosome_parameters *sparams_;

  friend Chromosome;

protected:
  double v_[3];
  double w_[3];

public:
  Chromatid(unsigned long seed) : BrRod(seed), kc(seed) {
    printf("  NEW chromatid\n");
    SetSID(species_id::chromosome);
  }
  void Init(chromosome_parameters *sparams);

  void GetBodyFrame() {
    BrRod::GetBodyFrame();
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      v_[i_dim] = body_frame_[i_dim];
      w_[i_dim] = body_frame_[3 + i_dim];
    }
  }
};

#endif