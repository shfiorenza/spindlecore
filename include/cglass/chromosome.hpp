#ifndef _CGLASS_CHROMOSOME_H_
#define _CGLASS_CHROMOSOME_H_

#include "chromatid.hpp"
#include "object.hpp"
#include "rng.hpp"
// #include "tracker.hpp"

// Safe arc cosine from -1 to 1
inline double safe_acos(double x) {
  if (x < -1.0)
    x = -1.0;
  else if (x > 1.0)
    x = 1.0;
  return acos(x);
}

#define SKIN_ 8.0

// By convention, the Object that 'Chromosome' represents
// is the centromere since it is a logical mid-way point
class Chromosome : public Object {
private:
  bool af_tip_crowd_{false};
  double af_xc_assemble_{0.0};
  double af_r0_{0.0};
  double af_tip_on_rate_assemble_{0.0};
  double rcutoff_{0.0};

  double delta_kmc_{0.0};

  chromosome_parameters *sparams_;

  std::vector<Chromatid> sisters_;
  std::vector<int> n_bound_;

public:
  int orientation_status_{0};

public:
  Chromosome(unsigned long seed);
  void Init(chromosome_parameters *sparams);

  void SetPosition(double const *const new_pos) {
    // printf("holla\n");
    Object::SetPosition(new_pos);
    // for (int i_sis{0}; i_sis < sisters_.size(); i_sis++) {
    //   double pos_sis[3]{new_pos[0], new_pos[1], new_pos[2]};
    //   for (int i{0}; i < 3; i++) {
    //     double sign{i_sis == 0 ? -1.0 : 1.0};
    //     pos_sis[i] += 0.8 * sign * diameter_;
    //   }
    //   sisters_[i_sis].SetPosition(pos_sis);
    // }
    // for (auto &&sis : sisters_) {
    //   sis.SetPosition(new_pos);
    // }
  }

  void Update_1_2_Probability();
  void DetermineAttachmentType();

  void KMC_1_2();
  void KMC_2_1();
  void KMC_2_1_FIndep();
  void KMC_2_1_FDep();
  void RunKMC();

  void Draw(std::vector<graph_struct *> &graph_array) {
    Object::Draw(graph_array);
    for (auto &&sis : sisters_) {
      sis.Draw(graph_array);
      // for (int i_dim{0}; i_dim < 3; i_dim++) {
      //   printf("U[%i] = %g\n", i_dim, sis.orientation_[i_dim]);
      // }
    }
  }
  void GetInteractors(std::vector<Object *> &ixors) {
    ixors.push_back(this);
    for (auto &&sis : sisters_) {
      ixors.push_back(&sis);
    }
  }
  void UpdatePosition();
};

#endif