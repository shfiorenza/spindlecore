#ifndef _CGLASS_KINETOCHORE_H_
#define _CGLASS_KINETOCHORE_H_

#include "neighbor_list.hpp"
#include "object.hpp"
#include "triangle_mesh.hpp"

class Chromatid;

class Kinetochore : public Object {
  friend Chromatid;

private:
  Chromatid *parent_;

public:
  bool trip_update_on_removal_{false};
  int attachment_status_{-1};
  double n_exp_tot_{0.0};
  int n_bound_{0};

  TriMesh af_; // attachment factors;

  NeighborList neighbors_;

public:
  Kinetochore(unsigned long seed) : Object(seed) {
    SetSID(species_id::chromosome);
  }
  void Init(Chromatid *parent);
  void CreateRefVectors() {}
  void CreateBindingSites() {}
  void CreateTiangulatedMesh() {}

  // focus on this
  void UpdateNeighbors() {}

  void Step() {}

  void Update_1_2();
  void Stage_1_2_Probability() {}
  void DetermineAttachmentStatus() {}

  bool Insert_1_2() { return false; }
  int Remove_2_1_Fdep() { return 0; }
};
#endif