#ifndef _CGLASS_KINETOCHORE_H_
#define _CGLASS_KINETOCHORE_H_

#include "object.hpp"

class Kinetochore : public Object {

private:
public:
  bool trip_update_on_removal_{false};
  int attachment_status_{-1};
  double n_exp_tot_{0.0};
  int n_bound_{0};

public:
  Kinetochore(unsigned long seed) : Object(seed) {
    SetSID(species_id::chromosome);
  }
  void Init();
  void CreateRefVectors() {}
  void CreateBindingSites() {}
  void CreateTiangulatedMesh() {}

  void UpdateNeighbors() {}

  void Step() {}

  void Update_1_2();
  void Stage_1_2_Probability() {}
  void DetermineAttachmentStatus() {}

  bool Insert_1_2() { return false; }
  int Remove_2_1_Fdep() { return 0; }
};
#endif