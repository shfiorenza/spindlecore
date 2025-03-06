#ifndef _CGLASS_CHROMOSOME_SPECIES_H_
#define _CGLASS_CHROMOSOME_SPECIES_H_

#include "KMC/kmc.hpp"
#include "chromosome.hpp"
#include "species.hpp"

class ChromosomeSpecies : public Species<Chromosome, species_id::chromosome> {
private:
  bool do_anaphase_{false};
  int n_bioriented_{0};
  int sac_bioriented_times_{0};
  int sac_delay_{0};
  int sac_status_{0};

  bool *update_;
  std::string checkpoint_file_;
  std::vector<Object *> *objs_;

  // Tracker *tracker_ = nullptr;
  LookupTable lut_;

public:
private:
public:
  ChromosomeSpecies(unsigned long seed);
  void Init(std::string spec_name, ParamsParser &parser);
  void PopMember();
  void AddMember();

  void Reserve();

  void CheckSAC();
  void UpdatePositions();

  void GetInteractors(std::vector<Object *> &ixors);

  void GetLastInteractors(std::vector<Object *> &ixors) {
    members_.back().GetInteractors(ixors);
  }
};

#endif