#include "cglass/chromosome_species.hpp"

ChromosomeSpecies::ChromosomeSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::chromosome);
}

void ChromosomeSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  if (GetNInsert() <= 0) {
    return;
  }
}

void ChromosomeSpecies::Reserve() {
  SpeciesBase::Reserve();
  members_.reserve(GetNInsert());
}

void ChromosomeSpecies::PopMember() { Species::PopMember(); }

void ChromosomeSpecies::AddMember() { Species::AddMember(); }

void ChromosomeSpecies::CheckSAC() {
  if (do_anaphase_) {
    //// Check if the SAC should be silenced
    //if (n_chromosomes_bioriented == nchromosomes_) {
    //    sac_bioriented_times_ += nkmc_;
    //} else {
    //    sac_bioriented_times_ = 0; // reset if there is a break
    //}
    if (n_bioriented_ == members_.size()) {
      // SF TODO incorporate
      // sac_bioriented_times_ += nkmc_; // Just accumulate the time
      sac_bioriented_times_ += 1; // Just accumulate the time
    }

    //// Make sure to turn off the Aurora B destabilization after anaphase
    //if (sac_status_ == 0) {
    //    for (int ic = 0; ic < nchromosomes_; ++ic) {
    //        chromosome_orientation_status_[ic] = 4; // Force into having ABK disabled
    //    }
    //}

    if ((sac_bioriented_times_ > sac_delay_) && (sac_status_ == 1)) {
      // Turn off the SAC
      sac_status_ = 0;
      std::cout << "Turning off the SAC, step: " << params_->i_step
                << std::endl;
      std::cout << "  SAC bi-orientation time: " << sac_bioriented_times_
                << std::endl;
      std::cout << "  SAC status: " << sac_status_ << std::endl;
      std::cout << "  bi-oriented: " << n_bioriented_ << std::endl;
    }
  }
}

void ChromosomeSpecies::UpdatePositions() {
  // Run brownian dynamics and update statistics
  n_bioriented_ = 0;
  for (auto &&chromo : members_) {
    chromo.UpdatePosition();
    chromo.Update_1_2_Probability();
    chromo.DetermineAttachmentType();
    if (chromo.orientation_status_ == 4) {
      n_bioriented_++;
    }
  }
  CheckSAC();
  // Do binding/unbinding
  // SF TODO include nkmc_ for chromo-specific timestep
  for (auto &&chromo : members_) {
    chromo.RunKMC();
  }
}

void ChromosomeSpecies::GetInteractors(std::vector<Object *> &ixors) {
  for (auto &&chromo : members_) {
    chromo.GetInteractors(ixors);
  }
}