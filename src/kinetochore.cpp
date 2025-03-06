#include "cglass/kinetochore.hpp"
#include <iostream>

void Kinetochore::Init() {}

void Kinetochore::Update_1_2() {

  /*
  // Already have neighbor list
  n_exp_tot_ = 0.0;
  for (int isite = 0; isite < nsites_; ++isite) {
    //std::cout << "  af[" << isite << "]\n";
    n_exp_[isite] = 0.0;
    if (attach_[isite] != -1) {
      continue;
    }

    // Each site has it's own neighbor list and everything!
    for (nl_list::iterator p = neighbors_[isite].begin();
         p != neighbors_[isite].end(); ++p) {
      if (parent_->af_tip_distance_ > 0.0) {
        // Tip enhancement works slightly differently
        double dtip = parent_->af_tip_distance_;
        if (!p->duplicate_) {
          // Side of the microtubule
          double lside = lbond[p->label] - dtip;
          double rside[3] = {0.0};
          for (int i = 0; i < ndim_; ++i) {
            rside[i] =
                rbond[p->label][i] -
                0.5 * dtip * ubond[p->label][i]; // set the new r distance
          }
          double dr[3] = {0.0};
          for (int i = 0; i < ndim_; ++i) {
            dr[i] = (*anchors_)[isite].pos[i] - rside[i];
          }
          double binding_affinity = 0.0;
          // The side also cares about the SAC being on or off and adjusts its rate accordingly
          if (parent_->sac_status_ == 1) {
            binding_affinity =
                parent_->af_side_eps_eff_ * parent_->af_side_on_rate_;
          } else {
            binding_affinity =
                parent_->af_side_eps_eff_ * parent_->anaphase_rate_;
          }
          // Add the aurora B contribution
          if (parent_->aurora_b_effects_ == 1 &&
              parent_->chromosome_orientation_status_[cidx_] != 4 &&
              parent_->sac_status_ == 1) {
            binding_affinity = parent_->aurora_b_factor_ * binding_affinity;
          } else if (parent_->aurora_b_effects_ == 3 &&
                     parent_->sac_status_ == 1) {
            // Inter-kinetochore force-dependent on-rate
            double fmag_interkc = 0.0;
            for (int i = 0; i < 3; ++i) {
              fmag_interkc += SQR(parent_->f_interkc_[idx_][i]);
            }
            fmag_interkc = std::sqrt(fmag_interkc);
            // Determine if the interkinetochore force is towards the direction of ueff of the kinetochore, or not
            double interkcforce_dot_ueff =
                dot_product(3, parent_->f_interkc_[idx_], ueff_);

            double aurora_b_force_factor = 0.0;
            if (fmag_interkc <= parent_->aurora_b_force_plateau_) {
              aurora_b_force_factor = 1.0;
            } else if (interkcforce_dot_ueff > 0.0) {
              aurora_b_force_factor = 1.0;
            } else {
              aurora_b_force_factor =
                  1.0 /
                  (1.0 + parent_->aurora_b_force_dependence_ *
                             (fmag_interkc - parent_->aurora_b_force_plateau_));
              //std::cout << "kc[" << idx_ << "] fmag = " << fmag_interkc << ", direction = " << interkcforce_dot_ueff << std::endl;
              //std::cout << "   aurora b stabilization factor (onrate) = " << aurora_b_force_factor << std::endl;
            }
            binding_affinity = binding_affinity * aurora_b_force_factor;
          }
          p->value = binding_affinity *
                     Stage_1_2_Probability(ndim_, p->label, false, dr,
                                           ubond[p->label], lside);
          //if (p->value > 0.0) {
          //    std::cout << "Side contribution[" << p->label << "] = " << p->value << std::endl;
          //}
        } else {
          // Tip region
          double ltip = dtip;
          double rtip[3] = {0.0};
          for (int i = 0; i < ndim_; ++i) {
            rtip[i] = rbond[p->label][i] +
                      0.5 * lbond[p->label] * ubond[p->label][i] -
                      0.5 * dtip * ubond[p->label][i]; // set the new r distance
          }
          double dr[3] = {0.0};
          for (int i = 0; i < ndim_; ++i) {
            dr[i] = (*anchors_)[isite].pos[i] - rtip[i];
          }
          double binding_affinity = 0.0;
          if (parent_->properties_->bonds.poly_state[p->label] == GROW) {
            // Change based on the SAC status
            if (parent_->sac_status_ == 1) {
              binding_affinity =
                  parent_->af_tip_eps_eff_ * parent_->af_tip_on_rate_assemble_;
            } else {
              binding_affinity =
                  parent_->af_tip_eps_eff_ * parent_->anaphase_rate_;
            }
          } else {
            // Change the lockdown procedure and amount based on the SAC status
            if (parent_->sac_status_ == 1) {
              binding_affinity =
                  parent_->af_tip_eps_eff_ * parent_->af_tip_on_rate_disemble_;
            } else {
              binding_affinity =
                  parent_->af_tip_eps_eff_ * parent_->anaphase_rate_;
            }
          }
          // Add the aurora B contribution
          if (parent_->aurora_b_effects_ == 1 &&
              parent_->chromosome_orientation_status_[cidx_] != 4 &&
              parent_->sac_status_ == 1) {
            binding_affinity = parent_->aurora_b_factor_ * binding_affinity;
          } else if (parent_->aurora_b_effects_ == 3 &&
                     parent_->sac_status_ == 1) {
            // Inter-kinetochore force-dependent on-rate
            double fmag_interkc = 0.0;
            for (int i = 0; i < 3; ++i) {
              fmag_interkc += SQR(parent_->f_interkc_[idx_][i]);
            }
            fmag_interkc = std::sqrt(fmag_interkc);
            // Determine if the interkinetochore force is towards the direction of ueff of the kinetochore, or not
            double interkcforce_dot_ueff =
                dot_product(3, parent_->f_interkc_[idx_], ueff_);

            double aurora_b_force_factor = 0.0;
            if (fmag_interkc <= parent_->aurora_b_force_plateau_) {
              aurora_b_force_factor = 1.0;
            } else if (interkcforce_dot_ueff > 0.0) {
              aurora_b_force_factor = 1.0;
            } else {
              aurora_b_force_factor =
                  1.0 /
                  (1.0 + parent_->aurora_b_force_dependence_ *
                             (fmag_interkc - parent_->aurora_b_force_plateau_));
              //std::cout << "kc[" << idx_ << "] fmag = " << fmag_interkc << ", direction = " << interkcforce_dot_ueff << std::endl;
              //std::cout << "   aurora b stabilization factor (onrate) = " << aurora_b_force_factor << std::endl;
            }
            binding_affinity = binding_affinity * aurora_b_force_factor;
          }
          p->value = binding_affinity *
                     Stage_1_2_Probability(ndim_, p->label, true, dr,
                                           ubond[p->label], ltip);
          //if (p->value > 0.0) {
          //    std::cout << "Tip contribution[" << p->label << "] = " << p->value << std::endl;
          //}
        }
      } else {
        double dr[3] = {0.0};
        for (int i = 0; i < ndim_; ++i) {
          dr[i] = (*anchors_)[isite].pos[i] - rbond[p->label][i];
        }
        double binding_affinity = 0.0;
        binding_affinity =
            parent_->af_side_eps_eff_ * parent_->af_side_on_rate_;
        // Add the aurora B contribution
        if (parent_->aurora_b_effects_ == 1 &&
            parent_->chromosome_orientation_status_[cidx_] != 4 &&
            parent_->sac_status_ == 1) {
          binding_affinity = parent_->aurora_b_factor_ * binding_affinity;
        } else if (parent_->aurora_b_effects_ == 3 &&
                   parent_->sac_status_ == 1) {
          // Inter-kinetochore force-dependent on-rate
          double fmag_interkc = 0.0;
          for (int i = 0; i < 3; ++i) {
            fmag_interkc += SQR(parent_->f_interkc_[idx_][i]);
          }
          fmag_interkc = std::sqrt(fmag_interkc);
          // Determine if the interkinetochore force is towards the direction of ueff of the kinetochore, or not
          double interkcforce_dot_ueff =
              dot_product(3, parent_->f_interkc_[idx_], ueff_);

          double aurora_b_force_factor = 0.0;
          if (fmag_interkc <= parent_->aurora_b_force_plateau_) {
            aurora_b_force_factor = 1.0;
          } else if (interkcforce_dot_ueff > 0.0) {
            aurora_b_force_factor = 1.0;
          } else {
            aurora_b_force_factor =
                1.0 /
                (1.0 + parent_->aurora_b_force_dependence_ *
                           (fmag_interkc - parent_->aurora_b_force_plateau_));
            //std::cout << "kc[" << idx_ << "] fmag = " << fmag_interkc << ", direction = " << interkcforce_dot_ueff << std::endl;
            //std::cout << "   aurora b stabilization factor (onrate) = " << aurora_b_force_factor << std::endl;
          }
          binding_affinity = binding_affinity * aurora_b_force_factor;
        }
        p->value = binding_affinity *
                   Stage_1_2_Probability(ndim_, p->label, false, dr,
                                         ubond[p->label], lbond[p->label]);
      }

      // Check to see if the probability of attaching will be greater than 1.0, this will lead to
      // problems down the road
      if (p->value * parent_->delta_kmc_ > 1.0) {
        std::cerr << "ERROR: kMC on-rate probability for chromosomes is > 1.0, "
                     "p->value = "
                  << p->value << std::endl;
        std::cerr << "   Not exiting, continuing on\n";
        //exit(1);
      }

      n_exp_[isite] += p->value;
      n_exp_tot_ += p->value;
    }
  }
  */
}