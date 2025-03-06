#include "cglass/chromosome.hpp"
#include <iostream>

Chromosome::Chromosome(unsigned long seed) : Object(seed) {
  printf("NEW chromosome\n");
  SetSID(species_id::chromosome);
  sisters_.emplace_back(Chromatid(seed));
  sisters_.emplace_back(Chromatid(seed));
}

void Chromosome::Init(chromosome_parameters *sparams) {

  sparams_ = sparams;
  color_ = sparams->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = 0.0; //sparams->length;
  diameter_ = sparams->diameter;
  double pos[3] = {0, 0, 0};
  double u[3] = {0, 0, 0};
  rng_.RandomCoordinate(space_, pos, diameter_);
  rng_.RandomUnitVector(n_dim_, u);
  InsertAt(pos, u);
  // SF TODO change this to GetBodyFrame here and in Centromere()
  // Get vectors that define plane perpendicular to u
  double v[3] = {0, 0, 0}; // aux vector; defines plane of SPB together with w
  double w[3] = {0, 0, 0}; // aux vector; defines plane of SPB together with v
  // SF: not sure if this is general or b/c we initialize aligned to zhat
  double xhat[3] = {1.0, 0.0, 0.0};
  double yhat[3] = {0.0, 1.0, 0.0};
  if (1.0 - ABS(u[0]) > 1e-2) {
    cross_product(u, xhat, v, 3);
  } else {
    cross_product(u, yhat, v, 3);
  }
  double norm_factor{sqrt(1.0 / dot_product(3, v, v))};
  for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
    v[i_dim] *= norm_factor;
  }
  cross_product(u, v, w, 3);
  double pos_sis_one[3];
  double pos_sis_two[3];
  for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
    pos_sis_one[i_dim] = pos[i_dim] + v[i_dim] * diameter_ / 2.0;
    pos_sis_two[i_dim] = pos[i_dim] - v[i_dim] * diameter_ / 2.0;
  }
  for (auto &&sis : sisters_) {
    sis.Init(sparams);
  }
  sisters_[0].InsertAt(pos_sis_one, u);
  sisters_[1].InsertAt(pos_sis_two, u);
}

void Chromosome::Update_1_2_Probability() {

  for (auto &&sis : sisters_) {
    sis.kc.Update_1_2();
  }
}

void Chromosome::DetermineAttachmentType() {

  // Ask the kinetochores to determine their attachment status
  //std::cout << "Step[" << properties->i_current_step << "] chromosomes determining attachment type\n";
  int kc0 = 0;
  int kc1 = 1;
  sisters_[kc0].kc.DetermineAttachmentStatus();
  sisters_[kc1].kc.DetermineAttachmentStatus();
  // Now, figure out what kind of attachment we happen to be
  if (sisters_[kc0].kc.attachment_status_ == -1 &&
      sisters_[kc1].kc.attachment_status_ == -1) {
    // Unattached means both are -1
    orientation_status_ = 0;
  } else if (sisters_[kc0].kc.attachment_status_ == 2 ||
             sisters_[kc1].kc.attachment_status_ == 2) {
    // Merotelic means either of them is merotelic
    orientation_status_ = 2;
  } else if (sisters_[kc0].kc.attachment_status_ ==
             sisters_[kc1].kc.attachment_status_) {
    // Syntelic (they are attached to the same pole)
    orientation_status_ = 3;
  } else {
    // Remaining options are to be amphitelic (attached to different poles), or monotelic (only one is attached)
    if (sisters_[kc0].kc.attachment_status_ == -1 ||
        sisters_[kc1].kc.attachment_status_ == -1) {
      // One of them is unattached, and so this must be a monotelic state
      orientation_status_ = 1;
    } else {
      // Only remaining option is to be in amphitelic attachment
      orientation_status_ = 4;
    }
  }
  //std::cout << "Attachment Status: " << chromosome_orientation_status_[ic] << std::endl;
}

void Chromosome::KMC_1_2() {

  for (auto &&sis : sisters_) {
    // How many attach?
    int Ntot = rng_.RandomPoisson(sis.kc.n_exp_tot_ * delta_kmc_);
    // Now do the attachment
    for (int itrial = 0; itrial < Ntot; ++itrial) {
      //if (properties->i_current_step >= 7577191) {
      //    std::cout << "(" << properties->i_current_step << ")";
      //    std::cout << "KC[" << ikc << "] attaching trial: " << itrial << std::endl;
      //    std::cout << "  nbound before new attach: " << n_bound_[ikc] << std::endl;
      //}
      bool successful_insert = sis.kc.Insert_1_2();
      if (successful_insert) {
        sis.kc.n_bound_++;
        //PrintFrame();
      }
    }
  }
}

void Chromosome::KMC_2_1() {
  if (af_xc_assemble_ == 0.0 && af_r0_ == 0.0) {
    KMC_2_1_FIndep();
  } else {
    KMC_2_1_FDep();
  }
}

void Chromosome::KMC_2_1_FIndep() {

  // SF NOTE: apparently this is never supposed to run. lol. remove?
  std::cout << "CM Shouldn't be here...\n";
  return;
  // exit(1);
  /*
  // Force independent, just detach via probability distribution
  unsigned int noff = gsl_ran_binomial(
      rng_, af_tip_on_rate_assemble_ * delta_kmc_, n_bound_[ikc]);

  for (unsigned int ioff = 0; ioff < noff; ++ioff) {
    //std::cout << "KC[" << ikc << "] detaching trial: "<< ioff << std::endl;
    Kinetochore *kc = &(kinetochores_[ikc]);
    bool did_remove = kc->Remove_2_1(parameters, properties);

    if (!did_remove) {
      std::cout << "Something funny happened during KC detachment\n";
      exit(1);
    }
  }

  // Assuming everythign went well, and we removed, rebulid the neighbor list
  // if need be
  if (noff > 0) {
    Kinetochore *kc = &(kinetochores_[ikc]);
    kc->trip_update_on_removal_ = true;
    //double r_cutoff2 = SQR(rcutoff_ + SKIN_);
    //kc->UpdateNeighbors(parameters->n_dim,
    //                    parameters->n_periodic,
    //                    properties->bonds.n_bonds,
    //                    properties->unit_cell.h,
    //                    properties->bonds.r_bond,
    //                    properties->bonds.s_bond,
    //                    properties->bonds.u_bond,
    //                    properties->bonds.length,
    //                    r_cutoff2);
  }
  */
}

void Chromosome::KMC_2_1_FDep() {

  for (auto &&sis : sisters_) {
    int noff = sis.kc.Remove_2_1_Fdep();

    // SF NOTE: this ensures that only ONE microtubule tip is bound to each AF
    // SF NOTE: I could put this in chromo_MGMT, but should rewrite using C-GLASS tools
    // Check for crowding effects
    if (af_tip_crowd_) {
      printf("Need to implement tip crowding FX\n");
      exit(1);
      /*
    for (int ikc2 = 0; ikc2 < nkcs_; ++ikc2) {
      Kinetochore *kc2 = &(kinetochores_[ikc2]);

      for (int isite1 = 0; isite1 < naf_; ++isite1) {
        for (int isite2 = 0; isite2 < naf_; ++isite2) {
          if (isite1 == isite2 && ikc == ikc2) {
            continue; // Don't worry about ourselves
          }
          // If they aren't both attached, continue
          if (kc->attach_[isite1] == -1 && kc2->attach_[isite2] == -1) {
            continue;
          }

          int ibond1 = kc->attach_[isite1];
          int ibond2 = kc2->attach_[isite2];

          // Check if on the same bond...
          if (ibond1 != ibond2) {
            continue;
          }

          double dstab1 =
              properties->bonds.length[ibond1] - kc->cross_pos_[isite1];
          double dstab2 =
              properties->bonds.length[ibond2] - kc2->cross_pos_[isite2];

          // Only check unbinding of myself (isite1 and ikc/kc)
          if (dstab2 < af_tip_distance_ && dstab1 < af_tip_distance_) {
            if (dstab1 > dstab2) {
              // Unbind me
              //std::cout << "KC[" << ikc << "]{" << isite1 << "} MT[" << ibond1 << "] position: " << kc->cross_pos_[isite1]
              //    << ", crowding against KC[" << ikc2 << "]{" << isite2 << "} position: " << kc2->cross_pos_[isite2] << std::endl;
              n_bound_[ikc]--;
              kc->attach_[isite1] = -1;
              kc->cross_pos_[isite1] = 0.0;
              noff++;
            }
          }
        }
      }
    }
  */
    }
    if (noff > 0) {
      //std::cout << "KC[" << ikc << "] detached trials: " << noff << std::endl;
      sis.kc.trip_update_on_removal_ = true;
      //PrintFrame();
      //double r_cutoff2 = SQR(rcutoff_ + SKIN_);
      //kc->UpdateNeighbors(parameters->n_dim,
      //                    parameters->n_periodic,
      //                    properties->bonds.n_bonds,
      //                    properties->unit_cell.h,
      //                    properties->bonds.r_bond,
      //                    properties->bonds.s_bond,
      //                    properties->bonds.u_bond,
      //                    properties->bonds.length,
      //                    r_cutoff2);
    }
  }
}

void Chromosome::RunKMC() {

  // Randomly decide to do detach->attach or reverse
  int g[2] = {0, 1};
  for (int i = 0; i < 2; ++i) {
    int j{(int)rng_.RandomInt(2)};
    // int j = gsl_rng_uniform_int(properties->rng.r, 2);
    int swap = g[i];
    g[i] = g[j];
    g[j] = swap;
  }
  for (int i = 0; i < 2; ++i) {
    switch (g[i]) {
    case 0:
      KMC_1_2();
      break;
    case 1:
      KMC_2_1();
      break;
    }
  }
  // Check if the kinetochores need to rebuild their neighbor list
  for (int ikc = 0; ikc < sisters_.size(); ++ikc) {
    Kinetochore *kc = (&sisters_[ikc].kc);
    if (kc->trip_update_on_removal_) {
      kc->trip_update_on_removal_ = false;
      double r_cutoff2 = SQR(rcutoff_ + SKIN_);
      kc->UpdateNeighbors();
    }
  }
}

void Chromosome::UpdatePosition() {

  // ! SF TODO check interKC forces; currently unaccounted for

  // SF TODO check:
  // orientation -> v_
  // v_ -> orientation

  // First: calculate and add force/torques to chromatin sisters
  // SF TODO should forces be zeroed here??
  for (auto &&sis : sisters_) {
    sis.ZeroForce();
    sis.GetBodyFrame();
  }
  int n_dim{params_->n_dim};
  double k, r0, kth, kv;
  //SF TODO put these in parameter structure/yaml file
  k = 6.0;
  r0 = 4.0;
  kth = 450.0;
  kv = 450.0;
  double r[3], rhat[3];
  // Find the initial r vector r = r_A - r_B
  double rmag2 = 0.0;
  for (int i = 0; i < n_dim; ++i) {
    r[i] = sisters_[0].position_[i] - sisters_[1].position_[i];
    rmag2 += SQR(r[i]);
  }
  double rmag = sqrt(rmag2);
  for (int i = 0; i < n_dim; ++i) {
    rhat[i] = r[i] / rmag;
  }

  // Calculate angular stuff
  // double dotA = dot_product(n_dim, sisters_[0].orientation_, rhat);
  // double dotB = dot_product(n_dim, sisters_[1].orientation_, rhat);
  double dotA = dot_product(n_dim, sisters_[0].v_, rhat);
  double dotB = dot_product(n_dim, sisters_[1].v_, rhat);
  double thetaA = safe_acos(dotA);
  double thetaB = safe_acos(dotB);
  double sinthetaA = sin(thetaA);
  double sinthetaB = sin(thetaB);

  // Linear factor
  double linearfactor = -k * (rmag - r0);

  // theta A
  double thetaAoversinA = 0.0;
  double rhatcrossuA[3] = {0.0};
  double rhatcrossrhatcrossuA[3] = {0.0};
  if (sinthetaA != 0.0) {
    thetaAoversinA = thetaA / sinthetaA;
  }
  // cross_product(rhat, sisters_[0].orientation_, rhatcrossuA,
  cross_product(rhat, sisters_[0].v_, rhatcrossuA,
                n_dim); // cross product takes pointers, so can't just use again
  cross_product(rhat, rhatcrossuA, rhatcrossrhatcrossuA, n_dim);

  double thetaAfactor = 0.0;
  if (rmag > 0.0) {
    thetaAfactor = -kth / rmag * thetaAoversinA;
  }

  // theta B
  double thetaBoversinB = 0.0;
  double rhatcrossuB[3] = {0.0};
  double rhatcrossrhatcrossuB[3] = {0.0};
  if (sinthetaB != 0.0) {
    thetaBoversinB = thetaB / sinthetaB;
  }
  // cross_product(rhat, sisters_[1].orientation_, rhatcrossuB,
  cross_product(rhat, sisters_[1].v_, rhatcrossuB,
                n_dim); // cross product takes pointers, so can't just use again
  cross_product(rhat, rhatcrossuB, rhatcrossrhatcrossuB, n_dim);

  double thetaBfactor = 0.0;
  if (rmag > 0.0) {
    thetaBfactor = -kth / rmag * thetaBoversinB;
  }

  // update the forces
  double force[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < n_dim; ++i) {
    force_[i] =
        (linearfactor * rhat[i] + thetaAfactor * rhatcrossrhatcrossuA[i] +
         thetaBfactor * rhatcrossrhatcrossuB[i]);
    // f_chromosome[2 * ic + 1][i] -=
    //     (linearfactor * rhat[i] + thetaAfactor * rhatcrossrhatcrossuA[i] +
    //      thetaBfactor * rhatcrossrhatcrossuB[i]);
    // f_interkc[2 * ic][i] +=
    //     (linearfactor * rhat[i] + thetaAfactor * rhatcrossrhatcrossuA[i] +
    //      thetaBfactor * rhatcrossrhatcrossuB[i]);
    // f_interkc[2 * ic + 1][i] -=
    //     (linearfactor * rhat[i] + thetaAfactor * rhatcrossrhatcrossuA[i] +
    //      thetaBfactor * rhatcrossrhatcrossuB[i]);
  }
  sisters_[0].AddForce(force);
  sisters_[1].SubForce(force);

  // Other stuff required for torque
  double dotV =
      // dot_product(n_dim, sisters_[0].v_, sisters_[1].v_);
      dot_product(n_dim, sisters_[0].orientation_, sisters_[1].orientation_);
  double thetaV = safe_acos(dotV);
  double sinthetaV = sin(thetaV);
  double thetaVoversinV = 0.0;
  if (sinthetaV != 0.0) {
    thetaVoversinV = thetaV / sinthetaV;
  }
  double vhatAcrossvhatB[3] = {0.0};
  // cross_product(sisters_[0].v_, sisters_[1].v_,
  cross_product(sisters_[0].orientation_, sisters_[1].orientation_,
                vhatAcrossvhatB, n_dim);

  // torques
  double torque_one[3];
  double torque_two[3];
  for (int i = 0; i < n_dim; ++i) {
    torque_one[i] = (-kth * thetaAoversinA * rhatcrossuA[i] +
                     kv * thetaVoversinV * vhatAcrossvhatB[i]);
    torque_two[i] = (-kth * thetaBoversinB * rhatcrossuB[i] -
                     kv * thetaVoversinV * vhatAcrossvhatB[i]);
  }
  sisters_[0].AddTorque(torque_one);
  sisters_[1].AddTorque(torque_two);

  // Next: integrate to update their position (and add thermal noise)
  for (auto &&sis : sisters_) {
    sis.Integrate();
  }
  // Final: update position of centromere to be exactly in between two sisters
  double r_com[3];
  for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
    r_com[i_dim] =
        (sisters_[0].position_[i_dim] + sisters_[1].position_[i_dim]) / 2.0;
  }
  SetPosition(r_com);
}