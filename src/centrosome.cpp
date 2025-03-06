#include "cglass/centrosome.hpp"
#include "cglass/filament.hpp"
#include <iostream>

size_t Centrosome::i_spb_ = 0;

Centrosome::Centrosome(unsigned long seed) : Object(seed) {
  SetSID(species_id::centrosome);
  printf("NEW centrosome\n");
}

// called during AddMember()
void Centrosome::Init(centrosome_parameters *sparams) {

  int n_dim{params_->n_dim};
  sparams_ = sparams;
  color_ = sparams->color;
  draw_ = draw_type::_from_string(sparams_->draw_type.c_str());
  length_ = sparams->length;
  diameter_ = sparams->diameter;
  diffusion_ = 0.0;     //noise_tr_ * sqrt(24.0 * diameter_ / delta_);
  diffusion_rot_ = 0.0; //noise_rot_ * sqrt(8.0 * CUBE(diameter_) / delta_);
  Logger::Trace("Inserting object %d randomly", GetOID());
  // SF FIXME this is bad lol
  printf("I_SPB: %zu\n", i_spb_);
  index_ = i_spb_++;

  double sign{index_ == 0 ? 1.0 : -1.0};
  double pos[3] = {0, params_->system_radius, 0}; // position of C.O.M.
  double u[3] = {0, -1.0, 0}; // normal; points to boundary center
  if (index_ > 0) {
    if (sparams_->insertion_type.compare("bioriented") == 0) {
      printf("BIORIENT\n");
      if (sparams_->num != 2) {
        printf("Error; need 2 SPBS to biorient\n");
        exit(1);
      }
      pos[1] *= -1.0;
      u[1] *= -1.0;
    } else if (sparams_->insertion_type.compare("normal") == 0) {
      printf("NORMAL\n");
      // Displace second SPB by 1.5 diameter to avoid overlap
      pos[0] = -1.5 * sparams->diameter;
      double r_mag2{0.0};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        r_mag2 += SQR(pos[i_dim]);
      }
      r_mag2 = sqrt(r_mag2);
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        u[i_dim] = -1.0 * pos[i_dim] / r_mag2;
        pos[i_dim] = -1.0 * u[i_dim] * params_->system_radius;
      }
    }
  }

  // Insert centrosome
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
  for (int i_dim{0}; i_dim < n_dim; i_dim++) {
    v[i_dim] *= norm_factor;
  }
  cross_product(u, v, w, 3);
  InsertAt(pos, u);
  // Update auxiliary vectors
  for (int i{0}; i < 3; i++) {
    r_[i] = position_[i] = pos[i];
    u_[i] = orientation_[i] = u[i];
    v_[i] = v[i];
    w_[i] = w[i];
  }

  // Initialize associated anchors
  n_filaments_ = sparams->num_anchors_ea; // FIXME: change to n_anchors
  // Initialize AnchorSite properties
  anchors_.resize(n_filaments_);
  for (int i_anchor{0}; i_anchor < n_filaments_; i_anchor++) {
    // Set basic physical properties -- SF TODO put in YAML files
    anchors_[i_anchor].k_ = 100;
    anchors_[i_anchor].kr_ = 1000;
    anchors_[i_anchor].r0_ = 1;
    // Centrosome surface is 2-D,so we can define position w/ R and phi
    // For now, distribute them uniformly along circumference of a circle
    // SF TODO make this adjustable, or even fully just randomly distributed
    double r_on_spb = 0.333 * sparams_->diameter;
    double delta_phi = 2.0 * M_PI / double(n_filaments_);
    double offset_phi{index_ == 0 ? 0 : delta_phi / 2.0};
    double base_pos[3];
    for (int i_dim{0}; i_dim < n_dim; i_dim++) {
      base_pos[i_dim] =
          GetPosition()[i_dim] +
          r_on_spb * sin(offset_phi + i_anchor * delta_phi) * v_[i_dim] +
          r_on_spb * cos(offset_phi + i_anchor * delta_phi) * w_[i_dim];
    }
    // Set anchor absolute position and orientation (same as SPB for now)
    // These are dynamic and are continuously updated
    for (int i_dim{0}; i_dim < n_dim; i_dim++) {
      anchors_[i_anchor].pos_[i_dim] = base_pos[i_dim];
      anchors_[i_anchor].u_[i_dim] = GetOrientation()[i_dim];
    }
    // Set anchor relative position and orientation
    // These are static and are relative to the local coord sys of SPBs
    double r_rel[3];
    for (int i_dim{0}; i_dim < n_dim; i_dim++) {
      r_rel[i_dim] = anchors_[i_anchor].pos_[i_dim] - GetPosition()[i_dim];
      printf("R_REL[%i] = %g\n", i_dim, r_rel[i_dim]);
    }
    double u_proj = dot_product(n_dim, r_rel, u_);
    double v_proj = dot_product(n_dim, r_rel, v_);
    double w_proj = dot_product(n_dim, r_rel, w_);
    for (int i_dim{0}; i_dim < n_dim; i_dim++) {
      anchors_[i_anchor].pos_rel_[i_dim] =
          u_[i_dim] * u_proj + v_[i_dim] * v_proj + w_[i_dim] * w_proj;
      printf("%g\n", anchors_[i_anchor].pos_rel_[i_dim]);
    }
    // no angular springs for now (mode 0 from NAB)
    double u_anchor[3];
    for (int i_dim{0}; i_dim < n_dim; i_dim++) {
      u_anchor[i_dim] = anchors_[i_anchor].u_[i_dim];
    }
    double u_proj_rel = dot_product(n_dim, u_anchor, u_);
    double v_proj_rel = dot_product(n_dim, u_anchor, v_);
    double w_proj_rel = dot_product(n_dim, u_anchor, w_);
    anchors_[i_anchor].u_rel_[0] = u_proj_rel;
    anchors_[i_anchor].u_rel_[1] = v_proj_rel;
    anchors_[i_anchor].u_rel_[2] = w_proj_rel;
    // Position and orientation set; still need to anchor w/ filament
  }
  // UpdatePosition(pos, u, v, w);
}

void Centrosome::Draw(std::vector<graph_struct *> &graph_array) {
  for (auto &&anch : anchors_) {
    if (anch.filament_ == nullptr) {
      continue;
    }
    double u[3] = {0.0, 0.0, 0.0};
    double r_teth{0.0};
    for (int i{0}; i < n_dim_; ++i) {
      r_teth += SQR(anch.pos_[i] - anch.filament_->GetTailPosition()[i]);
    }
    r_teth = sqrt(r_teth);
    for (int i{0}; i < n_dim_; ++i) {
      // u[i] = (anch.filament_->GetTailPosition()[i] - anch.pos_[i]) / r_teth;
      anch.g_.r[i] = anch.pos_[i] + 0.5 * anch.u_[i] * r_teth;
    }
    std::copy(anch.u_, anch.u_ + 3, anch.g_.u);
    anch.g_.color = 1.8 * M_PI;
    anch.g_.diameter = 0.3;
    anch.g_.length = r_teth;
    anch.g_.draw = draw_;
    graph_array.push_back(&anch.g_);
  }
}

void Centrosome::ApplyInteractions() {
  int i_anchor{0};
  for (auto &&anchor : anchors_) {
    double dr[3];
    // Now that anchor positions are updated, apply force to filaments
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      dr[i_dim] =
          anchor.pos_[i_dim] - anchor.filament_->GetTailPosition()[i_dim];
    }
    double dr_mag2 = dot_product(params_->n_dim, dr, dr);
    double factor{dr_mag2 > 0.0 ? anchor.k_ * (1.0 - anchor.r0_ / sqrt(dr_mag2))
                                : 0.0};
    double f_spring[3];
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      f_spring[i_dim] = factor * dr[i_dim];
    }
    // apply forces -- SF TODO: integrate this into Interacte() routine
    anchor.filament_->AddForceTail(f_spring);
    SubForce(f_spring);
    double r_spb[3];
    double r_bond[3];
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      r_spb[i_dim] = anchor.pos_[i_dim] - GetPosition()[i_dim];
      // SF TODO verify below calculates center of tail rod position
      r_bond[i_dim] = anchor.pos_[i_dim] -
                      (anchor.filament_->GetTailPosition()[i_dim] +
                       0.5 * anchor.filament_->GetTailOrientation()[i_dim] *
                           anchor.filament_->GetBondLength());
      // printf("r_spb[%i] = %g\n", i_dim, r_spb[i_dim]);
      // printf("r_bond[%i] = %g\n", i_dim, r_bond[i_dim]);
    }
    double tau_spb[3];
    double tau_bond[3];
    cross_product(r_spb, f_spring, tau_spb, params_->n_dim);
    cross_product(r_bond, f_spring, tau_bond, params_->n_dim);
    anchor.filament_->AddTorqueTail(tau_bond);
    SubTorque(tau_spb);
  }
}

void Centrosome::UpdatePosition(double *r_new, double *u_new, double *v_new,
                                double *w_new) {

  for (int i{0}; i < 3; i++) {
    r_[i] = position_[i] = r_new[i];
    u_[i] = orientation_[i] = u_new[i];
    v_[i] = v_new[i];
    w_[i] = w_new[i];
  }

  int i_anchor{0};
  for (auto &&anchor : anchors_) {
    /*
    double factor_u = dot_product(3, anchor.pos_rel_, u_);
    double factor_v = dot_product(3, anchor.pos_rel_, v_);
    double factor_w = dot_product(3, anchor.pos_rel_, w_);
    // Calculate new lab frame coordinate
    for (int i = 0; i < 3; ++i) {
      anchor.pos_[i] =
          r_[i] + factor_u * u_[i] + factor_v * v_[i] + factor_w * w_[i];

      // printf("%.2f = %.2f + %.2f*%.2f + %.2f*%.2f + %.2f*%.2f\n",
      //        anchor.pos_[i], r_[i], factor_u, u_[i], factor_v, v_[i], factor_w,
      //        w_[i]);
    }
    // Make sure that we properly attack the orientation of the anchor list
    // Orientations!
    // double original_u[3] = {0.0};
    for (int i = 0; i < 3; ++i) {
      // original_u[i] = anchor.u_[i];
      anchor.u_[i] = anchor.u_rel_[0] * u_[i] + anchor.u_rel_[1] * v_[i] +
                     anchor.u_rel_[2] * w_[i];
    }
    */
    // SF: trash (or NOT??) below from misinterpreting NAB algorithm

    // Centrosome surface is 2-D,so we can define position w/ R and phi
    // For now, distribute them uniformly along circumference of a circle
    double r_on_spb = 0.333 * sparams_->diameter;
    double delta_phi = 2.0 * M_PI / double(n_filaments_);
    double offset_phi{index_ == 0 ? 0 : delta_phi / 2.0};
    // printf("r = %g\n", r_on_spb);
    double base_pos[3];
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      base_pos[i_dim] =
          GetPosition()[i_dim] +
          r_on_spb * sin(offset_phi + i_anchor * delta_phi) * v_[i_dim] +
          r_on_spb * cos(offset_phi + i_anchor * delta_phi) * w_[i_dim];
    }
    i_anchor++;
    if (anchor.filament_ == nullptr) {
      printf("hi from anchor #%i\n", i_anchor - 1);
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        anchor.pos_[i_dim] = base_pos[i_dim];
        anchor.u_[i_dim] = GetOrientation()[i_dim];
      }
      continue;
    }
    double r_new[3] = {0, 0, 0};
    double r_mag{0.0};
    // calculate u first; from filament tip to anchor SITE
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      r_new[i_dim] =
          base_pos[i_dim] - anchor.filament_->GetTailPosition()[i_dim];
      r_mag += SQR(r_new[i_dim]);
    }
    r_mag = sqrt(r_mag);
    // printf("r_teth = %g\n", r_teth);
    // double u_new[3] = {0, 0, 0};
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      anchor.u_[i_dim] = -r_new[i_dim] / r_mag;
      anchor.pos_[i_dim] = base_pos[i_dim];
    }

    // Create a check if the position isn't working correctly
    for (int i = 0; i < 3; ++i) {
      if (std::isnan(anchor.pos_[i])) {
        printf("woah buddy\n");
        exit(1);
        /*
          // Print out the information of the centrosome and exit
          std::cerr << "Encountered an error in SPB update code, printing "
                       "information then exiting\n";
          std::cerr << "  step: " << properties->i_current_step << std::endl;
          print_centrosomes_information(idx, &(properties->anchors));
          std::cerr << "  Force = " << force.as_row();
          std::cerr << "  Torque = " << torque.as_row();
          exit(1);
          */
      }
    }
    /*
      // Now that anchor positions are updated, apply force to filaments
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        dr[i_dim] = anchor.filament_->GetPosition()[i_dim] -
                    anchor.pos_[i_dim] -
                    (0.5 * anchor.filament_->GetLength() *
                     anchor.filament_->GetOrientation()[i_dim]);
      }
      double dr_mag2 = dot_product(params_->n_dim, dr, dr);
      double factor = anchor.k_ * (1.0 - anchor.r0_ / sqrt(dr_mag2));
      double f_spring[3];
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        f_spring[i_dim] = factor * dr[i_dim];
        // printf("f[%i] = %g\n", i_dim, f_spring[i_dim]);
      }
      // apply forces -- SF TODO: integrate this into Interacte() routine
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      }
      anchor.filament_->SubForce(f_spring);
      AddForce(f_spring);
      // anchor.filament_->UpdatePosition();
      // UpdatePosition();
      */
  }
  // UpdatePosition();
}