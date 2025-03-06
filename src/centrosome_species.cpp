#include "cglass/centrosome_species.hpp"
#include "cglass/crosslink_species.hpp"
#include "cglass/filament_species.hpp"
#include "cglass/rigid_filament_species.hpp"

// TODO: add overlap detection, wca potential between centros

CentrosomeSpecies::CentrosomeSpecies(unsigned long seed) : Species(seed) {
  SetSID(species_id::centrosome);
}

void CentrosomeSpecies::Init(std::string spec_name, ParamsParser &parser) {
  Species::Init(spec_name, parser);
  if (GetNInsert() <= 0) {
    return;
  }
  // NAB uses 'anchor' to refer to SPBs. Anchors are unique objects in C-GLASS tho (FIXME)
  int n_anchors{GetNInsert()};
  // All of the new rotation information for this class
  q_ = new arma::vec4[n_anchors];
  A_ = new arma::mat33[n_anchors];
  mu_tb_ = new arma::mat33[n_anchors];
  mu_rb_ = new arma::mat33[n_anchors];
  sqrt_mu_tb_ = new arma::mat33[n_anchors];
  sqrt_mu_rb_ = new arma::mat33[n_anchors];
  bsub_ = new arma::mat::fixed<4, 3>[n_anchors];

  double attach_diameter[n_anchors];
  // initialize geometry stuff lol
  for (int i_spb{0}; i_spb < n_anchors; i_spb++) {
    double spb_diffusion = sparams_.diffusion;
    attach_diameter[i_spb] = sparams_.attach_diameter;
    // Now that we have the diffusion scale, we need to encode the proper SPB
    // translation and rotation matrices into the matrices (including sqrt)
    // for proper rigid body dynamics
    arma::mat33 translation;
    arma::mat33 rotation;
    translation.eye();
    rotation.eye();
    translation(0, 0) = 0.5 * spb_diffusion;
    translation(1, 1) = spb_diffusion;
    translation(2, 2) = spb_diffusion;
    rotation(0, 0) = spb_diffusion / SQR(0.5 * attach_diameter[i_spb]);
    rotation(1, 1) = 0.01 * spb_diffusion / SQR(0.5 * attach_diameter[i_spb]);
    rotation(2, 2) = 0.01 * spb_diffusion / SQR(0.5 * attach_diameter[i_spb]);
    // Set the mobility matrices (assuming they don't change) and the sqrt of same
    mu_tb_[i_spb] = translation;
    mu_rb_[i_spb] = rotation;
    sqrt_mu_tb_[i_spb] = arma::sqrtmat_sympd(translation);
    sqrt_mu_rb_[i_spb] = arma::sqrtmat_sympd(rotation);
  }
}

void CentrosomeSpecies::Reserve() {
  SpeciesBase::Reserve();
  members_.reserve(GetNInsert());
}

void CentrosomeSpecies::PopMember() { Species::PopMember(); }
void CentrosomeSpecies::AddMember() {
  Species::AddMember();
  int i_spb{(int)members_.size() - 1};
  double u[3] = {members_.back().GetU(0), members_.back().GetU(1),
                 members_.back().GetU(2)};
  // Create the rotation matrix from body-space coordinate transformation
  CreateRotationMatrixBodySpace(A_[i_spb], u, members_.back().v_,
                                members_.back().w_);
  // Create the associated quaternion describing initial orientation of the rigid body
  QuaternionFromRotationMatrix(q_[i_spb], A_[i_spb]);
}

void CentrosomeSpecies::AnchorFilaments(SpeciesBase *filas, bool flexible) {

  int n_filaments{filas->GetNMembers()};
  if (n_filaments != GetNMembers() * sparams_.num_anchors_ea) {
    printf("Error: mis-match between SPB anchor and filament number. \n");
    exit(1);
  }

  // SF TODO: generalize this for more than 1 SPB / anchor lol
  size_t i_spb{0};
  for (int i_fila{0}; i_fila < filas->GetNMembers(); i_fila++) {
    Mesh *fil{dynamic_cast<Mesh *>(filas->GetMember(i_fila))};
    // SF TODO baaaaaaaaaaaaad
    int i_anchor{i_fila < sparams_.num_anchors_ea
                     ? i_fila
                     : i_fila - sparams_.num_anchors_ea};
    if (i_fila == sparams_.num_anchors_ea) {
      i_spb++;
    }
    // FIXME for discrete anchor site locations
    const double *const anchor_u = members_[i_spb].anchors_[i_anchor].u_;
    const double *const anchor_pos = members_[i_spb].anchors_[i_anchor].pos_;
    double anchor_endpoint[3] = {0, 0, 0};
    double new_pos[3] = {0, 0, 0};
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      anchor_endpoint[i_dim] =
          anchor_pos[i_dim] +
          anchor_u[i_dim] * members_[i_spb].anchors_[i_anchor].r0_;
      printf("END: %g = %g + %g*%g\n", anchor_endpoint[i_dim],
             anchor_pos[i_dim], anchor_u[i_dim],
             members_[i_spb].anchors_[i_anchor].r0_);
      new_pos[i_dim] =
          anchor_endpoint[i_dim] + 0.5 * fil->GetLength() * anchor_u[i_dim];
      printf("%g = %g + %g\n", new_pos[i_dim], anchor_pos[i_dim],
             0.5 * fil->GetLength() * anchor_u[i_dim]);
    }
    if (flexible) {
      dynamic_cast<Filament *>(fil)->InsertAt(new_pos, anchor_u);
    } else {
      dynamic_cast<RigidFilament *>(fil)->InsertAt(new_pos, anchor_u);
      printf("DING\n");
    }
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      printf("(%g)\n", fil->GetTailPosition()[i_dim]);
    }
    if (flexible) {
      dynamic_cast<Filament *>(fil)->UpdatePosition();
    } else {
      dynamic_cast<RigidFilament *>(fil)->UpdatePosition();
      printf("DONG\n");
    }
    members_[i_spb].anchors_[i_anchor].filament_ = fil;
  }

  // int i_fila{0};
  // for (auto &&fil : filas->GetMembers()) {
  //   double u[3] = {0, 0, 1}; // Initially align with z
  //   fil.SetOrientation(u);
  //   printf("FILA #%i\n", i_fila++);
  // }

  // UPDATE POS AND U OF ANCHORSITE
  // UPDATE POS AND U OF FILAMENT
  // UPDATE REL_POS AND REL_U OF ANCHORSITE
}

void CentrosomeSpecies::UpdatePositions() {

  double kT = 1.0;
  const double four_epsilon = 4.0;
  const double r_cutoff_factor2 = pow(2, 1.0 / 3.0);
  const double thermal_diff_scale = 1;
  double sys_radius = space_->radius;
  double delta_t = params_->delta;
  const arma::vec3 nx = {1.0, 0.0, 0.0};
  const arma::vec3 ny = {0.0, 1.0, 0.0};
  const arma::vec3 nz = {0.0, 0.0, 1.0};
  for (int idx{0}; idx < members_.size(); idx++) {
    Centrosome *centro{&members_[idx]};
    // First, WCA potentials
    for (int jdx{idx + 1}; jdx < members_.size(); jdx++) {
      double dr_mag2 = 0.0;
      double dr[3];
      for (int i = 0; i < params_->n_dim; ++i) {
        dr[i] = members_[jdx].GetR(i) - members_[idx].GetR(i);
        dr_mag2 += SQR(dr[i]);
      }
      const double d_eff = 1.0;
      const double d_eff2 = SQR(d_eff);
      double dr_mag = sqrt(dr_mag2);
      double scale =
          (0.5 * (members_[idx].diameter_ + members_[jdx].diameter_) - d_eff) /
          dr_mag;
      dr_mag2 = 0.0;
      for (int i = 0; i < params_->n_dim; ++i) {
        dr[i] -= scale * dr[i];
        dr_mag2 += SQR(dr[i]);
      }

      if (dr_mag2 < r_cutoff_factor2 * d_eff2) {
        // Add contributions to potential energy and forces.
        double rho2 = d_eff2 / dr_mag2;
        double rho6 = CUBE(rho2);
        double rho12 = SQR(rho6);

        double factor = 6.0 * four_epsilon * (2.0 * rho12 - rho6) / dr_mag2;
        double f_wca[3];
        for (int i = 0; i < params_->n_dim; ++i) {
          f_wca[i] = factor * dr[i];
        }
        members_[idx].SubForce(f_wca);
        members_[jdx].AddForce(f_wca);
      }
    }
    // Next, handle deterministic forces/torques and diffusion
    // Construct the amount (and direction) the SPB is outside the radius
    // If outside the nucleus, negative
    double rmag = 0.0;
    for (int i = 0; i < 3; ++i) {
      rmag += SQR(centro->GetR(i));
    }
    rmag = std::sqrt(rmag);
    double rhat[3] = {0.0};
    for (int i = 0; i < 3; ++i) {
      rhat[i] = centro->GetR(i) / rmag;
    }
    // Construct the force vector
    double forcevec[3] = {0.0};
    // SF TODO implement alternate force versions
    /*
    // If we are using a harmonic force
    if (properties->anchors.spb_confinement_type == 0) {
      double k = properties->anchors.centrosome_confinement_radial_k_;

      // Guard against pathological behavior by capping the force at some level
      double factor = -k * (rmag - sys_radius);
      double fcutoff = 0.1 / parameters->delta /
                       MAX(properties->anchors.mu_tb_[idx](0, 0),
                           properties->anchors.mu_tb_[idx](1, 1));
      if (factor < 0.0) {
        if (-1.0 * factor > fcutoff) {
          factor = -1.0 * fcutoff;
          std::cerr << " *** Force exceeded fcutoff "
                       "spb_confinement_potential_harmonic_bd radial ***\n";
        }
      } else {
        if (factor > fcutoff) {
          factor = fcutoff;
          std::cerr << " *** Force exceeded fcutoff "
                       "spb_confinement_potential_harmonic_bd radial ***\n";
        }
      }

      for (int i = 0; i < 3; ++i) {
        forcevec[i] = factor * rhat[i];
      }
    } else if (properties->anchors.spb_confinement_type == 1) {
        */

    // properties->anchors.centrosome_confinement_f0_;
    double f0 = sparams_.wall_f0;
    // double f0 = 10.3406326034025;
    double delta_r = ABS(sys_radius - rmag);
    double ne_ratio = sparams_.wall_ne;
    // double ne_ratio{24.61538};
    double factor = CalcNonMonotonicWallForce(ne_ratio, f0, delta_r);
    // Check the sign of the force, want to move outwards if we are in the nucleoplasm
    if (rmag > sys_radius) {
      for (int i = 0; i < 3; ++i) {
        forcevec[i] = -1.0 * factor * rhat[i];
      }
    } else {
      for (int i = 0; i < 3; ++i) {
        forcevec[i] = 1.0 * factor * rhat[i];
      }
    }
    // }
    //Based on the potential, calculate the torques on the system
    double uhat_dot_rhat = dot_product(3, rhat, centro->GetOrientation());
    double uhat_cross_rhat[3] = {0.0};
    cross_product(centro->GetOrientation(), rhat, uhat_cross_rhat, 3);
    double kr = sparams_.wall_kr;
    // double kr = 999.0; //properties->anchors.centrosome_confinement_angular_k_;
    // The torque is thusly
    double torquevec[3] = {0.0};
    for (int i = 0; i < 3; ++i) {
      torquevec[i] = -kr * (uhat_dot_rhat + 1.0) * uhat_cross_rhat[i];
      // printf("torque[%i] = %g\n", i, torquevec[i]);
    }
    // Add the contribution to the total forces
    centro->AddForce(forcevec);
    centro->AddTorque(torquevec);

    // SF TODO set up case (and choice) for confined SPBs
    // BELOW CODE IS FOR 'FREE' SPB, I.E., NOT CONFINED TO MEMBRANE
    // Set up armadillo versions of all of our information of interest
    arma::vec3 r = {centro->GetR(0), centro->GetR(1), centro->GetR(2)};
    arma::vec3 u = {centro->GetU(0), centro->GetU(1), centro->GetU(2)};
    arma::vec3 v = {centro->GetV(0), centro->GetV(1), centro->GetV(2)};
    arma::vec3 w = {centro->GetW(0), centro->GetW(1), centro->GetW(2)};
    arma::vec3 force = {centro->GetForce()[0], centro->GetForce()[1],
                        centro->GetForce()[2]};
    arma::vec3 torque = {centro->GetTorque()[0], centro->GetTorque()[1],
                         centro->GetTorque()[2]};
    // Deterministic translation
    arma::vec3 dr = (((A_[idx] * mu_tb_[idx]) * A_[idx].t()) * force) * delta_t;
    // Random thermal translation motion
    arma::vec3 theta = {thermal_diff_scale * rng_.RandomNormal(1.0),
                        thermal_diff_scale * rng_.RandomNormal(1.0),
                        thermal_diff_scale * rng_.RandomNormal(1.0)};
    dr +=
        ((A_[idx] * sqrt_mu_tb_[idx]) * theta) * std::sqrt(2.0 * kT * delta_t);
    // Now handle the rotation part
    bsub_[idx](0, 0) = -q_[idx][1];
    bsub_[idx](0, 1) = -q_[idx][2];
    bsub_[idx](0, 2) = -q_[idx][3];
    bsub_[idx](1, 0) = q_[idx][0];
    bsub_[idx](1, 1) = -q_[idx][3];
    bsub_[idx](1, 2) = q_[idx][2];
    bsub_[idx](2, 0) = q_[idx][3];
    bsub_[idx](2, 1) = q_[idx][0];
    bsub_[idx](2, 2) = -q_[idx][1];
    bsub_[idx](3, 0) = -q_[idx][2];
    bsub_[idx](3, 1) = q_[idx][1];
    bsub_[idx](3, 2) = q_[idx][0];
    // Deterministic rotation part
    arma::vec4 dq =
        (((bsub_[idx] * mu_rb_[idx]) * A_[idx].t()) * torque) * delta_t;
    // Random reorientation
    arma::vec3 theta_rotation = {thermal_diff_scale * rng_.RandomNormal(1.0),
                                 thermal_diff_scale * rng_.RandomNormal(1.0),
                                 thermal_diff_scale * rng_.RandomNormal(1.0)};
    dq += ((bsub_[idx] * sqrt_mu_rb_[idx]) * theta_rotation) *
          std::sqrt(2.0 * kT * delta_t);
    // Compute the lagrange correction
    arma::vec4 qtilde = q_[idx] + dq;
    // Compute the corrected quaternion
    double lambdaq = ComputeLagrangeCorrection(q_[idx], qtilde);
    if (std::isnan(lambdaq)) {
      std::cerr << " quaternion correction error: " << lambdaq << std::endl;
    }
    // Actually update the information
    r += dr;
    q_[idx] = qtilde + lambdaq * q_[idx];
    // Update the rotation matrix and orientation vectors
    RotationMatrixFromQuaternion(A_[idx], q_[idx]);
    u = A_[idx] * nx;
    v = A_[idx] * ny;
    w = A_[idx] * nz;

    double r_new[3];
    double u_new[3];
    double v_new[3];
    double w_new[3];
    // Write back into the anchor structure the proper double variables
    for (int i = 0; i < 3; ++i) {
      r_new[i] = r(i);
      u_new[i] = u(i);
      v_new[i] = v(i);
      w_new[i] = w(i);
    }
    centro->UpdatePosition(r_new, u_new, v_new, w_new);
  }
}

double CentrosomeSpecies::CalcNonMonotonicWallForce(double ne_ratio, double f0,
                                                    double delta_r) {
  //Constants from paper
  double a1 = 0.746;
  double a2 = 0.726;
  double alpha1 = 0.347;
  double alpha2 = 3.691;
  double a = a1 * a2;
  double b = sqrt(ne_ratio);
  double c = alpha1 + alpha2;
  // Calculate location of maximum force
  double xm = b * (((7.0 / 4.0) * M_PI) - c); //value of dr that gives max force
  // Calculate wall force
  // Exponential prefactor added to make forces conform to boundary conditions.
  // Characteristic length chosen to be 1/30th the value of delta_r that gives max force
  double factor =
      (1 - exp(-delta_r * 30 / xm)) *
      (2.0 * f0 * a * exp(-1.0 * delta_r / b) * cos((delta_r / b) + c) + f0);
  return factor;
}

double CentrosomeSpecies::ComputeLagrangeCorrection(const arma::vec4 &q,
                                                    const arma::vec4 &qtilde) {
  double qqtilde = arma::dot(q, qtilde);
  double qtilde2 = arma::dot(qtilde, qtilde);
  // Solve the quadratic formulat for what is correct. Doesn't matter which
  // solution we take, as the rotations are equivalent
  double d = (2.0 * qqtilde) * (2.0 * qqtilde) - 4.0 * (qtilde2 - 1.0);
  double lambdaq = (-2.0 * qqtilde + std::sqrt(d)) / 2.0;
  return lambdaq;
}

void CentrosomeSpecies::CreateRotationMatrixBodySpace(arma::mat33 &A, double *u,
                                                      double *v, double *w) {
  arma::vec3 u_ = {u[0], u[1], u[2]};
  arma::vec3 v_ = {v[0], v[1], v[2]};
  arma::vec3 w_ = {w[0], w[1], w[2]};
  A.col(0) = u_;
  A.col(1) = v_;
  A.col(2) = w_;
}

void CentrosomeSpecies::QuaternionFromRotationMatrix(arma::vec4 &q,
                                                     const arma::mat33 R) {
  if ((R(1, 1) >= -R(2, 2)) && (R(0, 0) >= -R(1, 1)) && (R(0, 0) >= -R(2, 2))) {
    q[0] = std::sqrt(1 + R(0, 0) + R(1, 1) + R(2, 2));
    q[1] = (R(2, 1) - R(1, 2)) / std::sqrt(1 + R(0, 0) + R(1, 1) + R(2, 2));
    q[2] = (R(0, 2) - R(2, 0)) / std::sqrt(1 + R(0, 0) + R(1, 1) + R(2, 2));
    q[3] = (R(1, 0) - R(0, 1)) / std::sqrt(1 + R(0, 0) + R(1, 1) + R(2, 2));
  } else if ((R(1, 1) <= -R(2, 2)) && (R(0, 0) >= R(1, 1)) &&
             (R(0, 0) >= R(2, 2))) {
    q[0] = (R(2, 1) - R(1, 2)) / std::sqrt(1 + R(0, 0) - R(1, 1) - R(2, 2));
    q[1] = std::sqrt(1 + R(0, 0) - R(1, 1) - R(2, 2));
    q[2] = (R(1, 0) + R(0, 1)) / std::sqrt(1 + R(0, 0) - R(1, 1) - R(2, 2));
    q[3] = (R(2, 0) + R(0, 2)) / std::sqrt(1 + R(0, 0) - R(1, 1) - R(2, 2));
  } else if ((R(1, 1) >= R(2, 2)) && (R(0, 0) <= R(1, 1)) &&
             (R(0, 0) <= -R(2, 2))) {
    q[0] = (R(0, 2) - R(2, 0)) / std::sqrt(1 - R(0, 0) + R(1, 1) - R(2, 2));
    q[1] = (R(1, 0) + R(0, 1)) / std::sqrt(1 - R(0, 0) + R(1, 1) - R(2, 2));
    q[2] = std::sqrt(1 - R(0, 0) + R(1, 1) - R(2, 2));
    q[0] = (R(2, 1) + R(1, 2)) / std::sqrt(1 - R(0, 0) + R(1, 1) - R(2, 2));
  } else if ((R(1, 1) <= R(2, 2)) && (R(0, 0) <= -R(1, 1)) &&
             (R(0, 0) <= R(2, 2))) {
    q[0] = (R(1, 0) - R(0, 1)) / std::sqrt(1 - R(0, 0) - R(1, 1) + R(2, 2));
    q[1] = (R(2, 0) + R(0, 2)) / std::sqrt(1 - R(0, 0) - R(1, 1) + R(2, 2));
    q[2] = (R(2, 1) + R(1, 2)) / std::sqrt(1 - R(0, 0) - R(1, 1) + R(2, 2));
    q[3] = std::sqrt(1 - R(0, 0) - R(1, 1) + R(2, 2));
  } else {
    std::cerr
        << "ERROR: quaterion_from_rotation_matrix did not work correctly!\n";
    exit(1);
  }

  // Have to multiply by 1/2 on all cases
  for (auto i = 0; i < 4; ++i) {
    q[i] = 0.5 * q[i];
  }
}

void CentrosomeSpecies::RotationMatrixFromQuaternion(arma::mat33 &R,
                                                     const arma::vec4 &q) {
  R(0, 0) = 1. - 2. * (q[2] * q[2] + q[3] * q[3]);
  R(0, 1) = 2. * (q[1] * q[2] - q[3] * q[0]);
  R(0, 2) = 2. * (q[1] * q[3] + q[2] * q[0]);
  R(1, 0) = 2. * (q[1] * q[2] + q[3] * q[0]);
  R(1, 1) = 1 - 2. * (q[1] * q[1] + q[3] * q[3]);
  R(1, 2) = 2. * (q[2] * q[3] - q[1] * q[0]);
  R(2, 0) = 2. * (q[1] * q[3] - q[2] * q[0]);
  R(2, 1) = 2. * (q[2] * q[3] + q[1] * q[0]);
  R(2, 2) = 1 - 2. * (q[1] * q[1] + q[2] * q[2]);
}

void CentrosomeSpecies::GetInteractors(std::vector<Object *> &ixors) {
  for (auto &&centro : members_) {
    centro.GetInteractors(ixors);
  }
}
