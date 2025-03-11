#include <cglass/filament.hpp>
#include <cglass/triangle_mesh_boundary.hpp>
#include <unistd.h>

void TriMeshBoundary::Init(system_parameters *params) {
  // SF temp before integrating with output_manager
  std::string force_filename{params->run_name + "_membrane_forces.file"};
  std::string vrt_filename{params->run_name + "_membrane_vrt_positions.file"};
  std::string adj_filename{params->run_name + "_membrane_vrt_adjacency.file"};
  forces_ = fopen(force_filename.c_str(), "w");
  vertices_ = fopen(vrt_filename.c_str(), "w");
  adjacency_ = fopen(adj_filename.c_str(), "w");
  if (forces_ == nullptr or vertices_ == nullptr or adjacency_ == nullptr) {
    printf("Error generating mesh data files\n");
    exit(1);
  }
  params_ = params;
  SetParameters();
  MakeIcosphere();
  InitializeMesh();
}

void TriMeshBoundary::SetParameters() {
  long seed{params_ == nullptr ? 0 : params_->seed};
  rng_ = new RNG(seed);
  r_sys_ = params_->system_radius;
  kappa_B_ = params_->mesh_kB;
  kappa_ = params_->mesh_k;
  kappa_l_ = params_->mesh_kl;
  kappa_v_ = params_->mesh_kV;
  gamma_ = params_->node_gamma;
  if (params_->draw_centroid) {
    o_.diameter = r_sys_ / 10.0;
    o_.length = 0.0;
    o_.draw = draw_type::fixed;
    o_.color = 2 * M_PI;
  }
}

void TriMeshBoundary::MakeIcosphere() {
  // We start with 20 faces, and each subdivision algorithm quadruples the number
  size_t n_faces{size_t(20 * std::pow(4, params_->n_subdivisions))};
  // For a triangle mesh, 2E = 3F
  size_t n_edges{3 * n_faces / 2};
  // For a genus 0 shape, V = 2 - F + E = 2 + F/2
  size_t n_verts{2 + n_faces / 2};
  tris_.reserve(n_faces);
  edges_.reserve(n_edges);
  vrts_.reserve(n_verts);
  printf("Expect %zu faces, %zu edges, %zu verts\n", n_faces, n_edges, n_verts);
  MakeIcosahedron(); // 20 triangle faces initially
  // 80 -> 320 -> 1,280 -> 5,120 -> 20,480 -> 81,920 -> 327,680 -> 1,310,720 faces
  for (int i_divide{0}; i_divide < params_->n_subdivisions; i_divide++) {
    printf("Subdivision iteration %i:\n", i_divide + 1);
    DivideFaces();
  }
  if (n_faces != tris_.size() or n_edges != edges_.size() or
      n_verts != vrts_.size()) {
    printf("Error in TriMesh::MakeIcosphere\n");
    exit(1);
  }
  ProjectToUnitSphere();
}

void TriMeshBoundary::MakeIcosahedron() {
  printf("Initializing icosahedron (20 triangles; 12 vertices)\n");
  double phi = (1.0f + sqrt(5.0f)) * 0.5f; // golden ratio
  double a = r_sys_ * 1.0f;
  double b = r_sys_ * 1.0f / phi;
  // add vertices
  vrts_.emplace_back(0.0, b, -a);
  vrts_.emplace_back(b, a, 0.0);
  vrts_.emplace_back(-b, a, 0.0);
  vrts_.emplace_back(0.0, b, a);
  vrts_.emplace_back(0.0, -b, a);
  vrts_.emplace_back(-a, 0.0, b);
  vrts_.emplace_back(0.0, -b, -a);
  vrts_.emplace_back(a, 0.0, -b);
  vrts_.emplace_back(a, 0.0, b);
  vrts_.emplace_back(-a, 0.0, -b);
  vrts_.emplace_back(b, -a, 0.0);
  vrts_.emplace_back(-b, -a, 0.0);
  // add triangles
  tris_.emplace_back(&vrts_[2], &vrts_[1], &vrts_[0]);
  tris_.emplace_back(&vrts_[1], &vrts_[2], &vrts_[3]);
  tris_.emplace_back(&vrts_[5], &vrts_[4], &vrts_[3]);
  tris_.emplace_back(&vrts_[4], &vrts_[8], &vrts_[3]);
  tris_.emplace_back(&vrts_[7], &vrts_[6], &vrts_[0]);
  tris_.emplace_back(&vrts_[6], &vrts_[9], &vrts_[0]);
  tris_.emplace_back(&vrts_[11], &vrts_[10], &vrts_[4]);
  tris_.emplace_back(&vrts_[10], &vrts_[11], &vrts_[6]);
  tris_.emplace_back(&vrts_[9], &vrts_[5], &vrts_[2]);
  tris_.emplace_back(&vrts_[5], &vrts_[9], &vrts_[11]);
  tris_.emplace_back(&vrts_[8], &vrts_[7], &vrts_[1]);
  tris_.emplace_back(&vrts_[7], &vrts_[8], &vrts_[10]);
  tris_.emplace_back(&vrts_[2], &vrts_[5], &vrts_[3]);
  tris_.emplace_back(&vrts_[8], &vrts_[1], &vrts_[3]);
  tris_.emplace_back(&vrts_[9], &vrts_[2], &vrts_[0]);
  tris_.emplace_back(&vrts_[1], &vrts_[7], &vrts_[0]);
  tris_.emplace_back(&vrts_[11], &vrts_[9], &vrts_[6]);
  tris_.emplace_back(&vrts_[7], &vrts_[10], &vrts_[6]);
  tris_.emplace_back(&vrts_[5], &vrts_[11], &vrts_[4]);
  tris_.emplace_back(&vrts_[10], &vrts_[8], &vrts_[4]);
  // add edges
  edges_.emplace_back(&vrts_[0], &vrts_[1]);
  edges_.emplace_back(&vrts_[0], &vrts_[2]);
  edges_.emplace_back(&vrts_[0], &vrts_[6]);
  edges_.emplace_back(&vrts_[0], &vrts_[7]);
  edges_.emplace_back(&vrts_[0], &vrts_[9]);
  edges_.emplace_back(&vrts_[1], &vrts_[2]);
  edges_.emplace_back(&vrts_[1], &vrts_[3]);
  edges_.emplace_back(&vrts_[1], &vrts_[7]);
  edges_.emplace_back(&vrts_[1], &vrts_[8]);
  edges_.emplace_back(&vrts_[2], &vrts_[3]);
  edges_.emplace_back(&vrts_[2], &vrts_[5]);
  edges_.emplace_back(&vrts_[2], &vrts_[9]);
  edges_.emplace_back(&vrts_[3], &vrts_[4]);
  edges_.emplace_back(&vrts_[3], &vrts_[5]);
  edges_.emplace_back(&vrts_[3], &vrts_[8]);
  edges_.emplace_back(&vrts_[4], &vrts_[5]);
  edges_.emplace_back(&vrts_[4], &vrts_[8]);
  edges_.emplace_back(&vrts_[4], &vrts_[10]);
  edges_.emplace_back(&vrts_[4], &vrts_[11]);
  edges_.emplace_back(&vrts_[5], &vrts_[9]);
  edges_.emplace_back(&vrts_[5], &vrts_[11]);
  edges_.emplace_back(&vrts_[6], &vrts_[7]);
  edges_.emplace_back(&vrts_[6], &vrts_[9]);
  edges_.emplace_back(&vrts_[6], &vrts_[10]);
  edges_.emplace_back(&vrts_[6], &vrts_[11]);
  edges_.emplace_back(&vrts_[7], &vrts_[8]);
  edges_.emplace_back(&vrts_[7], &vrts_[10]);
  edges_.emplace_back(&vrts_[8], &vrts_[10]);
  edges_.emplace_back(&vrts_[9], &vrts_[11]);
  edges_.emplace_back(&vrts_[10], &vrts_[11]);
}

void TriMeshBoundary::ProjectToUnitSphere() {
  // p much just normalization i believe
  for (auto &&vrt : vrts_) {
    double norm{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      norm += SQR(vrt.pos_[i_dim]);
    }
    norm = sqrt(norm);
    double new_pos[3];
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      double val{vrt.pos_[i_dim]};
      new_pos[i_dim] = r_sys_ * val / norm;
    }
    vrt.SetPos(new_pos);
  }
}

void TriMeshBoundary::InitializeMesh() {
  TriMesh::InitializeMesh();
  double l_sum{0.0};
  for (auto const &edge : edges_) {
    l_sum += edge.length_;
  }
  double area_sum{0.0};
  for (auto const &tri : tris_) {
    area_sum += tri.area_;
  }
  // Using updated centroid position, calculate each triangle volume
  double vol_sum{0.0};
  for (auto &&tri : tris_) {
    tri.UpdateVolume(centroid_);
    vol_sum += tri.volume_;
  }
  // Report statistics
  size_t n_flawed{0};
  size_t n_gucci{0};
  for (auto &&vrt : vrts_) {
    if (vrt.n_neighbs_ != 6) {
      n_flawed++;
    } else if (vrt.n_neighbs_ == 6) {
      n_gucci++;
    } else {
      printf("ummm?? in TRIANGLE MESH INIT\n");
      exit(1);
    }
  }
  printf("Final mesh statistics:\n");
  printf("  %zu vrts (%zu flawed; %zu ideal)\n", vrts_.size(), n_flawed,
         n_gucci);
  l_avg_ = l_sum / edges_.size();
  A_prime_ = area_sum / tris_.size();
  V_prime_ = vol_sum / tris_.size();
  printf("  l_avg = %g\n", l_avg_);
  printf("  A_prime = %g\n", A_prime_);
  printf("  A_calc = %g\n", 4.0 * M_PI * SQR(r_sys_) / tris_.size());
  printf("  V_prime = %g\n", V_prime_);
  printf("  V_calc = %g\n", (4.0 / 3.0) * M_PI * CUBE(r_sys_) / tris_.size());
  printf("  V_calc_alt = %g\n", (1.0 / 3.0) * A_prime_ * r_sys_);
  l_c0_ = 1.2 * l_avg_;
  l_c1_ = 0.8 * l_avg_;
  // used in vutukuri et al
  // l_max_ = 1.4 * l_avg_;
  // l_min_ = 0.6 * l_avg_;
  // used when more edge flipping is desired
  l_max_ = 1.67 * l_avg_;
  l_min_ = 0.33 * l_avg_;
}

void TriMeshBoundary::FlipEdges() {
  // Shuffle edge indices them so that order of flipping is random
  size_t i_entries[edges_.size()];
  for (auto &&edge : edges_) {
    i_entries[edge.i_] = edge.i_;
    edge.just_flipped = false;
  }
  bool flipparino{false};
  rng_->Shuffle(i_entries, edges_.size());
  for (int i_edge{0}; i_edge < edges_.size(); i_edge++) {
    Edge *edge = &edges_[i_entries[i_edge]];
    if (edge->just_flipped) {
      continue;
    }
    Vertex *vrt1 = edge->vrts_[0];
    Vertex *vrt2 = edge->vrts_[1];
    // Find adjacent triangles
    Triangle *left{edge->tris_[0]};
    Triangle *right{edge->tris_[1]};
    if (left == nullptr or right == nullptr) {
      printf("issue at edge %zu\n", edge->i_);
      printf("ABORT\n");
      exit(1);
    }
    Vertex *vrt3{left->GetOtherVertex(vrt1, vrt2)};
    Vertex *vrt4{right->GetOtherVertex(vrt1, vrt2)};
    // Get pre- and post-flip edge lengths
    // double l_12{edge.length_};
    double l_34{0.0};
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      l_34 += SQR(vrt3->pos_[i_dim] - vrt4->pos_[i_dim]);
    }
    l_34 = sqrt(l_34);
    // If post-flip length is outside of allowed bounds, automatically discard it
    if (l_34 >= 0.9 * l_max_ or l_34 <= 1.1 * l_min_) {
      continue;
    }
    // check that angle between edges is >90 deg (edges cross post-flip otherwise)
    // manual test implementation
    double r_12[3]; // 2->1
    double r_32[3]; // 2->3
    double r_42[3]; // 2->4
    double r_21[3]; // 1->2
    double r_31[3]; // 1->4
    double r_41[3]; // 1->3
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      r_12[i_dim] = vrt1->pos_[i_dim] - vrt2->pos_[i_dim];
      r_32[i_dim] = vrt3->pos_[i_dim] - vrt2->pos_[i_dim];
      r_42[i_dim] = vrt4->pos_[i_dim] - vrt2->pos_[i_dim];
      r_21[i_dim] = vrt2->pos_[i_dim] - vrt1->pos_[i_dim];
      r_31[i_dim] = vrt3->pos_[i_dim] - vrt1->pos_[i_dim];
      r_41[i_dim] = vrt4->pos_[i_dim] - vrt1->pos_[i_dim];
    }
    double l_12{0.0};
    double l_32{0.0};
    double l_42{0.0};
    double l_21{0.0};
    double l_31{0.0};
    double l_41{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      l_12 += SQR(r_12[i_dim]);
      l_32 += SQR(r_32[i_dim]);
      l_42 += SQR(r_42[i_dim]);
      l_21 += SQR(r_21[i_dim]);
      l_31 += SQR(r_31[i_dim]);
      l_41 += SQR(r_41[i_dim]);
    }
    l_12 = sqrt(l_12);
    l_32 = sqrt(l_32);
    l_42 = sqrt(l_42);
    l_21 = sqrt(l_21);
    l_31 = sqrt(l_31);
    l_41 = sqrt(l_41);
    // we dont want any angle between to-be flipped edge to be greater than 90 deg
    // otherwise, you get eges that pass through each other
    double theta_31{acos(dot_product(3, r_21, r_31) / (l_21 * l_31))};
    double theta_32{acos(dot_product(3, r_12, r_32) / (l_12 * l_32))};
    double theta_41{acos(dot_product(3, r_21, r_41) / (l_21 * l_41))};
    double theta_42{acos(dot_product(3, r_12, r_42) / (l_12 * l_42))};
    double lim{0.75 * M_PI_2}; // dont want angles anywhere near 90 deg
    if (theta_31 >= lim or theta_32 >= lim or theta_41 >= lim or
        theta_42 >= lim) {
      continue;
    }

    Edge *edge_23{left->GetEdge(vrt2, vrt3)};
    Edge *edge_31{left->GetEdge(vrt3, vrt1)};
    Edge *edge_24{right->GetEdge(vrt2, vrt4)};
    Edge *edge_41{right->GetEdge(vrt4, vrt1)};
    double l_23{edge_23->length_};
    // double l_31{edge_31->length_};
    double l_24{edge_24->length_};
    // double l_41{edge_41->length_};
    // Get area of pre-flip triangles (1->2->3->1 and 1->2->4->1)
    double area_left{left->area_};
    double area_right{right->area_};
    // Calculate areas of post-flip triangles (1->3->4->1 and 2->3->4->2)
    double s1{0.5 * (l_34 + l_41 + l_31)};
    double s2{0.5 * (l_34 + l_24 + l_23)};
    double area_top{sqrt(s1 * (s1 - l_34) * (s1 - l_41) * (s1 - l_31))};
    double area_bot{sqrt(s2 * (s2 - l_34) * (s2 - l_24) * (s2 - l_23))};
    // If post-flip triangles are wonnky (NaN area), automatically discard it
    if (area_top != area_top or area_bot != area_bot) {
      continue;
    }
    // Calculate energy of current configuration
    double current_energy{0.0};
    if (edge->length_ > l_c0_) {
      current_energy += kappa_B_ * exp(1.0 / (l_c0_ - l_12)) / (l_max_ - l_12);
    } else if (edge->length_ < l_c1_) {
      current_energy += kappa_B_ * exp(1.0 / (l_12 - l_c1_)) / (l_12 - l_min_);
    }
    current_energy += 0.5 * kappa_l_ * SQR(A_prime_ - area_left) / A_prime_;
    current_energy += 0.5 * kappa_l_ * SQR(A_prime_ - area_right) / A_prime_;
    // Calculate energy of configuration post-flip
    double postflip_energy{0.0};
    if (l_34 > l_c0_) {
      postflip_energy += kappa_B_ * exp(1.0 / (l_c0_ - l_34)) / (l_max_ - l_34);
    } else if (l_34 < l_c1_) {
      postflip_energy += kappa_B_ * exp(1.0 / (l_34 - l_c1_)) / (l_34 - l_min_);
    }
    postflip_energy += 0.5 * kappa_l_ * SQR(A_prime_ - area_top) / A_prime_;
    postflip_energy += 0.5 * kappa_l_ * SQR(A_prime_ - area_bot) / A_prime_;
    bool do_flip{false};
    // bool do_flip{true};
    double p_flip_base{0.3}; // from vutukuri -- seems weird but OK why not
    if (rng_->RandomUniform() < p_flip_base) {
      if (postflip_energy <= current_energy) {
        do_flip = true;
      } else {
        double delta{postflip_energy - current_energy};
        double p_flip{exp(-delta)};
        double ran = rng_->RandomUniform();
        if (ran < p_flip) {
          do_flip = true;
        }
      }
    }
    // flipped triangles: 3->2->4 and 3->1->4
    if (do_flip) {
      flipparino = true;
      // printf("FLIP edge %zu of triangles %zu and %zu\n", edge->i_, left->i_,
      //        right->i_);
      // printf("angles: %g and %g\n", angle_left * 180 / M_PI,
      //        angle_right * 180 / M_PI);
      edge->vrts_[0] = vrt3;
      edge->vrts_[1] = vrt4;
      // turn left triangle into top triangle
      left->vrts_[0] = vrt1;
      left->vrts_[1] = vrt3;
      left->vrts_[2] = vrt4;
      left->edges_[0] = edge;
      left->edges_[1] = edge_31;
      left->edges_[2] = edge_41;
      // turn right triangle into bottom triangle
      right->vrts_[0] = vrt2;
      right->vrts_[1] = vrt3;
      right->vrts_[2] = vrt4;
      right->edges_[0] = edge;
      right->edges_[1] = edge_23;
      right->edges_[2] = edge_24;
      // need to mark all edges as just flipped; otherwise can get wonky triangles
      for (int i_edge{0}; i_edge < 3; i_edge++) {
        left->edges_[i_edge]->just_flipped = true;
        right->edges_[i_edge]->just_flipped = true;
      }
      // dynamically updating causes minor error in vol (~0.1%), but improves performance
      // (error is from the fact that edge flips change the centroid slightly)
      edge->Update();
      left->UpdateArea();
      right->UpdateArea();
      left->UpdateVolume(centroid_);
      right->UpdateVolume(centroid_);
      // vrt1->SetColor(0.0, draw_type::fixed);
      // vrt2->SetColor(0.0, draw_type::fixed);
      // vrt3->SetColor(0.0, draw_type::fixed);
      // vrt4->SetColor(0.0, draw_type::fixed);
    }
  }
  if (flipparino) {
    UpdateNeighbors();
    // uncomment below if you want to remove dynamic updating above
    // UpdateMesh();
  }
}

void TriMeshBoundary::ApplyMembraneForces() {

  // SF TODO all calculations are currently done redundantly
  // SF TODO i.e., we do not take advtange of newton's 3rd law
  // SF TODO once validated, increase computational efficiency by addressing this

  double f_teth_sum{0.0};
  size_t n_teth{0};
  double f_bend_sum{0.0};
  size_t n_bend{0};
  double f_area_sum{0.0};
  size_t n_area{0};
  double f_vol_sum{0.0};
  size_t n_vol{0};

  // SAF TODO -- make loop over edges (with pre-updated r vec+mag)
  // TETHER FORCE
  for (auto &&vrt : vrts_) {
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
      Vertex *neighb{vrt.neighbs_[i_neighb]};
      if (neighb == nullptr) {
        printf("what in TriMesh::UpdatePosition()\n");
        exit(1);
      }
      double r_ij[3];
      double rmag{0.0};
      for (int i{0}; i < 3; i++) {
        r_ij[i] = vrt.GetPosition()[i] - neighb->GetPosition()[i];
        rmag += SQR(r_ij[i]);
      }
      rmag = sqrt(rmag);
      // SF TODO need some way of getting back into range w/o infinite forces
      // use a RLY STRONG SPRING ???
      if (rmag >= l_max_ or rmag <= l_min_) {
        printf("CORRECTION!\n");
        double kspring{10};
        double mag{kspring * (l_avg_ - rmag)}; // use l_avg as rest length;
        double vec[3];
        for (int i{0}; i < 3; i++) {
          vec[i] = mag * r_ij[i] / rmag;
        }
        vrt.AddForce(vec);
        continue;
      }
      // free movement within a certain range
      if (rmag <= l_c0_ and rmag >= l_c1_) {
        continue;
      }
      double fmag{0.0};
      // attraction
      if (rmag > l_c0_) {
        fmag = (1 / (l_max_ - rmag)) *
               (SQR(1.0 / (rmag - l_c1_)) - (1.0 / (l_max_ - rmag))) *
               exp(1.0 / (rmag - l_c1_));
      }
      // repulsion
      if (rmag < l_c1_) {
        fmag = (1.0 / (rmag - l_min_)) *
               ((1.0 / (rmag - l_min_)) - SQR(1.0 / (l_c0_ - rmag))) *
               exp(1.0 / (l_c0_ - rmag));
      }
      fmag *= kappa_B_;
      double f[3];
      for (int i{0}; i < 3; i++) {
        f[i] = fmag * r_ij[i] / rmag;
      }
      vrt.AddForce(f);
      f_teth_sum += fmag;
      n_teth++;
    }
  }
  // DISCRETE BENDING FORCES
  for (auto &&vrt : vrts_) {
    double sum_lsqT{0.0};        // scalar
    double sum_del_lsqT[3]{{}};  // vec; scalar in each dim
    double sum_rT[3]{{}};        // vec; scalr in each dim
    double sum_del_rT[3][3]{{}}; // tensor; vec in each dim
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
      // this assumes neighbors are ordered in a ring with + oriented ccw
      int i_plus{i_neighb == (vrt.n_neighbs_ - 1) ? 0 : i_neighb + 1};
      int i_minus{i_neighb == 0 ? vrt.n_neighbs_ - 1 : i_neighb - 1};
      Vertex *neighb{vrt.neighbs_[i_neighb]};
      Vertex *neighb_plus{vrt.neighbs_[i_plus]};
      Vertex *neighb_minus{vrt.neighbs_[i_minus]};
      double r_ij[3]; // points from j (neighb) to i (vrt)
      double r_ij_plus[3];
      double r_ij_minus[3];
      double r_jj_plus[3];
      double r_jj_minus[3];
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        r_ij[i_dim] = vrt.pos_[i_dim] - neighb->pos_[i_dim];
        r_ij_plus[i_dim] = vrt.pos_[i_dim] - neighb_plus->pos_[i_dim];
        r_ij_minus[i_dim] = vrt.pos_[i_dim] - neighb_minus->pos_[i_dim];
        r_jj_plus[i_dim] = neighb->pos_[i_dim] - neighb_plus->pos_[i_dim];
        r_jj_minus[i_dim] = neighb->pos_[i_dim] - neighb_minus->pos_[i_dim];
      }
      // the edge l_ij = |r_ij| connects two triangles, get length of each
      double l_ij{0.0};
      double l_ij_plus{0.0};
      double l_ij_minus{0.0};
      double l_jj_plus{0.0};
      double l_jj_minus{0.0};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        l_ij += SQR(r_ij[i_dim]);
        l_ij_plus += SQR(r_ij_plus[i_dim]);
        l_ij_minus += SQR(r_ij_minus[i_dim]);
        l_jj_plus += SQR(r_jj_plus[i_dim]);
        l_jj_minus += SQR(r_jj_minus[i_dim]);
      }
      l_ij = sqrt(l_ij);
      l_ij_plus = sqrt(l_ij_plus);
      l_ij_minus = sqrt(l_ij_minus);
      l_jj_plus = sqrt(l_jj_plus);
      l_jj_minus = sqrt(l_jj_minus);
      double chi_minus{dot_product(3, r_ij_minus, r_jj_minus) /
                       (l_ij_minus * l_jj_minus)};
      double chi_plus{dot_product(3, r_ij_plus, r_jj_plus) /
                      (l_ij_plus * l_jj_plus)};
      double T_ij{(chi_minus / sqrt(1.0 - SQR(chi_minus))) +
                  (chi_plus / sqrt(1.0 - SQR(chi_plus)))};
      double grad_lsq[3];
      double grad_chi_plus[3];
      double grad_chi_minus[3];
      double grad_T[3];
      // c.f. Appendix A of Guckenberger et al. Comp. Phys. Comm (2016)
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        grad_lsq[i_dim] = 2 * r_ij[i_dim];
        grad_chi_plus[i_dim] =
            (1.0 / (l_ij_plus * l_jj_plus)) *
            (r_jj_plus[i_dim] -
             (l_jj_plus / l_ij_plus) * chi_plus * r_ij_plus[i_dim]);
        grad_chi_minus[i_dim] =
            (1.0 / (l_ij_minus * l_jj_minus)) *
            (r_jj_minus[i_dim] -
             (l_jj_minus / l_ij_minus) * chi_minus * r_ij_minus[i_dim]);
        grad_T[i_dim] =
            grad_chi_plus[i_dim] / std::pow((1.0 - SQR(chi_plus)), 1.5) +
            grad_chi_minus[i_dim] / std::pow((1.0 - SQR(chi_minus)), 1.5);
      }
      sum_lsqT += SQR(l_ij) * T_ij;
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        sum_rT[i_dim] += r_ij[i_dim] * T_ij;
        sum_del_lsqT[i_dim] +=
            (T_ij * grad_lsq[i_dim] + SQR(l_ij) * grad_T[i_dim]);
        for (int j_dim{0}; j_dim < 3; j_dim++) {
          sum_del_rT[i_dim][j_dim] += r_ij[j_dim] * grad_T[i_dim];
          if (i_dim == j_dim) {
            sum_del_rT[i_dim][j_dim] += T_ij;
          }
        }
      }
    }
    // then we can use these weights to find vector forces
    double f_bend[3]{{}};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      f_bend[i_dim] += 2 * dot_product(3, sum_rT, sum_rT) / SQR(sum_lsqT) *
                       sum_del_lsqT[i_dim];
      f_bend[i_dim] -= 4 * dot_product(3, sum_rT, sum_del_rT[i_dim]) / sum_lsqT;
      f_bend[i_dim] *= kappa_;
    }
    vrt.AddForce(f_bend);
    double fmag_bend{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      fmag_bend += SQR(f_bend[i_dim]);
    }
    f_bend_sum += sqrt(fmag_bend);
    n_bend++;
  }
  // SF TODO loop over pre-updated triangles?
  // AREA AND VOLUME CONSERVATION FORCE
  for (auto &&vrt : vrts_) {
    double f_area[3]{{}};
    double f_vol[3]{{}};
    for (int i_neighb{0}; i_neighb < vrt.n_neighbs_; i_neighb++) {
      // this assumes neighbors are ordered in a ring with + oriented ccw
      int i_plus{i_neighb == (vrt.n_neighbs_ - 1) ? 0 : i_neighb + 1};
      Vertex *neighb{vrt.neighbs_[i_neighb]};
      Vertex *neighb_plus{vrt.neighbs_[i_plus]};
      // only consider contribution from fwd triangle to avoid double counting
      Triangle *tri_plus{vrt.tris_[i_neighb]};
      double r_jj_plus[3]; // points from j_plus (fwd_neighb) to j (neighb)
      double r_oj[3];      // points from j (neighb) to o (centroid)
      double r_oj_plus[3]; // points from j_plus (fwd_neighb) to o (centroid)
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        r_jj_plus[i_dim] = neighb->pos_[i_dim] - neighb_plus->pos_[i_dim];
        r_oj_plus[i_dim] = centroid_[i_dim] - neighb_plus->GetPosition()[i_dim];
        r_oj[i_dim] = centroid_[i_dim] - neighb->GetPosition()[i_dim];
      }
      // the edge l_ij = |r_ij| connects two triangles, get length of each
      double l_jj_plus{0.0};
      double l_oj_plus{0.0};
      double l_oj{0.0};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        l_jj_plus += SQR(r_jj_plus[i_dim]);
        l_oj_plus += SQR(r_oj_plus[i_dim]);
        l_oj += SQR(r_oj[i_dim]);
      }
      l_jj_plus = sqrt(l_jj_plus);
      l_oj_plus = sqrt(l_oj_plus);
      l_oj = sqrt(l_oj);
      // calculate force from area conservation
      // cross product order intentionally reversed cuz dir of r_jj+ is wrong
      double f_area_vec[3];
      cross_product(r_jj_plus, tri_plus->nhat_, f_area_vec, 3);
      double area{tri_plus->area_};
      double f_area_mag{-0.25 * kappa_l_ * (area - A_prime_) / A_prime_};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        f_area[i_dim] += f_area_mag * f_area_vec[i_dim];
      }
      // need area of triangle that connects neighb and neighb_plus w/ polygon center
      double s_opp{0.5 * (l_jj_plus + l_oj + l_oj_plus)};
      double area_opp{sqrt(s_opp * (s_opp - l_jj_plus) * (s_opp - l_oj) *
                           (s_opp - l_oj_plus))};
      double V{tri_plus->volume_};
      // double epsilon{1e-1};
      // if (std::fabs(area - tri_plus->area_) > epsilon) {
      //   printf("bruh\n");
      //   printf("edge lengths: %g | %g | %g\n", l_ij, l_ij_plus, l_jj_plus);
      //   printf("area: %g vs %g\n", area, tri_plus->area_);
      // }
      // if (std::fabs(V - tri_plus->volume_) > epsilon) {
      //   printf("bruh\n");
      //   printf("edge lengths: %g | %g | %g\n", l_ij, l_ij_plus, l_jj_plus);
      //   printf("volume: %g vs %g\n", V, tri_plus->volume_);
      //   for (int i{0}; i < 3; i++) {
      //     printf("  edge %zu\n", tri_plus->edges_[i]->i_);
      //   }
      // }
      double f_vol_mag{(-1.0 / 3.0) * kappa_v_ * (V - V_prime_) * area_opp /
                       V_prime_};
      if (f_vol_mag != f_vol_mag) {
        vrt.SetColor(0.0, draw_type::fixed);
        neighb_plus->SetColor(1.0, draw_type::fixed);
        printf("V: %g\n", V);
        printf("V_prime: %g\n", V_prime_);
        printf("area_opp: %g\n", area_opp);
        do_not_pass_go_ = true;
        return;
      }
      // just the inward pointing normal of the face opposite to vrt in prism
      double f_vol_vec[3];
      cross_product(r_oj_plus, r_jj_plus, f_vol_vec, 3);
      double f_vol_vec_norm{0.0};
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        f_vol_vec_norm += SQR(f_vol_vec[i_dim]);
      }
      f_vol_vec_norm = sqrt(f_vol_vec_norm);
      if (f_vol_vec_norm == 0.0) {
        printf("r_oj+ = <%g, %g, %g>\n", r_oj_plus[0], r_oj_plus[1],
               r_oj_plus[2]);
        printf("r_jj+ = <%g, %g, %g>\n", r_jj_plus[0], r_jj_plus[1],
               r_jj_plus[2]);
        printf("NOPE\n");
        exit(1);
      }
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        f_vol[i_dim] += f_vol_mag * f_vol_vec[i_dim] / f_vol_vec_norm;
      }
    }
    vrt.AddForce(f_area);
    vrt.AddForce(f_vol);
    double fmag_area{0.0};
    double fmag_vol{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      fmag_area += SQR(f_area[i_dim]);
      fmag_vol += SQR(f_vol[i_dim]);
    }
    f_area_sum += sqrt(fmag_area);
    n_area++;
    f_vol_sum += sqrt(fmag_vol);
    n_vol++;
  }
  f_avgs_[0] = f_teth_sum / n_teth;
  f_avgs_[1] = f_bend_sum / n_bend;
  f_avgs_[2] = f_area_sum / n_area;
  f_avgs_[3] = f_vol_sum / n_vol;
}

void TriMeshBoundary::ApplyBoundaryForces() {

  if (boundary_neighbs_.empty()) {
    return;
  }
  // values taken directly from NAB
  double r_cutoff2 = SQR(pow(2, 1.0 / 6.0) * 0.5);
  double sigma2 = 0.25;
  double four_epsilon = 4.0;
  int n_bonds_tot{0};
  for (int i_neighb{0}; i_neighb < boundary_neighbs_.size(); i_neighb++) {
    Object *neighb{boundary_neighbs_[i_neighb]};
    Filament *fil{dynamic_cast<Filament *>(neighb)};
    n_bonds_tot += fil->GetNBonds();
  }
  int i_bond{0};
  f_mem_.resize(n_bonds_tot);
  for (int i_neighb{0}; i_neighb < boundary_neighbs_.size(); i_neighb++) {
    Object *neighb{boundary_neighbs_[i_neighb]};
    Filament *fil{dynamic_cast<Filament *>(neighb)};
    int n_bonds = fil->GetNBonds();
    for (int j_bond{0}; j_bond < n_bonds; j_bond++) {
      // printf("bond %i\n", j_bond);
      double rmin[3], r_min_mag2, rcontact[3], mu;
      Bond *bond{fil->GetBond(j_bond)};
      double r[3] = {bond->GetPosition()[0], bond->GetPosition()[1],
                     bond->GetPosition()[2]};
      double u[3] = {bond->GetOrientation()[0], bond->GetOrientation()[1],
                     bond->GetOrientation()[2]};
      double l{bond->GetLength()};
      // SF TODO incorporate periodic space (incorporate s and h)
      int j_tri = mindist_.SpheroPolygon(this, r, r, u, l, rmin, &r_min_mag2,
                                         rcontact, &mu);
      f_mem_[i_bond].color = 1.8 * M_PI;
      f_mem_[i_bond].diameter = 0.25;
      f_mem_[i_bond].draw = draw_type::fixed;
      double u_mindist[3];
      double l_mindist{sqrt(r_min_mag2)};
      // when distance is zero, apply cutoff force -- SF TODO validate
      if (l_mindist == 0) {
        f_mem_[i_bond].color = 0.0 * M_PI;
        for (int i_dim{0}; i_dim < 3; i_dim++) {
          f_mem_[i_bond].r[i_dim] = rcontact[i_dim];
        }
        f_mem_[i_bond].diameter = 3;
        f_mem_[i_bond].length = 0;
        double f_cutoff{0.1 / params_->delta * gamma_};
        double f_lj[3] = {0.0};
        for (int i = 0; i < params_->n_dim; ++i) {
          f_lj[i] = f_cutoff;
        }
        Triangle *tri{&tris_[j_tri]};
        for (auto &&vrt : tri->vrts_) {
          vrt->AddForce(f_lj);
        }
        fil->SubForce(f_lj);
        continue;
        // do_not_pass_go_ = true;
        // return;
      }
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        u_mindist[i_dim] = rmin[i_dim] / l_mindist;
      }
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        f_mem_[i_bond].r[i_dim] =
            rcontact[i_dim] - u_mindist[i_dim] * l_mindist / 2.0;
        f_mem_[i_bond].u[i_dim] = u_mindist[i_dim];
      }
      f_mem_[i_bond].length = l_mindist;
      if (r_min_mag2 < r_cutoff2) {
        f_mem_[i_bond].color = 0.0 * M_PI;
        f_mem_[i_bond].diameter = 0.5;
        // Calculate WCA potential and forces
        double rho2 = sigma2 / r_min_mag2;
        double rho6 = CUBE(rho2);
        double rho12 = SQR(rho6);

        // u += four_epsilon * (rho12 - rho6) + u_shift;
        double factor = 6.0 * four_epsilon * (2.0 * rho12 - rho6) / r_min_mag2;
        // printf("r_min_mag = %g\n", sqrt(r_min_mag2));
        // printf("factor = %g\n", factor);
        double f_cutoff{0.1 / params_->delta * gamma_};
        // Truncate the forces if necessary
        double r_min_mag = sqrt(r_min_mag2);
        if (factor * r_min_mag > f_cutoff) {
          // printf(" *** Force exceeded f_cutoff "
          //        "kinetochoremesh_mt_wca_potential_neighbors ***\n");
          // printf("r_min_mag: %g\n", r_min_mag);
          // printf("factor: %g\n", factor);
          factor = f_cutoff / r_min_mag;
        }
        double f_lj[3] = {0.0};
        for (int i = 0; i < params_->n_dim; ++i) {
          f_lj[i] = factor * rmin[i];
        }
        Triangle *tri{&tris_[j_tri]};
        for (auto &&vrt : tri->vrts_) {
          vrt->AddForce(f_lj);
        }
        fil->SubForce(f_lj);
        // SF TODO incorporate torques
        /*
      // Calculate torques
      double rcontact_kc[3] = {0.0};
      double rcontact_mt[3] = {0.0};
      for (int i = 0; i < ndim; ++i) {
        //rcontact_kc[i] = kc_iter->r_[i] - rcontact[i];
        rcontact_kc[i] = chromosomes->r_[ikc][i] - rcontact[i];
        rcontact_mt[i] = mu * u_bond[ibond][i];
      }
      //std::cout << "rcontact_kc(" << rcontact_kc[0] << ", " << rcontact_kc[1] << ", " << rcontact_kc[2] << ")\n";
      //std::cout << "rcontact_mt(" << rcontact_mt[0] << ", " << rcontact_mt[1] << ", " << rcontact_mt[2] << ")\n";
      double tau[3] = {0.0};
      cross_product(rcontact_kc, f_lj, tau, 3);
      for (int i = 0; i < 3; ++i) {
        t_kc[ikc][i] -= tau[i];
      }
      cross_product(rcontact_mt, f_lj, tau, 3);
      for (int i = 0; i < 3; ++i) {
        t_bond[ibond][i] -= tau[i];
      }
      */
      }
      i_bond++;
    }
  }
}

void TriMeshBoundary::WriteOutputs() {
  // SF TODO temp before integrating into output_manager
  if (params_->i_step % params_->mesh_steps_per_datapoint != 0) {
    return;
  }
  if (i_datapoint_ < params_->mesh_datapoints) {
    // write different membrane forces
    int n_written = fwrite(f_avgs_, sizeof(double), 4, forces_);
    i_datapoint_++;
    for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
      double pos[3];
      int adj[n_edges_max_];
      Vertex *vrt{&vrts_[i_vrt]};
      for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
        pos[i_dim] = vrt->GetPosition()[i_dim];
      }
      for (int i_neighb{0}; i_neighb < vrt->n_neighbs_; i_neighb++) {
        adj[i_neighb] = (int)vrt->neighbs_[i_neighb]->i_;
      }
      // padding for non-ideal vertices (less or more than 6 neighbs)
      for (int i_neighb{vrt->n_neighbs_}; i_neighb < n_edges_max_; i_neighb++) {
        adj[i_neighb] = -1;
      }
      fwrite(pos, sizeof(double), 3, vertices_);
      fwrite(adj, sizeof(int), n_edges_max_, adjacency_);
    }
  } else {
    printf("data collection done! (%i datapoints written)\n", i_datapoint_);
    fclose(forces_);
    fclose(vertices_);
    fclose(adjacency_);
    early_exit = true;
  }
}

void TriMeshBoundary::Draw(std::vector<graph_struct *> &graph_array) {
  TriMesh::Draw(graph_array);
  if (params_->draw_mindist) {
    for (int i{0}; i < f_mem_.size(); i++) {
      graph_array.push_back(&f_mem_[i]);
    }
  }
}

void TriMeshBoundary::UpdatePositions() {

  if (do_not_pass_go_) {
    return;
  }
  // shrink this jonny
  if (params_->mesh_shrink_rate > 0.0) {
    l_avg_ *= (1.0 - params_->mesh_shrink_rate);
    A_prime_ *= SQR(1.0 - params_->mesh_shrink_rate);
    V_prime_ *= CUBE(1.0 - params_->mesh_shrink_rate);
    l_c0_ = 1.2 * l_avg_;
    l_c1_ = 0.8 * l_avg_;
    l_max_ = 1.667 * l_avg_;
    l_min_ = 0.333 * l_avg_;
  }
  UpdateMesh(); // update edge lengths, triangle area/vol, etc.
  if (do_not_pass_go_) {
    return;
  }
  ApplyMembraneForces();
  if (do_not_pass_go_) {
    return;
  }
  ApplyBoundaryForces();
  if (do_not_pass_go_) {
    return;
  }
  if (params_->enable_flipping) {
    FlipEdges();
    if (do_not_pass_go_) {
      return;
    }
  }
  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    double sigma{sqrt(2 * params_->delta / gamma_)}; // kbT = 1
    double dr[3] = {0.0, 0.0, 0.0};
    double const *const r_prev = vrts_[i_vrt].GetPosition();
    double r_final[3];
    for (int i = 0; i < 3; ++i) {
      double vel{vrts_[i_vrt].GetForce()[i] / gamma_};
      double noise{rng_->RandomNormal(sigma)};
      dr[i] = vel * params_->delta + noise;
      r_final[i] = r_prev[i] + dr[i];
    }
    if (i_vrt == 0) {
      r_final[0] += 0.01;
    }
    vrts_[i_vrt].SetPos(r_final);
  }
  WriteOutputs();
}