#include <cglass/filament.hpp>
#include <cglass/triangle_mesh.hpp>
#include <unistd.h>

void TriMesh::Init(system_parameters *params) {
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

void TriMesh::SetParameters() {
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

void TriMesh::MakeIcosphere() {
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

void TriMesh::MakeIcosahedron() {
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

void TriMesh::DivideFaces() {
  std::vector<Triangle> new_faces;
  std::vector<Edge> new_edges;
  size_t n_vrts_pre{vrts_.size()};
  // loop over all triangles in mesh
  for (auto &&face : tris_) {
    Vertex *vrt0{face.vrts_[0]};
    Vertex *vrt1{face.vrts_[1]};
    Vertex *vrt2{face.vrts_[2]};
    // each edge is bisected to generate a finer mesh
    double pos10[3] = {(vrt0->pos_[0] + vrt1->pos_[0]) / 2.0,
                       (vrt0->pos_[1] + vrt1->pos_[1]) / 2.0,
                       (vrt0->pos_[2] + vrt1->pos_[2]) / 2.0};
    double pos12[3] = {(vrt2->pos_[0] + vrt1->pos_[0]) / 2.0,
                       (vrt2->pos_[1] + vrt1->pos_[1]) / 2.0,
                       (vrt2->pos_[2] + vrt1->pos_[2]) / 2.0};
    double pos20[3] = {(vrt0->pos_[0] + vrt2->pos_[0]) / 2.0,
                       (vrt0->pos_[1] + vrt2->pos_[1]) / 2.0,
                       (vrt0->pos_[2] + vrt2->pos_[2]) / 2.0};
    // ensure we do not add new vertices more than once
    Vertex *vrt10{nullptr}, *vrt12{nullptr}, *vrt20{nullptr};
    bool dup10{false}, dup12{false}, dup20{false};
    bool error{false};
    for (auto &&vrt : vrts_) {
      if (vrt.pos_[0] == pos10[0] and vrt.pos_[1] == pos10[1] and
          vrt.pos_[2] == pos10[2]) {
        if (vrt10 != nullptr) {
          error = true;
        }
        dup10 = true;
        vrt10 = &vrt;
      }
      if (vrt.pos_[0] == pos12[0] and vrt.pos_[1] == pos12[1] and
          vrt.pos_[2] == pos12[2]) {
        if (vrt12 != nullptr) {
          error = true;
        }
        dup12 = true;
        vrt12 = &vrt;
      }
      if (vrt.pos_[0] == pos20[0] and vrt.pos_[1] == pos20[1] and
          vrt.pos_[2] == pos20[2]) {
        if (vrt20 != nullptr) {
          error = true;
        }
        dup20 = true;
        vrt20 = &vrt;
      }
    }
    if (error) {
      printf("Critical error in TriMesh::DivideFaces()\n");
      exit(1);
    }
    if (!dup10) {
      vrts_.emplace_back(pos10[0], pos10[1], pos10[2]);
      vrt10 = &vrts_.back();
    }
    if (!dup12) {
      vrts_.emplace_back(pos12[0], pos12[1], pos12[2]);
      vrt12 = &vrts_.back();
    }
    if (!dup20) {
      vrts_.emplace_back(pos20[0], pos20[1], pos20[2]);
      vrt20 = &vrts_.back();
    }
    if (vrt10 == nullptr or vrt12 == nullptr or vrt20 == nullptr) {
      printf("nope\n");
      exit(1);
    }
    // 1 face turns into 4, and 3 edges turn into 9
    new_faces.emplace_back(vrt0, vrt10, vrt20);
    new_faces.emplace_back(vrt1, vrt10, vrt12);
    new_faces.emplace_back(vrt2, vrt12, vrt20);
    new_faces.emplace_back(vrt10, vrt12, vrt20);
    // only add outer edges if vertex is new -- otherwise will double count
    if (!dup10) {
      new_edges.emplace_back(vrt0, vrt10);
      new_edges.emplace_back(vrt1, vrt10);
    }
    if (!dup20) {
      new_edges.emplace_back(vrt0, vrt20);
      new_edges.emplace_back(vrt2, vrt20);
    }
    if (!dup12) {
      new_edges.emplace_back(vrt1, vrt12);
      new_edges.emplace_back(vrt2, vrt12);
    }
    // inner edges are always added since they cannot be double counted
    new_edges.emplace_back(vrt10, vrt20);
    new_edges.emplace_back(vrt10, vrt12);
    new_edges.emplace_back(vrt12, vrt20);
  }
  printf("  %zu -> %zu triangles\n", tris_.size(), new_faces.size());
  printf("  %zu -> %zu vertices\n", n_vrts_pre, vrts_.size());
  printf("  %zu -> %zu edges\n", edges_.size(), new_edges.size());
  tris_ = new_faces;
  edges_ = new_edges;
}

void TriMesh::InitializeMesh() {
  // Initialize vertex storage + auxiliary parameters
  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    vrts_[i_vrt].i_ = i_vrt;
    vrts_[i_vrt].tris_.resize(n_edges_max_);
    vrts_[i_vrt].edges_.resize(n_edges_max_);
    vrts_[i_vrt].neighbs_.resize(n_edges_max_);
    vrts_[i_vrt].SetDiameter(params_->node_diameter);
    vrts_[i_vrt].SetColor(1.5, draw_type::fixed);
  }
  // Initialize edge indices and update lengths
  double l_sum{0.0};
  for (int i_edge{0}; i_edge < edges_.size(); i_edge++) {
    edges_[i_edge].i_ = i_edge;
    edges_[i_edge].Update();
    l_sum += edges_[i_edge].length_;
  }
  // Initialize triangle indices and manually calculate areas (edges not assigned yet)
  double area_sum{0.0};
  for (int i_tri{0}; i_tri < tris_.size(); i_tri++) {
    tris_[i_tri].i_ = i_tri;
    Triangle *tri{&tris_[i_tri]};
    double l1{0.0};
    double l2{0.0};
    double l3{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      l1 += SQR(tri->vrts_[0]->pos_[i_dim] - tri->vrts_[1]->pos_[i_dim]);
      l2 += SQR(tri->vrts_[0]->pos_[i_dim] - tri->vrts_[2]->pos_[i_dim]);
      l3 += SQR(tri->vrts_[1]->pos_[i_dim] - tri->vrts_[2]->pos_[i_dim]);
    }
    l1 = sqrt(l1);
    l2 = sqrt(l2);
    l3 = sqrt(l3);
    double s{0.5 * (l1 + l2 + l3)};
    tri->area_ = sqrt(s * (s - l1) * (s - l2) * (s - l3));
    area_sum += tri->area_;
  }
  // Now that triangle areas have been updated, calculate centroid position
  UpdateCentroid();
  // Using updated centroid position, calculate each triangle volume
  double vol_sum{0.0};
  for (auto &&tri : tris_) {
    tri.UpdateVolume(centroid_);
    vol_sum += tri.volume_;
  }
  // Update neighbor lists and pointers for vertices, edges, and triangles
  UpdateNeighbors();
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

void TriMesh::UpdateCentroid() {
  for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
    centroid_[i_dim] = 0.0;
  }
  double area_sum{0.0};
  for (auto &&tri : tris_) {
    for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
      centroid_[i_dim] += tri.area_ * tri.GetCenterPos(i_dim);
    }
    area_sum += tri.area_;
  }
  for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
    centroid_[i_dim] /= area_sum;
  }
  o_.r[0] = centroid_[0];
  o_.r[1] = centroid_[1];
  o_.r[2] = centroid_[2];
}

void TriMesh::UpdateNeighbors() {
  UpdateTriangles();
  // NOTE: currently only works with convex objects
  // NOTE: an update to ensure neighboring normals are aligned would fix this
  // NOTE: alternatively, a ray-tracing algorithm to determine in/out would be best
  // This function is purposefully coded to be explicit yet inefficient
  // (it does not get called often if at all beyond initialization)

  // Reset storage of neighb and triangle ptrs in each vertex
  for (auto &&vrt : vrts_) {
    vrt.n_tris_ = 0;
    vrt.n_neighbs_ = 0;
    vrt.n_edges_ = 0;
  }
  // Update triangle ptrs stored by each vertex and reset triangle ptrs
  for (auto &&tri : tris_) {
    tri.vrts_[0]->tris_[tri.vrts_[0]->n_tris_++] = &tri;
    tri.vrts_[1]->tris_[tri.vrts_[1]->n_tris_++] = &tri;
    tri.vrts_[2]->tris_[tri.vrts_[2]->n_tris_++] = &tri;
    for (auto &&edge : tri.edges_) {
      edge = nullptr;
    }
  }
  // Update edge + neighb ptrs stored by each vertex
  for (auto &&edge : edges_) {
    edge.vrts_[0]->edges_[edge.vrts_[0]->n_edges_++] = &edge;
    edge.vrts_[1]->edges_[edge.vrts_[1]->n_edges_++] = &edge;
    edge.vrts_[0]->neighbs_[edge.vrts_[0]->n_neighbs_++] = edge.vrts_[1];
    edge.vrts_[1]->neighbs_[edge.vrts_[1]->n_neighbs_++] = edge.vrts_[0];
  }
  // Update triangles and edges so that they track each other
  for (auto &&edge : edges_) {
    int n_found{0};
    // the two triangles we are a part of will be contained by each vertex
    for (int i_tri{0}; i_tri < edge.vrts_[0]->n_tris_; i_tri++) {
      Triangle *tri{edge.vrts_[0]->tris_[i_tri]};
      if (tri->Contains(edge.vrts_[1])) {
        // store triangle pointer in edge
        edge.tris_[n_found++] = tri;
        // add edge pointer to said triangle
        for (int i_edge{0}; i_edge < 3; i_edge++) {
          // printf("query %i\n", i_edge);
          if (tri->edges_[i_edge] == nullptr) {
            // printf("edge %i added\n", i_edge);
            tri->edges_[i_edge] = &edge;
            break;
          } else if (i_edge == 2) {
            printf("error @ edge %zu\n", edge.i_);
            exit(1);
          }
        }
      }
    }
    if (n_found != 2) {
      edge.vrts_[0]->SetColor(2 * M_PI, draw_type::fixed);
      edge.vrts_[0]->SetDiameter(2);
      printf("Error: found %i triangles that contain edge #%zu\n", n_found,
             edge.i_);
      do_not_pass_go_ = true;
      return;
      // exit(1);
    }
  }
  // Order triangle ptrs so that they are in consecutive order (share edges)
  for (auto &&vrt : vrts_) {
    Triangle *tris_ordered[vrt.n_tris_]{{}};
    Edge *edges_ordered[vrt.n_edges_]{{}};
    Vertex *neighbs_ordered[vrt.n_neighbs_]{{}};
    if (vrt.n_tris_ != vrt.n_edges_ or vrt.n_tris_ != vrt.n_neighbs_) {
      printf("Error in neighbor lists\n");
      exit(1);
    }
    // starting point doesn't matter, so just use whichever entry is first
    neighbs_ordered[0] = vrt.neighbs_[0];
    for (int i_edge{0}; i_edge < vrt.n_edges_; i_edge++) {
      Edge *edge{vrt.edges_[i_edge]};
      if (edge->Contains(&vrt) and edge->Contains(vrt.neighbs_[0])) {
        edges_ordered[0] = edge;
        break;
      }
    }
    if (edges_ordered[0] == nullptr) {
      printf("Error: first edge not found in UpdateNeighbors\n");
      exit(1);
    }
    // once starting point is established, loop over other edges to sort
    for (int i_entry{1}; i_entry < vrt.n_edges_; i_entry++) {
      Edge *prev_edge{edges_ordered[i_entry - 1]};
      Vertex *prev_neighb{neighbs_ordered[i_entry - 1]};
      for (int i_edge{0}; i_edge < vrt.n_edges_; i_edge++) {
        Edge *next_edge{vrt.edges_[i_edge]};
        // disregard self-comparisons (free life advice 4 ya)
        if (prev_edge == next_edge) {
          continue;
        }
        Triangle *next_tri{nullptr};
        // SF TODO remove chicanery
        if (prev_edge == nullptr) {
          printf("no thx @ %i\n", i_entry);
          for (auto &&edge : vrt.edges_) {
            printf("  edge %zu\n", edge->i_);
          }
          printf("(%i entries total)\n", vrt.n_edges_);
          vrt.SetColor(2 * M_PI, draw_type::fixed);
          vrt.SetDiameter(2);
          do_not_pass_go_ = true;
          return;
          // exit(1);
        }
        for (auto &&tri_i : prev_edge->tris_) {
          for (auto &&tri_j : next_edge->tris_) {
            if (tri_i == tri_j) {
              if (next_tri != nullptr) {
                printf("can't have more than 1 match my guy\n");
                exit(1);
              }
              next_tri = tri_i;
            }
          }
        }
        // if triangle match wasn't found, these edges are not adjacent
        if (next_tri == nullptr) {
          continue;
        }
        // if triangle was found, need to check its ring-wise direction
        Vertex *next_neighb{next_edge->GetOtherEnd(&vrt)};
        // points from prev_neighb (j_prev) to vrt (i)
        double r_ij_prev[3];
        // ponts from next_neighb (j_next) to vrt (i);
        double r_ij_next[3];
        for (int i_dim{0}; i_dim < params_->n_dim; i_dim++) {
          r_ij_prev[i_dim] = vrt.pos_[i_dim] - prev_neighb->pos_[i_dim];
          r_ij_next[i_dim] = vrt.pos_[i_dim] - next_neighb->pos_[i_dim];
        }
        // for proper ccw order, r_ij_prev x r_ij_next should align with nhat
        // NOTE: this is what needs to be changed for concave objects
        double n_tri[3];
        cross_product(r_ij_prev, r_ij_next, n_tri, 3);
        if (dot_product(3, n_tri, next_tri->nhat_) > 0.0) {
          tris_ordered[i_entry - 1] = next_tri;
          edges_ordered[i_entry] = next_edge;
          neighbs_ordered[i_entry] = next_neighb;
          // since we are adding tris "behind" edges/neighbs, last tri must be added manually
          if (i_entry == vrt.n_edges_ - 1) {
            tris_ordered[i_entry] = next_edge->GetOtherTriangle(next_tri);
          }
        }
      }
    }
    for (int i_entry{0}; i_entry < vrt.n_neighbs_; i_entry++) {
      vrt.neighbs_[i_entry] = neighbs_ordered[i_entry];
      vrt.edges_[i_entry] = edges_ordered[i_entry];
      vrt.tris_[i_entry] = tris_ordered[i_entry];
    }
  }
}

void TriMesh::FlipEdges() {
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

void TriMesh::UpdateTriangles() {

  UpdateCentroid();
  for (auto &&tri : tris_) {
    double r1[3], r2[3]; // two edges from this triangle
    double r_origin[3];  // points from origin to center of triangle
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      r1[i_dim] = tri.vrts_[0]->pos_[i_dim] - tri.vrts_[1]->pos_[i_dim];
      r2[i_dim] = tri.vrts_[1]->pos_[i_dim] - tri.vrts_[2]->pos_[i_dim];
      r_origin[i_dim] = tri.GetCenterPos(i_dim) - centroid_[i_dim];
    }
    double n[3];
    cross_product(r1, r2, n, 3);
    if (dot_product(3, n, r_origin) < 0.0) {
      for (int i_dim{0}; i_dim < 3; i_dim++) {
        n[i_dim] *= -1;
      }
      if (dot_product(3, n, r_origin) < 0.0) {
        printf("what the fuc\n");
      }
    }
    double n_mag{0.0};
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      n_mag += SQR(n[i_dim]);
    }
    n_mag = sqrt(n_mag);
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      tri.nhat_[i_dim] = n[i_dim] / n_mag;
    }
  }
}

void TriMesh::ProjectToUnitSphere() {
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

void TriMesh::Draw(std::vector<graph_struct *> &graph_array) {
  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    vrts_[i_vrt].Draw(graph_array);
  }
  if (params_->draw_mindist) {
    for (int i{0}; i < f_mem_.size(); i++) {
      graph_array.push_back(&f_mem_[i]);
    }
  }
  if (params_->draw_centroid) {
    graph_array.push_back(&o_); // plots centroid pos
  }
}

void TriMesh::UpdateMesh() {
  // Zero out all forces
  for (auto &&vrt : vrts_) {
    vrt.ZeroForce();
  }
  // Update edge lengths (used to calculate triangle areas)
  for (auto &&edge : edges_) {
    edge.Update();
  }
  // Update triangle areas (used to calculate origin position)
  for (auto &&tri : tris_) {
    tri.UpdateArea();
  }
  // Update mesh origin (used to calculate triangle volumes)
  UpdateCentroid();
  // Update triangle volumes
  for (auto &&tri : tris_) {
    tri.UpdateVolume(centroid_);
  }
  // Update triangle nhats SF TODO make per-triangle nhat update
  UpdateTriangles();

  // mathematical jonnies from NAB used in calculating boundary forces
  for (int itri{0}; itri < tris_.size(); itri++) {
    Triangle *tri{&tris_[itri]};
    double aVector[3] = {0.0};
    double bVector[3] = {0.0};
    double cVector[3] = {0.0};

    // Compute rotation to a plane and rotated coordinates
    for (int i = 0; i < 3; ++i) {
      aVector[i] =
          tri->vrts_[1]->GetPosition()[i] - tri->vrts_[0]->GetPosition()[i];
      bVector[i] =
          tri->vrts_[2]->GetPosition()[i] - tri->vrts_[1]->GetPosition()[i];
    }
    cross_product(aVector, bVector, cVector, 3);
    double cLenInXY = sqrt(cVector[0] * cVector[0] + cVector[1] * cVector[1]);
    double cLength = sqrt(cVector[0] * cVector[0] + cVector[1] * cVector[1] +
                          cVector[2] * cVector[2]);

    double sinGam, cosGam, sinBet, cosBet;
    if (cLenInXY < 1.0e-10) {
      sinGam = 0.0;
      cosGam = 1.0;
    } else {
      sinGam = -cVector[1] / cLenInXY;
      cosGam = cVector[0] / cLenInXY;
    }
    sinBet = cLenInXY / cLength;
    cosBet = cVector[2] / cLength;

    tri->cosGamma_ = cosGam;
    tri->sinGamma_ = sinGam;
    tri->cosBeta_ = cosBet;
    tri->sinBeta_ = sinBet;

    // Calculate and save rotated coordiantes
    double zRotSum = 0;
    for (int ivert = 0; ivert < 3; ++ivert) {
      int iv = ivert; //tri->indVert[ivert][itri];
      tri->XYrot_[0][ivert] = (tri->vrts_[iv]->GetPosition()[0] * cosGam -
                               tri->vrts_[iv]->GetPosition()[1] * sinGam) *
                                  cosBet -
                              tri->vrts_[iv]->GetPosition()[2] * sinBet;
      tri->XYrot_[1][ivert] = tri->vrts_[iv]->GetPosition()[0] * sinGam +
                              tri->vrts_[iv]->GetPosition()[1] * cosGam;
      zRotSum = zRotSum +
                (tri->vrts_[iv]->GetPosition()[0] * cosGam -
                 tri->vrts_[iv]->GetPosition()[1] * sinGam) *
                    sinBet +
                tri->vrts_[iv]->GetPosition()[2] * cosBet;
    }
    tri->Zrot_ = zRotSum / 3.;
  }
}

void TriMesh::ApplyMembraneForces() {

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

void TriMesh::ApplyBoundaryForces() {

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

void TriMesh::WriteOutputs() {
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

void TriMesh::UpdatePositions() {

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
