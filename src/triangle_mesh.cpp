#include <cglass/filament.hpp>
#include <cglass/triangle_mesh.hpp>
#include <unistd.h>

void TriMesh::Init() {}

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
  for (int i_edge{0}; i_edge < edges_.size(); i_edge++) {
    edges_[i_edge].i_ = i_edge;
    edges_[i_edge].Update();
  }
  // Initialize triangle indices and manually calculate areas (edges not assigned yet)
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
  }
  // Now that triangle areas have been updated, calculate centroid position
  UpdateCentroid();
  // Update neighbor lists and pointers for vertices, edges, and triangles
  UpdateNeighbors();
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

void TriMesh::Draw(std::vector<graph_struct *> &graph_array) {
  for (int i_vrt{0}; i_vrt < vrts_.size(); i_vrt++) {
    vrts_[i_vrt].Draw(graph_array);
  }
  if (params_->draw_centroid) {
    graph_array.push_back(&o_); // plots centroid pos
  }
}