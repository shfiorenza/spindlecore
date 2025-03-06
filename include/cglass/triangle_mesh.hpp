#ifndef _CGLASS_TRIANGLE_MESH_
#define _CGLASS_TRIANGLE_MESH_

// #include "common_libs.hpp"
// #include "definitions.hpp"
#include "minimum_distance.hpp"
#include "rng.hpp"
#include "site.hpp"

// TODO add param storage and sync site seeds

struct Triangle;
struct Edge;
struct Vertex : public Site {
  size_t i_{0}; // index in master vrts_ list
  // SF TODO link
  int seed{0};
  double pos_[3];

  int n_tris_ = 0;
  std::vector<Triangle *> tris_; // triangles this vertex is a part of
  int n_neighbs_ = 0;
  std::vector<Vertex *> neighbs_;
  int n_edges_ = 0; // should be replaced by n_neigbs when all is done
  std::vector<Edge *> edges_;

  Vertex() : Site(seed) {
    pos_[0] = pos_[1] = pos_[2] = 0;
    Site::SetPositionXYZ(0, 0, 0);
  }
  Vertex(double x, double y, double z) : Site(seed) {
    pos_[0] = x;
    pos_[1] = y;
    pos_[2] = z;
    Site::SetPositionXYZ(x, y, z);
  }
  friend bool operator==(const Vertex &lhs, const Vertex &rhs) {
    return (lhs.pos_[0] == rhs.pos_[0] and lhs.pos_[1] == rhs.pos_[1] and
            lhs.pos_[2] == rhs.pos_[2]);
  }
  friend bool operator!=(const Vertex &lhs, const Vertex &rhs) {
    return !(lhs == rhs);
  }
  void SetPos(const double *const new_pos) {
    pos_[0] = pos[0] = position_[0] = new_pos[0];
    pos_[1] = pos[1] = position_[1] = new_pos[1];
    pos_[2] = pos[2] = position_[2] = new_pos[2];
  }
};

struct Edge {
  size_t i_{0}; // index in master edges_ list
  bool just_flipped{false};

  double length_{0.0};
  double vector_[3]; // points from vrt 0 to vrt 1

  Vertex *vrts_[2];   // each endpoint
  Triangle *tris_[2]; // each adjacent triangle
  Edge() {}
  Edge(Vertex *start, Vertex *end) {
    vrts_[0] = start;
    vrts_[1] = end;
  }
  friend bool operator==(const Edge &lhs, const Edge &rhs) {
    return ((lhs.vrts_[0] == rhs.vrts_[0] and lhs.vrts_[1] == rhs.vrts_[1]) or
            (lhs.vrts_[0] == rhs.vrts_[1] and lhs.vrts_[1] == rhs.vrts_[0]));
  }
  friend bool operator!=(const Edge &lhs, const Edge &rhs) {
    return !(lhs == rhs);
  }
  bool Contains(Vertex *vrt) { return vrts_[0] == vrt or vrts_[1] == vrt; }
  Vertex *GetOtherEnd(Vertex *vrt) {
    return vrt == vrts_[0] ? vrts_[1] : vrts_[0];
  }
  Triangle *GetOtherTriangle(Triangle *tri) {
    return tri == tris_[0] ? tris_[1] : tris_[0];
  }
  void Update() {
    length_ = 0.0;
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      vector_[i_dim] = vrts_[1]->pos_[i_dim] - vrts_[0]->pos_[i_dim];
      length_ += SQR(vector_[i_dim]);
    }
    length_ = sqrt(length_);
  }
};

struct Triangle {

  size_t i_{0}; // index in master tris_ lit

  bool flipped_{false};
  double area_;
  double volume_;
  double nhat_[3];
  double color_[3];

  // euler angles
  double cosGamma_;
  double sinGamma_;
  double cosBeta_;
  double sinBeta_;
  double Zrot_;
  double XYrot_[2][3]; // [dim][i_vrt]

  Vertex *vrts_[3]{{}}; // vertices that compose this triangle
  Edge *edges_[3]{{}};  // edges that compose this triangle
  Triangle *neighbs_[3]{{}};

  Triangle(Vertex *v1, Vertex *v2, Vertex *v3) {
    vrts_[0] = v1;
    vrts_[1] = v2;
    vrts_[2] = v3;
    color_[0] = rand() % 255;
    color_[1] = rand() % 255;
    color_[2] = rand() % 255;
  }
  double GetCenterPos(int i_dim) {
    return (vrts_[0]->pos_[i_dim] + vrts_[1]->pos_[i_dim] +
            vrts_[2]->pos_[i_dim]) /
           3.0;
  }
  bool Contains(Vertex *vrt) {
    return vrts_[0] == vrt or vrts_[1] == vrt or vrts_[2] == vrt;
  }

  Vertex *GetOtherVertex(Vertex *vrt1, Vertex *vrt2) {
    if (vrts_[0] != vrt1 and vrts_[0] != vrt2) {
      return vrts_[0];
    } else if (vrts_[1] != vrt1 and vrts_[1] != vrt2) {
      return vrts_[1];
    } else if (vrts_[2] != vrt1 and vrts_[2] != vrt2) {
      return vrts_[2];
    } else {
      printf("Error finding other vertex in triangle %zu\n", i_);
      exit(1);
    }
  }
  Edge *GetEdge(Vertex *vrt1, Vertex *vrt2) {
    for (auto &&edge : edges_) {
      if (edge->Contains(vrt1) and edge->Contains(vrt2)) {
        return edge;
      }
    }
    printf("Error finding edge in triangle %zu\n", i_);
    exit(1);
  }
  void Update(double origin[]) {
    UpdateArea();
    UpdateVolume(origin);
  }
  void UpdateArea() {
    double a{edges_[0]->length_};
    double b{edges_[1]->length_};
    double c{edges_[2]->length_};
    double s{0.5 * (a + b + c)};
    area_ = sqrt(s * (s - a) * (s - b) * (s - c));
  }
  void UpdateVolume(double origin[]) {
    // update nhat
    // update volume
    double A[3];
    double B[3];
    double C[3];
    for (int i_dim{0}; i_dim < 3; i_dim++) {
      A[i_dim] = vrts_[0]->pos_[i_dim] - origin[i_dim];
      B[i_dim] = vrts_[1]->pos_[i_dim] - origin[i_dim];
      C[i_dim] = vrts_[2]->pos_[i_dim] - origin[i_dim];
    }
    double BxC[3];
    cross_product(B, C, BxC, 3);
    volume_ = std::fabs(dot_product(3, A, BxC) / 6.0);
  }
};

class TriMesh {

private:
  static const size_t n_edges_min_{3};  // true for any connected graph
  static const size_t n_edges_max_{10}; // arbitrary choice

  bool do_not_pass_go_{false};
  int i_datapoint_{0};

  FILE *forces_{nullptr};    // average force from each type of potential
  FILE *vertices_{nullptr};  // position (3D per vrt per step)
  FILE *adjacency_{nullptr}; // adjacency matrix  (2D per vrt per step)
  system_parameters *params_{nullptr};

  double f_avgs_[4]; // indices 0-4: tether, bend, area, vol

  double l_avg_{0.0};
  double gamma_{0.0};

  // params for radial force
  double kappa_B_{0.0};
  double l_max_{0.0};
  double l_min_{0.0};
  double l_c0_{0.0};
  double l_c1_{0.0};
  // params for bending force
  double kappa_{0.0};
  // params for area conservation force
  double kappa_l_{0.0};
  double A_prime_{0.0};
  // params for volume conservation force
  double kappa_v_{0.0};
  double V_prime_{0.0};
  RNG *rng_; // SF TODO link with system RNG
  MinimumDistance mindist_;

public:
  double r_sys_{0.0};
  double centroid_[3];

  std::vector<Object *> boundary_neighbs_;

  std::vector<Vertex> vrts_;
  std::vector<Triangle> tris_;
  std::vector<Edge> edges_;

  std::vector<graph_struct> f_mem_;
  graph_struct o_;

private:
  void SetParameters();
  void MakeIcosphere();
  void MakeIcosahedron();
  void DivideFaces();
  void ProjectToUnitSphere();
  void InitializeMesh();
  void FlipEdges();
  void UpdateCentroid();
  void UpdateTriangles();
  void UpdateNeighbors();
  void UpdateMesh();
  void ApplyMembraneForces();
  void ApplyBoundaryForces();

public:
  TriMesh() {}
  void Init(system_parameters *params);
  void Draw(std::vector<graph_struct *> &graph_array);
  void UpdatePositions();
  void WriteOutputs();
};

#endif