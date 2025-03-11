#ifndef _CGLASS_TRIANGLE_MESH_BOUNDARY_
#define _CGLASS_TRIANGLE_MESH_BOUNDARY_

#include "triangle_mesh.hpp"

class TriMeshBoundary : public TriMesh {

private:
  int i_datapoint_{0};

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

  RNG *rng_;                 // SF TODO link with system RNG
  FILE *forces_{nullptr};    // average force from each type of potential
  FILE *vertices_{nullptr};  // position (3D per vrt per step)
  FILE *adjacency_{nullptr}; // adjacency matrix  (2D per vrt per step)

public:
  double r_sys_{0.0};

  std::vector<Object *> boundary_neighbs_;

  std::vector<graph_struct> f_mem_;

private:
  void SetParameters();
  void MakeIcosphere();
  void MakeIcosahedron();
  void ProjectToUnitSphere();
  void InitializeMesh();

  void FlipEdges();
  void ApplyMembraneForces();
  void ApplyBoundaryForces();

public:
  TriMeshBoundary() {}
  void Init(system_parameters *params);
  void WriteOutputs();
  void Draw(std::vector<graph_struct *> &graph_array);

  void UpdatePositions();
};

#endif