#ifndef _CGLASS_CENTROSOME_H_
#define _CGLASS_CENTROSOME_H_

#include "object.hpp"
// #include "filament.hpp"
// #include "receptor.hpp"
#include <armadillo>

// class Filament;
class Mesh;
class CentrosomeSpecies;

struct AnchorSite {
  double k_{0.0};
  double kr_{0.0};
  double r0_{0.0};
  double pos_[3];
  double pos_rel_[3];
  double u_[3];
  double u_rel_[3];
  Mesh *filament_;
  graph_struct g_;
};

// struct AnchorTether {

//   double k_{0.0};
//   double kr_{0.0};
//   double r_{0.0};
//   double r0_{0.0};
//   Filament *filament_;
//   AnchorSite *site_;
//   graph_struct g_;
// };

class Centrosome : public Object {
protected:
  int n_filaments_;

  double noise_tr_ = 0.0;
  double noise_rot_ = 0.0;
  double diffusion_ = 0.0;
  double diffusion_rot_ = 0.0;
  double body_frame_[6];

  // FIXME; introduce 'bioriented' and 'natural' lol
  static size_t i_spb_; // lol this is bad
  size_t index_;        // consequence of above badness

  double r_[3];
  double u_[3];
  double v_[3];
  double w_[3];

  std::vector<AnchorSite> anchors_;

  centrosome_parameters *sparams_;

  // FIXME this is probably unnecessary and a sloppy middleman fix
  friend CentrosomeSpecies;

public:
protected:
  void ApplyForcesTorques();
  void ApplyBoundaryForces();
  void InsertCentrosome();
  void GenerateAnchorSites();
  void SetDiffusion();
  void Translate();
  void Rotate();
  void Integrate();
  void RandomizeAnchorPosition(int i_fil);

public:
  Centrosome(unsigned long seed);
  void Init(centrosome_parameters *sparams);

  void SetPosition(double const *const new_pos) {}
  void UpdatePosition() {}

  void GetInteractors(std::vector<Object *> &ixors) {}
  double GetR(int i_dim) { return r_[i_dim]; }
  double GetU(int i_dim) { return u_[i_dim]; }
  double GetV(int i_dim) { return v_[i_dim]; }
  double GetW(int i_dim) { return w_[i_dim]; }

  void Draw(std::vector<graph_struct *> &graph_array);
  void ApplyInteractions();
  void UpdatePosition(double *r_new, double *u_new, double *v_new,
                      double *w_new);
};

#endif // _CGLASS_CENTROSOME_H_
