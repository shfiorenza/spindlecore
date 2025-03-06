#ifndef _CGLASS_CENTROSOME_SPECIES_H_
#define _CGLASS_CENTROSOME_SPECIES_H_

#include "centrosome.hpp"
#include "species.hpp"

class FilamentSpecies;
class RigidFilamentSpecies;

class CentrosomeSpecies : public Species<Centrosome, species_id::centrosome> {

private:
  arma::vec4 *q_;  // Quaternion describing current rotation of the SPB
  arma::mat33 *A_; // Rotation matrix describing SPB body-frame to fixed-frame
  arma::mat33 *mu_tb_;           // Translation mobility matrix
  arma::mat33 *mu_rb_;           // Rotation mobility matrix
  arma::mat33 *sqrt_mu_tb_;      // Sqrt translation mobility matrix
  arma::mat33 *sqrt_mu_rb_;      // Sqrt rotation mobility matrix
  arma::mat::fixed<4, 3> *bsub_; // Transform matrix for quaternion rotation

protected:
  bool midstep_;

public:
  CentrosomeSpecies(unsigned long seed);
  ~CentrosomeSpecies() {
    delete q_;
    delete A_;
    delete mu_tb_;
    delete mu_rb_;
    delete sqrt_mu_tb_;
    delete sqrt_mu_rb_;
    delete bsub_;
  }
  void Init(std::string spec_name, ParamsParser &parser);
  void PopMember();
  void AddMember();
  void Reserve();
  void ArrangeMembers() {} // blank for now
  std::string GetFilamentSpeciesName() {
    return sparams_.filament_species_name;
  }

  void AnchorFilaments(SpeciesBase *filas, bool flexible);

  void UpdatePositions();

  void ApplyInteractions() {
    // SF TODO put calculations from UpdatePositions() here?
    for (auto &&centro : members_) {
      centro.ApplyInteractions();
    }
  }

  double CalcNonMonotonicWallForce(double ne_ratio, double f0, double delta_r);

  double ComputeLagrangeCorrection(const arma::vec4 &q,
                                   const arma::vec4 &qtilde);

  void CreateRotationMatrixBodySpace(arma::mat33 &A, double *u, double *v,
                                     double *w);

  void QuaternionFromRotationMatrix(arma::vec4 &q, const arma::mat33 R);
  void RotationMatrixFromQuaternion(arma::mat33 &R, const arma::vec4 &q);

  void GetInteractors(std::vector<Object *> &ixors);
  void GetLastInteractors(std::vector<Object *> &ixors) {}
};

#endif
