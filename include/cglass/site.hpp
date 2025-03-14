#ifndef _CGLASS_SITE_H_
#define _CGLASS_SITE_H_

#include "object.hpp"

class Bond; // Forward declaration
enum directed_type { OUTGOING, INCOMING, NONE };
typedef std::pair<Bond *, directed_type> directed_bond;
typedef std::vector<directed_bond>::iterator db_iterator;

// Sites, ie graph vertices
class Site : public Object {
protected:
  std::vector<directed_bond> bonds_;
  double tangent_[3]; // if one or two bonds, vector tangent to bonds at site
  double random_force_[3]; // random forces for filaments
  int n_bonds_;
  double theta_ = 0;
  double phi_ = 0;

public:
  Site(unsigned long seed);
  void AddBond(Bond *bond, directed_type dir);
  void SetPositionXYZ(double x, double y, double z) {
    position_[0] = x;
    position_[1] = y;
    position_[2] = z;
  }
  void Report();
  void ReportBonds();
  Bond *GetBond(int i);
  Bond *GetOtherBond(int bond_oid);
  void CalcTangent();
  void SetRandomForce(double *f_rand);
  void AddRandomForce();
  double const *const GetRandomForce();
  double const *const GetTangent();
  directed_bond GetDirectedBond(int i);
  directed_bond GetOutgoingBond();
  directed_bond GetOtherDirectedBond(int bond_oid);
  void RemoveOutgoingBonds();
  void RemoveBond(int bond_oid);
  virtual bool HasNeighbor(int other_oid);
  void WriteSpec(std::fstream &op);
  void ReadSpec(std::fstream &ip);
  void SetTheta(const double theta) { theta_ = theta; }
  void SetPhi(const double phi) { phi_ = phi; }
  const double GetTheta() const { return theta_; }
  const double GetPhi() const { return phi_; }
};

#endif // _CGLASS_SITE_H_
