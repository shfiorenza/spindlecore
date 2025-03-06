#ifndef _CGLASS_MINIMUM_DISTANCE_H_
#define _CGLASS_MINIMUM_DISTANCE_H_

#include "interaction.hpp"
#include "object.hpp"

class TriMesh;
class Triangle;

class MinimumDistance {
private:
  static int n_dim_;
  static int n_periodic_;
  static double *unit_cell_;
  static double boundary_cut2_;
  static space_struct *space_;

public:
  void PointPoint(double const *const r1, double const *const s1,
                  double const *const r2, double const *const s2, double *dr,
                  double *dr_mag2, double *midpoint);
  void PointCarrierLineInf(double *r_point, double *s_point, double *r_line,
                           double *s_line, double *u_line, double length,
                           double *dr, double *mu);
  void PointCarrierLine(double *r_point, double *s_point, double *r_line,
                        double *s_line, double *u_line, double length,
                        double *dr, double *r_contact, double *mu_ret);
  void Sphero(double const *const r_1, double const *const s_1,
              double const *const u_1, double const length_1,
              double const *const r_2, double const *const s_2,
              double const *const u_2, double const length_2, double *r_min,
              double *r_min_mag2, double *contact1, double *contact2);
  void SpheroDr(double *r_1, double *s_1, double *u_1, double length_1,
                double *r_2, double *s_2, double *u_2, double length_2,
                double *dr, double *r_min, double *r_min_mag2, double *lambda,
                double *mu);
  void SphereSphero(double const *const r_1, double const *const s_1,
                    double const *const r_2, double const *const s_2,
                    double const *const u_2, double const length_2,
                    double *r_min, double *r_min_mag2, double *contact2);
  void SpheroPlane(double *r_bond, double *u_bond, double length,
                   double *r_plane, double *n_plane, double *lambda,
                   double *r_min_mag2, double *r_min);
  void CarrierLines(double *r_1, double *s_1, double *u_1, double *r_2,
                    double *s_2, double *u_2, double *r_min, double *r_min_mag2,
                    double *lambda, double *mu);
  void PointSphereBC(double const *const r, double *dr, double *dr_mag2,
                     double buffer);
  void SpheroSphereBC(double const *const r, double const *const u,
                      double const length, double *r_contact, double *dr,
                      double *dr_mag2, double buffer);
  void PointBuddingBC(double const *const r, double *dr, double *dr_mag2,
                      double buffer);
  void SpheroBuddingBC(double const *const r, double const *const u,
                       double const length, double *dr, double *dr_mag2,
                       double *r_contact, double buffer);

  void PointWallBC(double const *const r, double const *const s, double *dr,
                   double *dr_mag2, double buffer);
  void SpheroWallBC(double const *const r, double const *const s,
                    double const *const u, double const length,
                    double *r_contact, double *dr, double *dr_mag2,
                    double buffer);

  // algorithms for triangulated meshes (below)
  int SpheroPolygon(TriMesh *polygon, double *r_1, double *s_1, double *u_1,
                    double length_1, double *rmin, double *rminmag2,
                    double *rcontact, double *mu);
  int SegmentToPolygon(TriMesh *polygon, double xStart, double yStart,
                       double zStart, double xEnd, double yEnd, double zEnd,
                       double *t, double *xRot, double *yRot, double *dist,
                       double *rcontact);
  void SegmentToTriangle(Triangle *tri, double xStart, double yStart,
                         double zStart, double xEnd, double yEnd, double zEnd,
                         double *t, double *xRot, double *yRot, double *dist);
  void SegmentDist(double xStart1, double yStart1, double zStart1, double xEnd1,
                   double yEnd1, double zEnd1, double xStart2, double yStart2,
                   double zStart2, double xEnd2, double yEnd2, double zEnd2,
                   double *t1, double *t2, double *dist);
  bool InsideTriangle(double *bx, double *by, double xTest, double yTest);
  // Returns the XOR between two booleans a, b
  inline bool b3dxor(bool a, bool b) { return ((a && !b) || (b && !a)); }
  // Inlined version of distance calc for ease
  // Function returns the parameter t for the point along the line
  // closest to x, y, z
  /**
 * @brief Inlined version of point distance to segment
 * @param (x) Point to find distance from
 * @param (a, b, c) Equation of line
 * @return parameter t for linear distance along line
 */
  inline double tri_point_to_line(double xStart, double yStart, double zStart,
                                  double a, double b, double c, double x,
                                  double y, double z) {
    return (a * (x - xStart) + b * (y - yStart) + c * (z - zStart)) /
           (a * a + b * b + c * c);
  }

  MinimumDistance() {}
  static void Init(space_struct *space, double boundary_cutoff_sq);
  void ObjectObject(Interaction &ix);
  bool CheckBoundaryInteraction(Interaction &ix);
  bool CheckOutsideBoundary(Object &o1);

  // XXX
  int GetNDim() { return n_dim_; }
};

#endif
