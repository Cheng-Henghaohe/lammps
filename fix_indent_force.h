#ifdef FIX_CLASS
// clang-format off
FixStyle(indent/force,FixIndentForce);
// clang-format on
#else

#ifndef LMP_FIX_INDENT_FORCE_H
#define LMP_FIX_INDENT_FORCE_H

#include "../fix.h"

namespace LAMMPS_NS {

class FixIndentForce : public Fix {
 public:
  FixIndentForce(class LAMMPS *, int, char **);
  ~FixIndentForce() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

 private:
  int istyle, scaleflag, side;
  double k, k3;
  char *xstr, *ystr, *zstr, *rstr, *pstr, *fxstr, *fystr, *fzstr;
  int xvar, yvar, zvar, rvar, pvar;
  double fxvalue, fyvalue, fzvalue;
  double xvalue, yvalue, zvalue, rvalue, pvalue;
  int indenter_flag, planeside;
  double indenter[4], indenter_all[7];  // [energy,fx,fy,fz] and [energy,fx,fy,fz,x,y,z]
  int cdim, varflag;
  int ilevel_respa;

  char *rlostr, *rhistr, *lostr, *histr;
  int rlovar, rhivar, lovar, hivar;
  double rlovalue, rhivalue, lovalue, hivalue;

  // Force control parameters
  double target_fx, target_fy, target_fz;    // Target forces
  double tolerance;                          // Force tolerance
  int max_iterations;                        // Max iterations per timestep
  double adjustment_factor;                  // Position adjustment factor
  
  // History tracking for improved iteration
  double prev_fx_error, prev_fy_error, prev_fz_error;  // Previous iteration errors
  double prev_adjustment[3];                            // Previous position adjustments
  double prev_position[3];                             // Previous indenter position
  double prev_force[3];                                // Previous total forces
  bool first_iteration;                                // Flag for first iteration

  // Aggressive control parameters
  int slow_convergence_count;                          // Count of slow convergence iterations
  double step_multiplier;                              // Dynamic step size multiplier
  int oscillation_count;                               // Count consecutive oscillations

  // methods for argument parsing
  int geometry(int, char **);
  int options(int, char **);
  void force_arg(int, char **);
  void force_control_params(int, char **);

  // Force control methods
  void apply_indenter_forces(bool apply_to_atoms = true);
  bool force_converged(double fx_error, double fy_error, double fz_error);
  void adjust_position(double fx_error, double fy_error, double fz_error);
  void handle_initial_contact(double *global_forces);
  
  // Improved iteration algorithms
  double get_adaptive_factor(double error, double target, int iteration);
  bool detect_oscillation(double current_adj, double prev_adj);
  double estimate_gradient(double current_force, double prev_force, double current_pos, double prev_pos);
  double get_dynamic_limit(double error, double target, int iteration);
  double get_aggressive_multiplier(int iteration);

  // methods for conical indenter
  bool PointInsideCone(int, double *, double, double, double, double, double *);
  void DistanceExteriorPoint(int, double *, double, double, double, double, double &, double &,
                             double &);
  void DistanceInteriorPoint(int, double *, double, double, double, double, double &, double &,
                             double &);
  void point_on_line_segment(double *, double *, double *, double *);
  double closest(double *, double *, double *, double);
};

}    // namespace LAMMPS_NS

#endif
#endif