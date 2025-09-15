#include "fix_indent_force.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "math_extra.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "utils.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NONE, SPHERE, CYLINDER, PLANE, CONE };
enum { INSIDE, OUTSIDE };

/* ---------------------------------------------------------------------- */

FixIndentForce::FixIndentForce(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), xstr(nullptr), ystr(nullptr), zstr(nullptr), rstr(nullptr), pstr(nullptr),
    rlostr(nullptr), rhistr(nullptr), lostr(nullptr), histr(nullptr), fxstr(nullptr),
    fystr(nullptr), fzstr(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix indent/force", error);

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 6;    // fx, fy, fz, x, y, z (or p for plane)
  energy_global_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  k = utils::numeric(FLERR, arg[3], false, lmp);
  if (k < 0.0) error->all(FLERR, "Illegal fix indent/force force constant: {}", k);
  k3 = k / 3.0;

  // Parse arguments: geometry, force args, options
  istyle = NONE;
  int iarg = geometry(narg - 4, &arg[4]) + 4;
  force_arg(3, &arg[iarg]);    // Force args are always 3 parameters (fx, fy, fz)
  iarg += 3;
  if (iarg < narg) options(narg - iarg, &arg[iarg]);

  // Set default force control parameters
  tolerance = 0.01;        // 1% tolerance
  max_iterations = 1000;    // Maximum iterations per timestep (increased for better convergence)
  adjustment_factor =
      0.01;    // Position adjustment factor 控制每次迭代中压头位置调整的步长，但未使用

  // Initialize history tracking variables
  // prev_fx_error / prev_fy_error / prev_fz_error:
  //   记录上一迭代在 x/y/z 方向上的力误差，用于收敛判断与步长自适应
  prev_fx_error = prev_fy_error = prev_fz_error = 0.0;

  // prev_adjustment[0..2]:
  //   记录上一迭代对压头位置在 x/y/z 方向的调整量，用于检测振荡并施加阻尼
  prev_adjustment[0] = prev_adjustment[1] = prev_adjustment[2] = 0.0;

  // prev_position[0..2]:
  //   记录上一迭代使用的压头参考位置（x/y/z），用于估算力-位移梯度（数值导数）
  prev_position[0] = prev_position[1] = prev_position[2] = 0.0;

  // prev_force[0..2]:
  //   记录上一迭代得到的合力分量（x/y/z），配合 prev_position 进行梯度估计与牛顿步
  prev_force[0] = prev_force[1] = prev_force[2] = 0.0;

  // first_iteration:
  //   标记当前时间步是否为第一次迭代；首迭代时关闭振荡检测/牛顿加速等历史相关策略
  first_iteration = true;

  // Initialize aggressive control parameters
  slow_convergence_count = 0;    // 慢收敛计数器，用于检测长时间未收敛
  step_multiplier = 1.0;         // 动态步长倍数，随未收敛时间递增
  oscillation_count = 0;         // 连续振荡计数器，用于智能振荡处理

  // Setup scaling
  // 说明：
  // - scaleflag 由 options() 的 "units" 选项设置：units lattice -> scaleflag=1，units box -> scaleflag=0
  // - 当使用 units lattice 时，输入的几何参数（球心坐标、半径、平面位置等）以晶格单位给出，
  //   需按晶格常数 xlattice/ylattice/zlattice 转换到盒子坐标单位；否则不缩放（因已是盒子单位）
  // - 后续会根据几何体类型（sphere/plane/…）在未绑定变量的情形下用这些缩放因子对 x/y/z/r/p 等进行一次性缩放
  const double xscale{scaleflag ? domain->lattice->xlattice
                                : 1.0};    // x 方向缩放因子（晶格->盒子），否则为 1
  const double yscale{scaleflag ? domain->lattice->ylattice
                                : 1.0};    // y 方向缩放因子（晶格->盒子），否则为 1
  const double zscale{scaleflag ? domain->lattice->zlattice
                                : 1.0};    // z 方向缩放因子（晶格->盒子），否则为 1

  // Apply scaling factors to geometry
  if (istyle == SPHERE || istyle == CYLINDER) {
    if (!xstr) xvalue *= xscale;
    if (!ystr) yvalue *= yscale;
    if (!zstr) zvalue *= zscale;
    if (!rstr) rvalue *= xscale;
  } else if (istyle == CONE) {
    if (!xstr) xvalue *= xscale;
    if (!ystr) yvalue *= yscale;
    if (!zstr) zvalue *= zscale;

    double scaling_factor = 1.0;
    switch (cdim) {
      case 0:
        scaling_factor = xscale;
        break;
      case 1:
        scaling_factor = yscale;
        break;
      case 2:
        scaling_factor = zscale;
        break;
    }

    if (!rlostr) rlovalue *= scaling_factor;
    if (!rhistr) rhivalue *= scaling_factor;
    if (!lostr) lovalue *= scaling_factor;
    if (!histr) hivalue *= scaling_factor;
  } else if (istyle == PLANE) {
    if (cdim == 0 && !pstr)
      pvalue *= xscale;
    else if (cdim == 1 && !pstr)
      pvalue *= yscale;
    else if (cdim == 2 && !pstr)
      pvalue *= zscale;
  } else {
    error->all(FLERR, "Unknown fix indent/force geometry: {}", istyle);
  }

  varflag = 0;
  if (xstr || ystr || zstr || rstr || pstr || rlostr || rhistr || lostr || histr || fxstr ||
      fystr || fzstr)
    varflag = 1;

  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;
  for (int i = 0; i < 7; i++) indenter_all[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixIndentForce::~FixIndentForce()
{
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] rstr;
  delete[] pstr;
  delete[] rlostr;
  delete[] rhistr;
  delete[] lostr;
  delete[] histr;
  delete[] fxstr;
  delete[] fystr;
  delete[] fzstr;
}

/* ---------------------------------------------------------------------- */

int FixIndentForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIndentForce::init()
{
  // Initialize position variables
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) error->all(FLERR, "Variable {} for fix indent/force does not exist", xstr);
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR, "Variable {} for fix indent/force is invalid style", xstr);
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all(FLERR, "Variable {} for fix indent/force does not exist", ystr);
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR, "Variable {} for fix indent/force is invalid style", ystr);
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all(FLERR, "Variable {} for fix indent/force does not exist", zstr);
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR, "Variable {} for fix indent/force is invalid style", zstr);
  }
  if (rstr) {
    rvar = input->variable->find(rstr);
    if (rvar < 0) error->all(FLERR, "Variable {} for fix indent/force does not exist", rstr);
    if (!input->variable->equalstyle(rvar))
      error->all(FLERR, "Variable {} for fix indent/force is invalid style", rstr);
  }
  if (pstr) {
    pvar = input->variable->find(pstr);
    if (pvar < 0) error->all(FLERR, "Variable {} for fix indent/force does not exist", pstr);
    if (!input->variable->equalstyle(pvar))
      error->all(FLERR, "Variable {} for fix indent/force is invalid style", pvar);
  }

  // Initialize cone variables
  if (rlostr) {
    rlovar = input->variable->find(rlostr);
    if (rlovar < 0) error->all(FLERR, "Variable {} for fix indent/force does not exist", rlostr);
    if (!input->variable->equalstyle(rlovar))
      error->all(FLERR, "Variable {} for fix indent/force is invalid style", rlostr);
  }
  if (rhistr) {
    rhivar = input->variable->find(rhistr);
    if (rhivar < 0) error->all(FLERR, "Variable {} for fix indent/force does not exist", rhistr);
    if (!input->variable->equalstyle(rhivar))
      error->all(FLERR, "Variable {} for fix indent/force is invalid style", rhivar);
  }
  if (lostr) {
    lovar = input->variable->find(lostr);
    if (lovar < 0) error->all(FLERR, "Variable {} for fix indent/force does not exist", lostr);
    if (!input->variable->equalstyle(lovar))
      error->all(FLERR, "Variable {} for fix indent/force is invalid style", lovar);
  }
  if (histr) {
    hivar = input->variable->find(histr);
    if (hivar < 0) error->all(FLERR, "Variable {} for fix indent/force does not exist", histr);
    if (!input->variable->equalstyle(hivar))
      error->all(FLERR, "Variable {} for fix indent/force is invalid style", hivar);
  }

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixIndentForce::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixIndentForce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIndentForce::post_force(int /*vflag*/)
{
  // Reset iteration tracking for new timestep
  first_iteration = true;

  // Reset aggressive control parameters for new timestep
  slow_convergence_count = 0;
  step_multiplier = 1.0;
  oscillation_count = 0;

  // Get target forces (initialize force variables first if needed)
  if (fxstr) {
    int fxvar = input->variable->find(fxstr);
    if (fxvar >= 0) target_fx = input->variable->compute_equal(fxvar);
  } else {
    target_fx = fxvalue;
  }
  if (fystr) {
    int fyvar = input->variable->find(fystr);
    if (fyvar >= 0) target_fy = input->variable->compute_equal(fyvar);
  } else {
    target_fy = fyvalue;
  }
  if (fzstr) {
    int fzvar = input->variable->find(fzstr);
    if (fzvar >= 0) target_fz = input->variable->compute_equal(fzvar);
  } else {
    target_fz = fzvalue;
  }

  // Wrap variable evaluations with clear/add
  if (varflag) modify->clearstep_compute();

  indenter_flag = 0;

  // Force control iteration loop - adjust indenter position until target load is achieved
  int iteration = 0;
  double fx_error, fy_error, fz_error;
  bool converged = false;

  do {
    iteration++;    // Increment at start so we count all attempts

    // Calculate indenter forces only (don't apply to atoms during iteration)
    apply_indenter_forces(false);

    // MPI reduce indenter forces across all processors
    double indenter_global[4];
    MPI_Allreduce(indenter, indenter_global, 4, MPI_DOUBLE, MPI_SUM, world);

    // Calculate force errors using global forces
    fx_error = target_fx - indenter_global[1];
    fy_error = target_fy - indenter_global[2];
    fz_error = target_fz - indenter_global[3];

    // Check convergence
    converged = force_converged(fx_error, fy_error, fz_error);

    // Track slow convergence for aggressive stepping
    if (!converged) {
      slow_convergence_count++;

      // Increase step multiplier for long-term non-convergence
      if (slow_convergence_count > 50) {
        step_multiplier = fmin(3.0, 1.0 + (slow_convergence_count - 50) * 0.02);
      }

      // Handle initial contact if needed
      handle_initial_contact(indenter_global);
    }

    // Always adjust position based on errors to maintain target load
    adjust_position(fx_error, fy_error, fz_error);

  } while (!converged && iteration < max_iterations);

  // Final force application - now apply forces to atoms with converged position
  apply_indenter_forces(true);

  // Force control completed - particles now have target force applied

  // Force control completed silently - forces can be accessed via thermo_style custom

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixIndentForce::apply_indenter_forces(bool apply_to_atoms)
{
  // Reset indenter forces
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;

  // Current indenter center
  double ctr[3] = {xvalue, yvalue, zvalue};
  if (xstr) {
    ctr[0] = input->variable->compute_equal(xvar);
    xvalue = ctr[0];  // Update xvalue for position output
  }
  if (ystr) {
    ctr[1] = input->variable->compute_equal(yvar);
    yvalue = ctr[1];  // Update yvalue for position output
  }
  if (zstr) {
    ctr[2] = input->variable->compute_equal(zvar);
    zvalue = ctr[2];  // Update zvalue for position output
  }

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delx, dely, delz, r, dr, fmag, fx, fy, fz;

  // Apply forces based on indenter geometry (same as original code)
  if (istyle == SPHERE) {
    domain->remap(ctr);
    double radius = rstr ? input->variable->compute_equal(rvar) : rvalue;
    if (radius < 0.0) error->all(FLERR, "Illegal fix indent/force sphere radius: {}", radius);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        delx = x[i][0] - ctr[0];
        dely = x[i][1] - ctr[1];
        delz = x[i][2] - ctr[2];
        domain->minimum_image(FLERR, delx, dely, delz);
        r = sqrt(delx * delx + dely * dely + delz * delz);
        if (side == OUTSIDE) {
          dr = r - radius;
          fmag = k * dr * dr;
        } else {
          dr = radius - r;
          fmag = -k * dr * dr;
        }
        if (dr >= 0.0) continue;
        fx = delx * fmag / r;
        fy = dely * fmag / r;
        fz = delz * fmag / r;
        
        // Only apply forces to atoms if requested (after convergence)
        if (apply_to_atoms) {
          f[i][0] += fx;
          f[i][1] += fy;
          f[i][2] += fz;
        }
        
        indenter[0] -= k3 * dr * dr * dr;
        indenter[1] -= fx;
        indenter[2] -= fy;
        indenter[3] -= fz;
      }
  } else if (istyle == PLANE) {
    // planar indenter
    double plane{pstr ? input->variable->compute_equal(pvar) : pvalue};
    if (pstr) pvalue = plane;  // Update pvalue for position output

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dr = planeside * (plane - x[i][cdim]);
        if (dr >= 0.0) continue;
        fmag = -planeside * k * dr * dr;

        // Only apply forces to atoms if requested (after convergence)
        if (apply_to_atoms) {
          f[i][cdim] += fmag;
        }

        indenter[0] -= k3 * dr * dr * dr;
        indenter[cdim + 1] -= fmag;
      }
  }
}

/* ---------------------------------------------------------------------- */

bool FixIndentForce::force_converged(double fx_error, double fy_error, double fz_error)
{
  // Only check convergence for directions with non-zero target forces
  // Zero target force means no control in that direction

  bool fx_converged = true;  // Default to converged for uncontrolled directions
  bool fy_converged = true;
  bool fz_converged = true;

  // Check x-direction only if target force is non-zero
  if (fabs(target_fx) > 1e-10) {
    double fx_tol = tolerance * fabs(target_fx);
    fx_converged = (fabs(fx_error) <= fx_tol);
  }

  // Check y-direction only if target force is non-zero
  if (fabs(target_fy) > 1e-10) {
    double fy_tol = tolerance * fabs(target_fy);
    fy_converged = (fabs(fy_error) <= fy_tol);
  }

  // Check z-direction only if target force is non-zero
  if (fabs(target_fz) > 1e-10) {
    double fz_tol = tolerance * fabs(target_fz);
    fz_converged = (fabs(fz_error) <= fz_tol);
  }

  return (fx_converged && fy_converged && fz_converged);
}

/* ---------------------------------------------------------------------- */

void FixIndentForce::adjust_position(double fx_error, double fy_error, double fz_error)
{
  // Store current position and forces for gradient estimation
  double current_pos[3] = {xvalue, yvalue, zvalue};
  if (istyle == PLANE) { current_pos[0] = current_pos[1] = current_pos[2] = pvalue; }

  double current_forces[3] = {target_fx - fx_error, target_fy - fy_error, target_fz - fz_error};

  // Adjust indenter position based on force errors with improved algorithms
  if (istyle == SPHERE) {
    // For sphere, adjust center position
    // Only adjust position if target force is non-zero in that direction
    if (!xstr && fabs(target_fx) > 1e-10) {
      double adaptive_factor = get_adaptive_factor(fx_error, target_fx, slow_convergence_count);
      double adjustment = adaptive_factor * fx_error * step_multiplier;

      // Smart oscillation handling - reduce damping for persistent convergence issues
      bool is_oscillating = false;
      if (!first_iteration && detect_oscillation(adjustment, prev_adjustment[0])) {
        oscillation_count++;
        is_oscillating = true;

        // Reduce oscillation damping if we've been oscillating for too long
        double damping_factor = (oscillation_count > 10) ? 0.8 : 0.5;
        adjustment *= damping_factor;
      } else {
        oscillation_count = 0;  // Reset oscillation count if not oscillating
      }

      // Try gradient-based acceleration (more aggressive)
      if (!first_iteration && !is_oscillating) {
        double gradient =
            estimate_gradient(current_forces[0], prev_force[0], current_pos[0], prev_position[0]);
        if (fabs(gradient) > 1e-10) {
          double newton_step = fx_error / gradient;
          // Use more aggressive Newton scaling, especially for slow convergence
          double newton_scale = (slow_convergence_count > 20) ? 1.2 : 0.8;
          if (fabs(newton_step * newton_scale) < fabs(adjustment) * 2.0) {
            adjustment = newton_scale * newton_step;
          }
        }
      }

      // Dynamic adjustment limit based on error magnitude (more aggressive)
      double max_adjustment = get_dynamic_limit(fx_error, target_fx, slow_convergence_count);
      if (fabs(adjustment) > max_adjustment) {
        adjustment = (adjustment > 0) ? max_adjustment : -max_adjustment;
      }

      xvalue += adjustment;
      prev_adjustment[0] = adjustment;
    }

    if (!ystr && fabs(target_fy) > 1e-10) {
      double adaptive_factor = get_adaptive_factor(fy_error, target_fy, slow_convergence_count);
      double adjustment = adaptive_factor * fy_error * step_multiplier;

      // Smart oscillation handling
      bool is_oscillating = false;
      if (!first_iteration && detect_oscillation(adjustment, prev_adjustment[1])) {
        oscillation_count++;
        is_oscillating = true;
        double damping_factor = (oscillation_count > 10) ? 0.8 : 0.5;
        adjustment *= damping_factor;
      } else {
        oscillation_count = 0;
      }

      // Aggressive gradient acceleration
      if (!first_iteration && !is_oscillating) {
        double gradient =
            estimate_gradient(current_forces[1], prev_force[1], current_pos[1], prev_position[1]);
        if (fabs(gradient) > 1e-10) {
          double newton_step = fy_error / gradient;
          double newton_scale = (slow_convergence_count > 20) ? 1.2 : 0.8;
          if (fabs(newton_step * newton_scale) < fabs(adjustment) * 2.0) {
            adjustment = newton_scale * newton_step;
          }
        }
      }

      // More aggressive dynamic limit
      double max_adjustment = get_dynamic_limit(fy_error, target_fy, slow_convergence_count);
      if (fabs(adjustment) > max_adjustment) {
        adjustment = (adjustment > 0) ? max_adjustment : -max_adjustment;
      }

      yvalue += adjustment;
      prev_adjustment[1] = adjustment;
    }

    if (!zstr && fabs(target_fz) > 1e-10) {
      double adaptive_factor = get_adaptive_factor(fz_error, target_fz, slow_convergence_count);
      double adjustment = adaptive_factor * fz_error * step_multiplier;

      // Smart oscillation handling
      bool is_oscillating = false;
      if (!first_iteration && detect_oscillation(adjustment, prev_adjustment[2])) {
        oscillation_count++;
        is_oscillating = true;
        double damping_factor = (oscillation_count > 10) ? 0.8 : 0.5;
        adjustment *= damping_factor;
      } else {
        oscillation_count = 0;
      }

      // Aggressive gradient acceleration
      if (!first_iteration && !is_oscillating) {
        double gradient =
            estimate_gradient(current_forces[2], prev_force[2], current_pos[2], prev_position[2]);
        if (fabs(gradient) > 1e-10) {
          double newton_step = fz_error / gradient;
          double newton_scale = (slow_convergence_count > 20) ? 1.2 : 0.8;
          if (fabs(newton_step * newton_scale) < fabs(adjustment) * 2.0) {
            adjustment = newton_scale * newton_step;
          }
        }
      }

      // More aggressive dynamic limit
      double max_adjustment = get_dynamic_limit(fz_error, target_fz, slow_convergence_count);
      if (fabs(adjustment) > max_adjustment) {
        adjustment = (adjustment > 0) ? max_adjustment : -max_adjustment;
      }

      zvalue += adjustment;
      prev_adjustment[2] = adjustment;
    }

  } else if (istyle == PLANE) {
    // For plane, adjust position along the normal direction
    if (!pstr) {
      double total_error = 0.0;
      double target_force = 0.0;

      if (cdim == 0) {
        total_error = fx_error;
        target_force = target_fx;
      } else if (cdim == 1) {
        total_error = fy_error;
        target_force = target_fy;
      } else if (cdim == 2) {
        total_error = fz_error;
        target_force = target_fz;
      }

      double adaptive_factor = get_adaptive_factor(total_error, target_force, slow_convergence_count);
      double adjustment = adaptive_factor * total_error * planeside * step_multiplier;

      // Smart oscillation handling
      bool is_oscillating = false;
      if (!first_iteration && detect_oscillation(adjustment, prev_adjustment[cdim])) {
        oscillation_count++;
        is_oscillating = true;
        double damping_factor = (oscillation_count > 10) ? 0.8 : 0.5;
        adjustment *= damping_factor;
      } else {
        oscillation_count = 0;
      }

      // Aggressive gradient acceleration
      if (!first_iteration && !is_oscillating) {
        double current_plane_force = target_force - total_error;
        double gradient =
            estimate_gradient(current_plane_force, prev_force[cdim], pvalue, prev_position[cdim]);
        if (fabs(gradient) > 1e-10) {
          double newton_step = total_error / gradient;
          double newton_scale = (slow_convergence_count > 20) ? 1.2 : 0.8;
          if (fabs(newton_step * newton_scale * planeside) < fabs(adjustment) * 2.0) {
            adjustment = newton_scale * newton_step * planeside;
          }
        }
      }

      // More aggressive dynamic limit
      double max_adjustment = get_dynamic_limit(total_error, target_force, slow_convergence_count);
      if (fabs(adjustment) > max_adjustment) {
        adjustment = (adjustment > 0) ? max_adjustment : -max_adjustment;
      }

      pvalue += adjustment;
      prev_adjustment[cdim] = adjustment;
    }
  }

  // Update history for next iteration
  prev_position[0] = current_pos[0];
  prev_position[1] = current_pos[1];
  prev_position[2] = current_pos[2];
  prev_force[0] = current_forces[0];
  prev_force[1] = current_forces[1];
  prev_force[2] = current_forces[2];
  prev_fx_error = fx_error;
  prev_fy_error = fy_error;
  prev_fz_error = fz_error;
  first_iteration = false;
}

/* ---------------------------------------------------------------------- */

void FixIndentForce::handle_initial_contact(double *global_forces)
{
  // Handle initial contact when current forces are near zero but target is not
  double current_force_mag =
      sqrt(global_forces[1] * global_forces[1] + global_forces[2] * global_forces[2] +
           global_forces[3] * global_forces[3]);
  double target_force_mag =
      sqrt(target_fx * target_fx + target_fy * target_fy + target_fz * target_fz);

  if (current_force_mag < 0.01 * target_force_mag && target_force_mag > 1e-10) {
    // Move indenter towards contact
    if (istyle == SPHERE) {
      // Calculate global center of mass with MPI consistency
      double cm_local[3] = {0.0, 0.0, 0.0};
      double cm_global[3] = {0.0, 0.0, 0.0};
      int count_local = 0, count_global = 0;
      double **x = atom->x;
      int *mask = atom->mask;
      int nlocal = atom->nlocal;

      // Calculate local center of mass
      for (int i = 0; i < nlocal; i++) {
        if (mask[i] & groupbit) {
          cm_local[0] += x[i][0];
          cm_local[1] += x[i][1];
          cm_local[2] += x[i][2];
          count_local++;
        }
      }

      // MPI reduce to get global center of mass
      MPI_Allreduce(cm_local, cm_global, 3, MPI_DOUBLE, MPI_SUM, world);
      MPI_Allreduce(&count_local, &count_global, 1, MPI_INT, MPI_SUM, world);

      if (count_global > 0) {
        cm_global[0] /= count_global;
        cm_global[1] /= count_global;
        cm_global[2] /= count_global;

        double move_factor = 0.1;

        // Move indenter ONLY in directions where target force is non-zero
        // Move direction depends on force sign (positive=compression, negative=tension)

        // X direction: only move if target_fx != 0
        if (!xstr && fabs(target_fx) > 1e-10) {
          double direction = (target_fx > 0) ? 1.0 : -1.0;    // positive force = move toward cm
          xvalue += direction * move_factor * (cm_global[0] - xvalue);
        }

        // Y direction: only move if target_fy != 0
        if (!ystr && fabs(target_fy) > 1e-10) {
          double direction = (target_fy > 0) ? 1.0 : -1.0;    // positive force = move toward cm
          yvalue += direction * move_factor * (cm_global[1] - yvalue);
        }

        // Z direction: only move if target_fz != 0
        if (!zstr && fabs(target_fz) > 1e-10) {
          double direction = (target_fz > 0) ? 1.0 : -1.0;    // positive force = move toward cm
          zvalue += direction * move_factor * (cm_global[2] - zvalue);
        }
      }

    } else if (istyle == PLANE) {
      // For plane indenter, only move along the normal direction
      double target_force_plane = 0.0;

      if (cdim == 0)
        target_force_plane = target_fx;
      else if (cdim == 1)
        target_force_plane = target_fy;
      else if (cdim == 2)
        target_force_plane = target_fz;

      // Only move if target force in plane direction is non-zero
      if (!pstr && fabs(target_force_plane) > 1e-10) {
        // Calculate global center of mass in plane direction
        double cm_local_sum = 0.0, cm_global_sum = 0.0;
        int count_local = 0, count_global = 0;
        double **x = atom->x;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;

        for (int i = 0; i < nlocal; i++) {
          if (mask[i] & groupbit) {
            cm_local_sum += x[i][cdim];
            count_local++;
          }
        }

        MPI_Allreduce(&cm_local_sum, &cm_global_sum, 1, MPI_DOUBLE, MPI_SUM, world);
        MPI_Allreduce(&count_local, &count_global, 1, MPI_INT, MPI_SUM, world);

        if (count_global > 0) {
          double cm_plane = cm_global_sum / count_global;
          double move_factor = 0.1;

          // Move direction depends on force sign and planeside
          // positive force means compression (plane pushes toward atoms)
          // planeside: -1 (lo) or +1 (hi)
          double direction = (target_force_plane > 0) ? planeside : -planeside;
          pvalue += direction * move_factor * fabs(cm_plane - pvalue) * 0.1;
        }
      }
    }
  }
}

// Include other geometry implementations and utility functions here
// (copying from original WANG implementation)

/* ---------------------------------------------------------------------- */

void FixIndentForce::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIndentForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixIndentForce::compute_scalar()
{
  if (indenter_flag == 0) {
    MPI_Allreduce(indenter, indenter_all, 4, MPI_DOUBLE, MPI_SUM, world);
    indenter_flag = 1;
  }

  // Always update position information (no MPI needed as it's same on all processors)
  indenter_all[4] = xvalue;    // x position
  indenter_all[5] = yvalue;    // y position
  if (istyle == PLANE) {
    indenter_all[6] = pvalue;    // plane position
  } else {
    indenter_all[6] = zvalue;    // z position
  }

  return indenter_all[0];
}

/* ---------------------------------------------------------------------- */

double FixIndentForce::compute_vector(int n)
{
  if (indenter_flag == 0) {
    MPI_Allreduce(indenter, indenter_all, 4, MPI_DOUBLE, MPI_SUM, world);
    indenter_flag = 1;
  }

  // Always update position information (no MPI needed as it's same on all processors)
  indenter_all[4] = xvalue;    // x position
  indenter_all[5] = yvalue;    // y position
  if (istyle == PLANE) {
    indenter_all[6] = pvalue;    // plane position
  } else {
    indenter_all[6] = zvalue;    // z position
  }

  // Return vector components: [fx, fy, fz, x, y, z/p]
  if (n < 3) {
    return indenter_all[n + 1];    // Forces: fx, fy, fz
  } else {
    return indenter_all[n + 1];    // Positions: x, y, z/p
  }
}

/* ----------------------------------------------------------------------
   parse input args for geometry of indenter
------------------------------------------------------------------------- */

int FixIndentForce::geometry(int narg, char **arg)
{
  if (narg < 0) utils::missing_cmd_args(FLERR, "fix indent/force", error);

  xstr = ystr = zstr = rstr = pstr = nullptr;
  xvalue = yvalue = zvalue = rvalue = pvalue = 0.0;

  // sphere
  if (strcmp(arg[0], "sphere") == 0) {
    if (5 > narg) utils::missing_cmd_args(FLERR, "fix indent/force sphere", error);

    if (utils::strmatch(arg[1], "^v_")) {
      xstr = utils::strdup(arg[1] + 2);
    } else
      xvalue = utils::numeric(FLERR, arg[1], false, lmp);
    if (utils::strmatch(arg[2], "^v_")) {
      ystr = utils::strdup(arg[2] + 2);
    } else
      yvalue = utils::numeric(FLERR, arg[2], false, lmp);
    if (utils::strmatch(arg[3], "^v_")) {
      zstr = utils::strdup(arg[3] + 2);
    } else
      zvalue = utils::numeric(FLERR, arg[3], false, lmp);
    if (utils::strmatch(arg[4], "^v_")) {
      rstr = utils::strdup(arg[4] + 2);
    } else
      rvalue = utils::numeric(FLERR, arg[4], false, lmp);

    istyle = SPHERE;
    return 5;
  }

  // plane
  if (strcmp(arg[0], "plane") == 0) {
    if (4 > narg) utils::missing_cmd_args(FLERR, "fix indent/force plane", error);
    if (strcmp(arg[1], "x") == 0)
      cdim = 0;
    else if (strcmp(arg[1], "y") == 0)
      cdim = 1;
    else if (strcmp(arg[1], "z") == 0)
      cdim = 2;
    else
      error->all(FLERR, "Unknown fix indent/force plane argument: {}", arg[1]);

    if (utils::strmatch(arg[2], "^v_")) {
      pstr = utils::strdup(arg[2] + 2);
    } else
      pvalue = utils::numeric(FLERR, arg[2], false, lmp);

    if (strcmp(arg[3], "lo") == 0)
      planeside = -1;
    else if (strcmp(arg[3], "hi") == 0)
      planeside = 1;
    else
      error->all(FLERR, "Unknown fix indent/force plane argument: {}", arg[3]);
    istyle = PLANE;
    return 4;
  }

  // invalid istyle arg
  error->all(FLERR, "Unknown fix indent/force geometry: {}", arg[0]);
  return 0;
}

/* ----------------------------------------------------------------------
   parse optional input args
------------------------------------------------------------------------- */

int FixIndentForce::options(int narg, char **arg)
{
  scaleflag = 1;
  side = OUTSIDE;

  int iarg = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "units") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix indent/force units", error);
      if (strcmp(arg[iarg + 1], "box") == 0)
        scaleflag = 0;
      else if (strcmp(arg[iarg + 1], "lattice") == 0)
        scaleflag = 1;
      else
        error->all(FLERR, "Unknown fix indent/force units argument: {}", arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "side") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix indent/force side", error);
      if (strcmp(arg[iarg + 1], "in") == 0)
        side = INSIDE;
      else if (strcmp(arg[iarg + 1], "out") == 0)
        side = OUTSIDE;
      else
        error->all(FLERR, "Unknown fix indent/force side argument: {}", arg[iarg + 1]);
      iarg += 2;

    } else
      error->all(FLERR, "Unknown fix indent/force argument: {}", arg[iarg]);
  }
  return 2;
}

/* ----------------------------------------------------------------------
   parse force arguments
------------------------------------------------------------------------- */

void FixIndentForce::force_arg(int narg, char **arg)
{
  if (narg < 3) utils::missing_cmd_args(FLERR, "fix indent/force", error);
  fxstr = fystr = fzstr = nullptr;
  fxvalue = fyvalue = fzvalue = 0.0;

  if (utils::strmatch(arg[0], "^v_")) {
    fxstr = utils::strdup(arg[0] + 2);
  } else
    fxvalue = utils::numeric(FLERR, arg[0], false, lmp);
  if (utils::strmatch(arg[1], "^v_")) {
    fystr = utils::strdup(arg[1] + 2);
  } else
    fyvalue = utils::numeric(FLERR, arg[1], false, lmp);
  if (utils::strmatch(arg[2], "^v_")) {
    fzstr = utils::strdup(arg[2] + 2);
  } else
    fzvalue = utils::numeric(FLERR, arg[2], false, lmp);
}

/* ----------------------------------------------------------------------
   Improved iteration algorithms
------------------------------------------------------------------------- */

double FixIndentForce::get_adaptive_factor(double error, double target, int iteration)
{
  double relative_error = fabs(error / (fabs(target) + 1e-10));

  // More aggressive adaptive step sizes with iteration-based scaling
  double base_factor;
  if (relative_error > 2.0) {
    base_factor = 0.008;    // Much larger step for very large errors
  } else if (relative_error > 1.0) {
    base_factor = 0.005;    // Larger step for large errors
  } else if (relative_error > 0.5) {
    base_factor = 0.003;    // Medium-large error, medium step
  } else if (relative_error > 0.1) {
    base_factor = 0.002;    // Medium error, reasonable step
  } else if (relative_error > 0.01) {
    base_factor = 0.001;    // Small error, small step
  } else {
    base_factor = 0.0005;   // Very small error, minimal step
  }

  // Apply aggressive multiplier for slow convergence
  double aggressive_multiplier = get_aggressive_multiplier(iteration);
  return base_factor * aggressive_multiplier;
}

bool FixIndentForce::detect_oscillation(double current_adj, double prev_adj)
{
  // Detect if consecutive adjustments are in opposite directions
  return (current_adj * prev_adj < 0.0) && (fabs(prev_adj) > 1e-10);
}

double FixIndentForce::estimate_gradient(double current_force, double prev_force,
                                         double current_pos, double prev_pos)
{
  double pos_diff = current_pos - prev_pos;
  if (fabs(pos_diff) < 1e-10) return 0.0;

  double force_diff = current_force - prev_force;
  return force_diff / pos_diff;
}

double FixIndentForce::get_dynamic_limit(double error, double target, int iteration)
{
  double relative_error = fabs(error / (fabs(target) + 1e-10));

  // More aggressive dynamic adjustment limit with iteration scaling
  double base_limit = 0.1;    // Larger base limit
  double dynamic_limit = base_limit * (1.0 + relative_error);    // More aggressive scaling

  // Apply aggressive multiplier for slow convergence
  double aggressive_multiplier = get_aggressive_multiplier(iteration);
  dynamic_limit *= aggressive_multiplier;

  // Higher cap for more aggressive control
  return fmin(dynamic_limit, 0.5);    // Much higher maximum limit
}

double FixIndentForce::get_aggressive_multiplier(int iteration)
{
  // Progressive scaling based on iteration count for slow convergence scenarios
  if (iteration < 20) {
    return 1.0;    // Normal behavior for first 20 iterations
  } else if (iteration < 50) {
    return 1.0 + (iteration - 20) * 0.05;    // Gradually increase: 1.0 to 2.5
  } else if (iteration < 100) {
    return 2.5 + (iteration - 50) * 0.03;    // More aggressive: 2.5 to 4.0
  } else if (iteration < 200) {
    return 4.0 + (iteration - 100) * 0.02;   // Very aggressive: 4.0 to 6.0
  } else {
    return fmin(8.0, 6.0 + (iteration - 200) * 0.01);  // Cap at 8.0x
  }
}

/* ----------------------------------------------------------------------
   Placeholder implementations for geometry functions (add more as needed)
------------------------------------------------------------------------- */

bool FixIndentForce::PointInsideCone(int /*dir*/, double * /*center*/, double /*lo*/, double /*hi*/,
                                     double /*rlo*/, double /*rhi*/, double * /*x*/)
{
  error->all(FLERR, "Cone geometry not implemented in fix indent/force");
  return false;
}

void FixIndentForce::DistanceExteriorPoint(int /*dir*/, double * /*center*/, double /*lo*/,
                                           double /*hi*/, double /*rlo*/, double /*rhi*/,
                                           double & /*x*/, double & /*y*/, double & /*z*/)
{
  error->all(FLERR, "Cone geometry not implemented in fix indent/force");
}

void FixIndentForce::DistanceInteriorPoint(int /*dir*/, double * /*center*/, double /*lo*/,

                                           double /*hi*/, double /*rlo*/, double /*rhi*/,

                                           double & /*x*/, double & /*y*/, double & /*z*/)
{

  error->all(FLERR, "Cone geometry not implemented in fix indent/force");
}

void FixIndentForce::point_on_line_segment(double * /*a*/, double * /*b*/, double * /*c*/,
                                           double * /*d*/)
{
  error->all(FLERR, "Cone geometry not implemented in fix indent/force");
}

double FixIndentForce::closest(double * /*x*/, double * /*near*/, double * /*nearest*/, double dsq)
{
  error->all(FLERR, "Cone geometry not implemented in fix indent/force");
  return dsq;
}