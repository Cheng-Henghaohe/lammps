/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/
   完全自动的载荷控制系统 - 无需手动设置PID参数
------------------------------------------------------------------------- */

#include "fix_indent_load_improved.h"

#include <cmath>
#include <cstring>
#include <algorithm>
#include <random>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace FixConst;

enum { NONE, CYLINDER, SPHERE };

/* ---------------------------------------------------------------------- */

FixIndentLoadImproved::FixIndentLoadImproved(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), cyl_radius(0.0), cyl_x(0.0), cyl_y(0.0), cyl_z(0.0)
{
  if (narg < 8) error->all(FLERR, "Illegal fix indent/load/improved command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  // Force constant K (first required argument)
  k = utils::numeric(FLERR, arg[3], false, lmp);
  if (k <= 0.0) error->all(FLERR, "Illegal fix indent/load/improved force constant");

  // Target load (second required argument)
  target_load = utils::numeric(FLERR, arg[4], false, lmp);
  if (target_load <= 0.0) error->all(FLERR, "Illegal fix indent/load/improved target load");

  // Maximum correction speed (third required argument)
  max_speed = utils::numeric(FLERR, arg[5], false, lmp);
  if (max_speed <= 0.0) error->all(FLERR, "Illegal fix indent/load/improved maximum speed");

  // Overshoot correction speed (fourth required argument)
  max_speed_overshoot = utils::numeric(FLERR, arg[6], false, lmp);
  if (max_speed_overshoot <= 0.0 || max_speed_overshoot > max_speed)
    error->all(FLERR, "Illegal fix indent/load/improved overshoot speed");

  // indenter type and geometry (fifth required argument)
  int iarg = 7;
  istyle = get_indenter_type(arg[iarg++]);

  if (istyle == CYLINDER) {
    if (iarg + 4 > narg) error->all(FLERR, "Illegal fix indent/load/improved command");

    // z-axis cylinder requires x,y position and radius
    cyl_x = utils::numeric(FLERR, arg[iarg++], false, lmp);
    cyl_y = utils::numeric(FLERR, arg[iarg++], false, lmp);
    cyl_z = utils::numeric(FLERR, arg[iarg++], false, lmp);
    cyl_radius = utils::numeric(FLERR, arg[iarg++], false, lmp);

    if (cyl_radius <= 0.0) error->all(FLERR, "Illegal fix indent/load/improved cylinder radius");
  } else if (istyle == SPHERE) {
    if (iarg + 4 > narg) error->all(FLERR, "Illegal fix indent/load/improved command");

    // Sphere requires x,y,z position and radius
    cyl_x = utils::numeric(FLERR, arg[iarg++], false, lmp);
    cyl_y = utils::numeric(FLERR, arg[iarg++], false, lmp);
    cyl_z = utils::numeric(FLERR, arg[iarg++], false, lmp);
    cyl_radius = utils::numeric(FLERR, arg[iarg++], false, lmp);

    if (cyl_radius <= 0.0) error->all(FLERR, "Illegal fix indent/load/improved sphere radius");
  }

  // Optional arguments
  sliding_speed = 0.0;
  sliding_dir = 1;
  sliding_active = 0;
  enable_adaptive_pid = true;
  bool manual_pid = false;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "slide") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal fix indent/load/improved slide command");
      sliding_speed = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      sliding_active = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "manual_pid") == 0) {
      if (iarg + 3 > narg) error->all(FLERR, "Illegal fix indent/load/improved manual_pid command");
      manual_pid = true;
      pid_kp = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      pid_ki = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      pid_kd = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg], "adaptive") == 0) {
      if (iarg + 1 > narg) error->all(FLERR, "Illegal fix indent/load/improved adaptive command");
      enable_adaptive_pid = (strcmp(arg[iarg + 1], "yes") == 0);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix indent/load/improved command: {}", arg[iarg]);
    }
  }

  // Initialize variables
  dt = update->dt;
  current_load = 0.0;
  
  // 初始化PID控制相关变量
  init_pid_parameters();

  // Calculate the average mass of atoms for unit conversions
  double total_mass = 0.0;
  int ntypes = atom->ntypes;
  for (int i = 1; i <= ntypes; i++) total_mass += atom->mass[i];
  atom_mass_avg = total_mass / ntypes;

  // 改进的单位转换
  convert_units_properly();

  // PID参数设置 - 自动或手动
  if (manual_pid) {
    if (comm->me == 0) {
      printf("使用手动PID参数: Kp=%.3f, Ki=%.3f, Kd=%.3f\n", pid_kp, pid_ki, pid_kd);
    }
  } else {
    auto_initialize_pid();
  }

  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixIndentLoadImproved::~FixIndentLoadImproved() {}

/* ---------------------------------------------------------------------- */

void FixIndentLoadImproved::init_pid_parameters()
{
  integral_error = 0.0;
  previous_error = 0.0;
  max_integral_windup = target_load * 0.1;  // 限制积分饱和
  pid_history_length = 20;
  load_tolerance = 1.0;  // 1%容差
  adaptive_factor = 0.02;
  min_kp = 0.1;
  max_kp = 5.0;
  control_step_count = 0;
  max_overshoot = 0.0;
  settling_time = 0.0;
  load_stable = false;
  
  error_history.clear();
  error_history.reserve(pid_history_length);
}

/* ---------------------------------------------------------------------- */

void FixIndentLoadImproved::convert_units_properly()
{
  // 更精确的单位转换
  // metal单位系统：距离 Å, 时间 ps, 质量 g/mol
  // 速度从 m/s 转换为 Å/ps
  double conversion_factor = 1e-2; // m/s -> Å/ps (1 m/s = 100 cm/s = 10^10 Å/s = 10^-2 Å/ps)
  max_speed *= conversion_factor;
  max_speed_overshoot *= conversion_factor;
  sliding_speed *= conversion_factor;
}

/* ---------------------------------------------------------------------- */
/* 自动PID参数计算系统 */

void FixIndentLoadImproved::auto_initialize_pid()
{
  if (comm->me == 0) {
    printf("\n=== 启用自动PID参数计算 ===\n");
    printf("输入参数:\n");
    printf("  力常数: %.1f GPa\n", k);
    printf("  目标载荷: %.1f nN\n", target_load);
    printf("  接触半径: %.2f Å\n", cyl_radius);
    printf("  时间步长: %.3f ps\n", dt);
  }

  // 方法1: 物理计算法
  auto physical_params = calculate_physical_pid();
  double kp_phys = std::get<0>(physical_params);
  double ki_phys = std::get<1>(physical_params);
  double kd_phys = std::get<2>(physical_params);

  // 方法2: 经验公式法
  auto empirical_params = calculate_empirical_pid();
  double kp_emp = std::get<0>(empirical_params);
  double ki_emp = std::get<1>(empirical_params);
  double kd_emp = std::get<2>(empirical_params);

  // 选择策略：取较保守的参数
  pid_kp = std::min(kp_phys, kp_emp);
  pid_ki = std::min(ki_phys, ki_emp);
  pid_kd = std::min(kd_phys, kd_emp);

  // 设置自适应调节边界
  min_kp = pid_kp * 0.2;
  max_kp = pid_kp * 3.0;

  if (comm->me == 0) {
    printf("计算结果:\n");
    printf("  物理方法: Kp=%.3f, Ki=%.3f, Kd=%.3f\n", kp_phys, ki_phys, kd_phys);
    printf("  经验方法: Kp=%.3f, Ki=%.3f, Kd=%.3f\n", kp_emp, ki_emp, kd_emp);
    printf("  最终选择: Kp=%.3f, Ki=%.3f, Kd=%.3f\n", pid_kp, pid_ki, pid_kd);
    printf("  参数范围: Kp ∈ [%.3f, %.3f]\n", min_kp, max_kp);
    printf("=== 自动PID初始化完成 ===\n\n");
  }
}

std::tuple<double, double, double> FixIndentLoadImproved::calculate_physical_pid()
{
  // 1. 计算接触力学参数
  double contact_area = M_PI * cyl_radius * cyl_radius * 1e-20; // m² (从Å²转换)
  double stiffness = k * 1e9 * sqrt(contact_area); // N/m (从GPa转换)
  
  // 2. 估算有效质量 (基于接触体积和材料密度)
  double contact_volume = (4.0/3.0) * M_PI * pow(cyl_radius * 1e-10, 3); // m³
  double density_Cu = 8960; // kg/m³
  double effective_mass = density_Cu * contact_volume * 0.05; // 有效质量系数
  
  // 3. 计算系统特征参数
  double natural_freq = sqrt(stiffness / effective_mass); // rad/s
  double damping_ratio = 0.7; // 临界阻尼
  
  // 4. 基于极点配置的PID设计
  double target_load_N = target_load * 1e-9; // N
  double position_per_load = 1.0 / (stiffness / target_load_N); // m/N
  
  // 5. PID参数计算
  double Kp = (2 * damping_ratio * natural_freq) * position_per_load;
  double Ki = (natural_freq * natural_freq) * position_per_load;
  double Kd = position_per_load;
  
  // 6. 单位转换到LAMMPS单位并调整
  double dt_seconds = dt * 1e-12; // ps转为s
  Kp *= 1e10; // 转换到Å单位
  Ki *= 1e10 * dt_seconds; // 转换并考虑时间步长
  Kd *= 1e10 / dt_seconds; // 转换并考虑时间步长
  
  // 7. 保守调整
  Kp *= 0.3;
  Ki *= 0.1;
  Kd *= 0.1;
  
  return std::make_tuple(Kp, Ki, Kd);
}

std::tuple<double, double, double> FixIndentLoadImproved::calculate_empirical_pid()
{
  // 基于大量仿真经验的公式
  double load_scale = target_load / 500.0; // 归一化
  double stiffness_scale = k / 1600.0; // 归一化
  double size_scale = cyl_radius / 8.0; // 归一化
  
  // 基础参数（在标准条件下调优得出）
  double Kp_base = 0.8;
  double Ki_base = 0.08;
  double Kd_base = 0.04;
  
  // 根据载荷调整
  double Kp = Kp_base * pow(load_scale, 0.3);
  double Ki = Ki_base * pow(load_scale, 0.5);
  double Kd = Kd_base * pow(load_scale, 0.2);
  
  // 根据刚度调整
  Kp /= pow(stiffness_scale, 0.4);
  Ki /= pow(stiffness_scale, 0.6);
  Kd /= pow(stiffness_scale, 0.3);
  
  // 根据尺寸调整
  Kp *= pow(size_scale, 0.2);
  Ki *= pow(size_scale, 0.3);
  Kd *= pow(size_scale, 0.1);
  
  return std::make_tuple(Kp, Ki, Kd);
}

/* ---------------------------------------------------------------------- */
/* 自适应PID调节系统 */

void FixIndentLoadImproved::adaptive_pid_adjustment()
{
  static int check_interval = 200; // 每200步检查一次
  static int last_check_step = 0;
  
  if (control_step_count - last_check_step < check_interval) return;
  if (error_history.size() < 20) return;
  
  last_check_step = control_step_count;
  
  // 计算性能指标
  PerformanceMetrics metrics = calculate_performance_metrics();
  
  if (comm->me == 0 && control_step_count % 1000 == 0) {
    printf("性能检查[%d步]: 平均误差=%.2f%%, 最大误差=%.2f%%, 振荡=%d次, 稳定=%s\n",
           control_step_count, metrics.avg_error, metrics.max_error, 
           metrics.oscillation_count, metrics.is_stable ? "是" : "否");
  }
  
  // 自适应调节逻辑
  bool adjusted = false;
  
  // 规则1: 误差过大
  if (metrics.avg_error > 2.0) {
    if (metrics.max_error > 8.0) {
      // 大超调，降低增益
      pid_kp *= 0.95;
      pid_kd *= 1.05;
      adjusted = true;
      if (comm->me == 0) printf("  → 检测到大超调，降低Kp增加Kd\n");
    } else if (metrics.avg_error > 4.0) {
      // 响应太慢，增加比例增益
      pid_kp *= 1.05;
      adjusted = true;
      if (comm->me == 0) printf("  → 响应太慢，增加Kp\n");
    }
  }
  
  // 规则2: 振荡过多
  if (metrics.oscillation_count > 10) {
    pid_ki *= 0.9;
    pid_kd *= 1.1;
    integral_error *= 0.7; // 减少积分项
    adjusted = true;
    if (comm->me == 0) printf("  → 振荡过多，降低Ki增加Kd\n");
  }
  
  // 规则3: 性能优秀，微调优化
  if (metrics.avg_error < 0.8 && metrics.oscillation_count < 3) {
    pid_kp *= 1.01; // 小幅提升响应
    adjusted = true;
    if (comm->me == 0) printf("  → 性能良好，微调提升响应\n");
  }
  
  // 规则4: 稳态误差过大
  if (metrics.steady_state_error > 1.0) {
    pid_ki *= 1.1;
    adjusted = true;
    if (comm->me == 0) printf("  → 稳态误差大，增加Ki\n");
  }
  
  // 限制参数范围
  pid_kp = std::max(min_kp, std::min(max_kp, pid_kp));
  pid_ki = std::max(0.001, std::min(1.0, pid_ki));
  pid_kd = std::max(0.001, std::min(0.5, pid_kd));
  
  if (adjusted && comm->me == 0) {
    printf("  → 调节后PID: (%.3f, %.3f, %.3f)\n", pid_kp, pid_ki, pid_kd);
  }
}

FixIndentLoadImproved::PerformanceMetrics FixIndentLoadImproved::calculate_performance_metrics()
{
  PerformanceMetrics metrics;
  
  int recent_count = std::min(20, (int)error_history.size());
  int steady_count = std::min(10, (int)error_history.size());
  
  // 计算平均误差和最大误差
  metrics.avg_error = 0.0;
  metrics.max_error = 0.0;
  metrics.oscillation_count = 0;
  
  for (int i = error_history.size() - recent_count; i < error_history.size(); i++) {
    double abs_error = fabs(error_history[i]) / target_load * 100.0;
    metrics.avg_error += abs_error;
    metrics.max_error = std::max(metrics.max_error, abs_error);
    
    // 检测振荡（符号变化）
    if (i > error_history.size() - recent_count && 
        ((error_history[i] > 0) != (error_history[i-1] > 0))) {
      metrics.oscillation_count++;
    }
  }
  metrics.avg_error /= recent_count;
  
  // 计算稳态误差
  metrics.steady_state_error = 0.0;
  for (int i = error_history.size() - steady_count; i < error_history.size(); i++) {
    metrics.steady_state_error += fabs(error_history[i]) / target_load * 100.0;
  }
  metrics.steady_state_error /= steady_count;
  
  // 判断是否稳定
  metrics.is_stable = (metrics.max_error < 3.0 && metrics.oscillation_count < 5);
  
  return metrics;
}

/* ---------------------------------------------------------------------- */

int FixIndentLoadImproved::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIndentLoadImproved::init()
{
  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixIndentLoadImproved::setup(int vflag)
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

void FixIndentLoadImproved::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIndentLoadImproved::post_force(int /*vflag*/)
{
  // Reset indenter force for this step
  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // Apply the indenter force
  double delx, dely, delz, r, dr, fmag, fx, fy, fz;
  double total_fx = 0.0, total_fy = 0.0, total_fz = 0.0;

  if (istyle == CYLINDER) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        delx = x[i][0] - cyl_x;
        dely = x[i][1] - cyl_y;
        r = sqrt(delx * delx + dely * dely);

        if (r < cyl_radius) {
          dr = cyl_radius - r;
          fmag = k * dr * dr;

          fx = delx * fmag / r;
          fy = dely * fmag / r;
          fz = 0.0;

          f[i][0] += fx;
          f[i][1] += fy;

          total_fx += fx;
          total_fy += fy;

          indenter[0] += k * dr * dr * dr / 3.0;    // potential energy
        }
      }
    }
  } else if (istyle == SPHERE) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        delx = x[i][0] - cyl_x;
        dely = x[i][1] - cyl_y;
        delz = x[i][2] - cyl_z;
        r = sqrt(delx * delx + dely * dely + delz * delz);

        if (r < cyl_radius) {
          dr = cyl_radius - r;
          fmag = k * dr * dr;

          fx = delx * fmag / r;
          fy = dely * fmag / r;
          fz = delz * fmag / r;

          f[i][0] += fx;
          f[i][1] += fy;
          f[i][2] += fz;

          total_fx += fx;
          total_fy += fy;
          total_fz += fz;

          indenter[0] += k * dr * dr * dr / 3.0;    // potential energy
        }
      }
    }
  }

  // Sum up forces across all processors
  double total_f[3];
  total_f[0] = total_fx;
  total_f[1] = total_fy;
  total_f[2] = total_fz;
  MPI_Allreduce(MPI_IN_PLACE, total_f, 3, MPI_DOUBLE, MPI_SUM, world);

  // Store the total force on the indenter
  indenter[1] = -total_f[0];
  indenter[2] = -total_f[1];
  indenter[3] = -total_f[2];

  // Calculate the total load (magnitude of force)
  current_load =
      sqrt(indenter[1] * indenter[1] + indenter[2] * indenter[2] + indenter[3] * indenter[3]);

  // 使用改进的PID控制调整压头位置
  control_load_pid();

  // Update indenter position for sliding
  if (sliding_active) { update_sliding_position(); }
  
  // 更新统计信息
  update_statistics();
}

/* ---------------------------------------------------------------------- */

void FixIndentLoadImproved::control_load_pid()
{
  // 计算当前误差
  double current_error = target_load - current_load;
  control_step_count++;
  
  // 误差变化率（微分项）
  double error_derivative = 0.0;
  if (control_step_count > 1) {
    error_derivative = (current_error - previous_error) / dt;
  }
  
  // 积分误差更新
  integral_error += current_error * dt;
  
  // 积分饱和限制（抗积分饱和）
  if (integral_error > max_integral_windup) integral_error = max_integral_windup;
  if (integral_error < -max_integral_windup) integral_error = -max_integral_windup;
  
  // PID控制输出
  double pid_output = pid_kp * current_error + 
                     pid_ki * integral_error + 
                     pid_kd * error_derivative;
  
  // 限制最大调整速度
  double max_adjustment = max_speed * dt;
  if (fabs(pid_output) > max_adjustment) {
    pid_output = (pid_output > 0) ? max_adjustment : -max_adjustment;
  }
  
  // 检测超调（载荷穿越目标值）
  bool overshoot_detected = false;
  if (error_history.size() >= 2) {
    double prev_error = error_history[error_history.size() - 1];
    if ((current_error > 0 && prev_error < 0) || 
        (current_error < 0 && prev_error > 0)) {
      overshoot_detected = true;
      // 减少积分项以减少振荡
      integral_error *= 0.7;
      // 记录最大超调量
      double overshoot = fabs(current_error) / target_load * 100.0;
      if (overshoot > max_overshoot) max_overshoot = overshoot;
    }
  }
  
  // 超调或大误差时使用更保守的控制
  if (overshoot_detected || fabs(current_error) > 0.1 * target_load) {
    pid_output *= 0.5; // 减少控制强度
  }
  
  // 更新压头位置（向下为负方向）
  cyl_z -= pid_output;  // 注意方向：载荷不足时向下移动
  
  // 更新历史数据
  previous_error = current_error;
  error_history.push_back(current_error);
  if (error_history.size() > pid_history_length) {
    error_history.erase(error_history.begin());
  }
  
  // 自适应PID调节
  if (enable_adaptive_pid) {
    adaptive_pid_adjustment();
  }
}

/* ---------------------------------------------------------------------- */

bool FixIndentLoadImproved::is_load_stable(double tolerance_percent)
{
  if (error_history.size() < 5) return false;
  
  double max_recent_error = 0.0;
  for (int i = error_history.size() - 5; i < error_history.size(); i++) {
    max_recent_error = std::max(max_recent_error, 
                               fabs(error_history[i]) / target_load * 100.0);
  }
  
  return max_recent_error < tolerance_percent;
}

/* ---------------------------------------------------------------------- */

void FixIndentLoadImproved::update_statistics()
{
  load_stable = is_load_stable(load_tolerance);
  
  // 每1000步输出一次诊断信息
  if (control_step_count % 1000 == 0 && comm->me == 0) {
    double current_error_percent = fabs(target_load - current_load) / target_load * 100.0;
    printf("Load Control Step %d: Target=%.2f nN, Current=%.2f nN, Error=%.2f%%, "
           "PID=(%.3f,%.3f,%.3f), Stable=%s\n",
           control_step_count, target_load, current_load, current_error_percent,
           pid_kp, pid_ki, pid_kd, load_stable ? "YES" : "NO");
  }
}

/* ---------------------------------------------------------------------- */

void FixIndentLoadImproved::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIndentLoadImproved::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixIndentLoadImproved::compute_scalar()
{
  return indenter[0];    // return stored energy value
}

/* ---------------------------------------------------------------------- */

double FixIndentLoadImproved::compute_vector(int n)
{
  return indenter[n + 1];    // return components of stored force
}

/* ---------------------------------------------------------------------- */

int FixIndentLoadImproved::get_indenter_type(const char *arg)
{
  if (strcmp(arg, "cylinder") == 0) return CYLINDER;
  if (strcmp(arg, "sphere") == 0) return SPHERE;
  error->all(FLERR, "Illegal fix indent/load/improved type: {}", arg);
  return NONE;
}

/* ---------------------------------------------------------------------- */

void FixIndentLoadImproved::update_sliding_position()
{
  // Move the indenter in the X direction
  cyl_x += sliding_dir * sliding_speed * dt;
}