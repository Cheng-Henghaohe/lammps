/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(indent/load/improved,FixIndentLoadImproved);
// clang-format on
#else

#ifndef LMP_FIX_INDENT_LOAD_IMPROVED_H
#define LMP_FIX_INDENT_LOAD_IMPROVED_H

#include "fix.h"
#include <vector>

namespace LAMMPS_NS {

class FixIndentLoadImproved : public Fix {
 public:
  FixIndentLoadImproved(class LAMMPS *, int, char **);
  ~FixIndentLoadImproved() override;
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
  int istyle;                    // indenter style (cylinder, sphere, etc.)
  double k;                      // force constant
  double target_load;            // target load for indenter (nN)
  double current_load;           // current measured load (nN)
  double max_speed;              // maximum speed for correction (m/s)
  double max_speed_overshoot;    // reduced speed during overshoot (m/s)
  double atom_mass_avg;          // average mass of atoms (for unit conversion)
  double sliding_speed;          // sliding speed in X direction (m/s)
  int sliding_dir;               // current sliding direction (+1 or -1)
  int sliding_active;            // 1 if sliding is active, 0 otherwise
  double cyl_radius;             // cylinder radius
  double cyl_x, cyl_y, cyl_z;    // cylinder position
  double indenter[4];            // indenter force and energy storage
  int indenter_flag;             // flag for whether indenter values are current
  int ilevel_respa;              // rRESPA level to apply force
  double dt;                     // timestep

  // PID控制相关变量
  double pid_kp, pid_ki, pid_kd;          // PID参数
  double integral_error;                   // 积分误差
  double previous_error;                   // 上一步误差
  double max_integral_windup;             // 最大积分饱和值
  std::vector<double> error_history;      // 误差历史
  int pid_history_length;                  // 历史长度
  double load_tolerance;                   // 载荷容差百分比
  
  // 自适应PID参数
  double adaptive_factor;                  // 自适应调整因子
  double min_kp, max_kp;                  // Kp的范围
  bool enable_adaptive_pid;               // 是否启用自适应PID
  
  // 统计和诊断
  int control_step_count;                 // 控制步数
  double max_overshoot;                   // 最大超调量
  double settling_time;                   // 稳定时间
  bool load_stable;                       // 载荷是否稳定
  
  // 性能指标结构
  struct PerformanceMetrics {
    double avg_error;
    double max_error;
    double steady_state_error;
    int oscillation_count;
    bool is_stable;
  };

  int get_indenter_type(const char *);
  void control_load_pid();                // PID控制函数
  void update_indenter_position();
  void update_sliding_position();
  void init_pid_parameters();             // 初始化PID参数
  void tune_pid_parameters();             // 自动调节PID参数
  bool is_load_stable(double tolerance_percent = 1.0);
  void update_statistics();               // 更新统计信息
  void convert_units_properly();          // 改进的单位转换
  
  // 自动PID计算功能
  void auto_initialize_pid();             // 自动初始化PID参数
  std::tuple<double, double, double> calculate_physical_pid();    // 物理计算法
  std::tuple<double, double, double> calculate_empirical_pid();   // 经验公式法
  void adaptive_pid_adjustment();         // 自适应调节
  PerformanceMetrics calculate_performance_metrics();  // 性能评估
};

}    // namespace LAMMPS_NS

#endif
#endif