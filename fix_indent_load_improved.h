/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/
   
   FixIndentLoadImproved - 纳米压痕载荷控制系统
   作者：基于LAMMPS框架开发
   功能：实现自适应PID控制的纳米压痕载荷控制
   特点：
   - 自动计算最优PID参数，无需手动调参
   - 支持圆柱形和球形压头
   - 实时载荷监控和自适应控制
   - 防止超调和振荡的保护机制
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(indent/load/improved,FixIndentLoadImproved);
// clang-format on
#else

#ifndef LMP_FIX_INDENT_LOAD_IMPROVED_H
#define LMP_FIX_INDENT_LOAD_IMPROVED_H

#include "fix.h"       // LAMMPS基础Fix类
#include <vector>      // STL向量容器，用于存储历史数据
#include <tuple>       // STL元组，用于返回多个值

namespace LAMMPS_NS {

/**
 * FixIndentLoadImproved类：实现纳米压痕仿真中的自动载荷控制
 * 
 * 继承关系：Fix (LAMMPS基础Fix类)
 * 
 * 主要功能：
 * 1. 使用PID控制算法自动调节压头位置以达到目标载荷
 * 2. 支持圆柱形和球形压头几何形状
 * 3. 自动计算最优PID参数，避免手动调参
 * 4. 实时监控载荷稳定性和控制性能
 * 5. 提供自适应参数调节功能
 * 
 * 使用方法：
 * fix ID group indent/load/improved k target_load max_speed overshoot_speed type geometry_params [options]
 * 
 * 例如：
 * fix 1 atoms indent/load/improved 1600 500 0.1 0.05 cylinder 0 0 50 8
 */
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
  // ======================== 压头几何参数 ========================
  int istyle;                    // 压头类型：CYLINDER(圆柱形)或SPHERE(球形)
  double cyl_radius;             // 压头半径 (单位：Å)
  double cyl_x, cyl_y, cyl_z;    // 压头位置坐标 (单位：Å)
  
  // ======================== 物理控制参数 ========================
  double k;                      // 接触力常数 (单位：GPa，控制压头与原子相互作用强度)
  double target_load;            // 目标载荷值 (单位：nN，我们要达到的载荷)
  double current_load;           // 当前实际载荷 (单位：nN，由力计算得出)
  double max_speed;              // 最大调整速度 (单位：m/s，限制压头移动速度)
  double max_speed_overshoot;    // 超调时的最大速度 (单位：m/s，防止过度调节)
  
  // ======================== 系统参数 ========================
  double dt;                     // 时间步长 (单位：ps，从LAMMPS获取)
  double atom_mass_avg;          // 原子平均质量 (用于单位转换和物理估算)
  int ilevel_respa;              // RESPA积分层级 (多尺度时间积分用)
  
  // ======================== 滑动功能参数 ========================
  double sliding_speed;          // 滑动速度 (单位：m/s，压头在x方向的滑动)
  int sliding_dir;               // 滑动方向 (+1正方向，-1负方向)
  int sliding_active;            // 滑动开关 (1开启，0关闭)
  
  // ======================== 输出数据存储 ========================
  double indenter[4];            // 压头数据：[势能, fx, fy, fz]
  int indenter_flag;             // 压头数据有效性标志

  // ======================== PID控制核心参数 ========================
  double pid_kp, pid_ki, pid_kd;          // PID控制器增益参数
                                          // pid_kp: 比例增益，响应当前误差
                                          // pid_ki: 积分增益，消除稳态误差
                                          // pid_kd: 微分增益，预测未来趋势，提高稳定性
  
  double integral_error;                   // 误差积分累积值 (消除稳态误差用)
  double previous_error;                   // 上一时间步的误差值 (计算微分用)
  double max_integral_windup;             // 积分饱和限制值 (防止积分项过大)
  
  // ======================== 历史数据管理 ========================
  std::vector<double> error_history;      // 误差历史记录 (用于性能分析和自适应调节)
  int pid_history_length;                  // 保留的历史数据长度
  double load_tolerance;                   // 载荷稳定判断容差 (百分比)
  
  // ======================== 自适应控制参数 ========================
  double adaptive_factor;                  // 自适应调节强度因子
  double min_kp, max_kp;                  // 比例增益的允许范围 (防止参数偏离过远)
  bool enable_adaptive_pid;               // 自适应功能开关
  
  // ======================== 性能监控与统计 ========================
  int control_step_count;                 // 控制步数计数器
  double max_overshoot;                   // 记录的最大超调量 (性能评估用)
  double settling_time;                   // 系统稳定时间
  bool load_stable;                       // 当前载荷稳定状态标志
  
  /**
   * 性能指标结构体：用于评估PID控制系统的性能
   * 
   * 这个结构体包含了评估控制系统质量的关键指标：
   * - avg_error: 平均误差百分比，反映整体控制精度
   * - max_error: 最大误差百分比，反映最坏情况
   * - steady_state_error: 稳态误差，反映长期控制精度
   * - oscillation_count: 振荡次数，反映系统稳定性
   * - is_stable: 整体稳定性判断
   */
  struct PerformanceMetrics {
    double avg_error;          // 平均误差百分比
    double max_error;          // 最大误差百分比
    double steady_state_error; // 稳态误差百分比
    int oscillation_count;     // 振荡次数计数
    bool is_stable;           // 系统稳定性标志
  };

  // ======================== 私有成员函数 ========================
  
  // 基础功能函数
  int get_indenter_type(const char *);   // 解析压头类型字符串
  void update_sliding_position();        // 更新滑动位置
  void convert_units_properly();          // 单位转换处理
  void update_statistics();               // 更新统计信息和输出诊断
  
  // PID控制核心函数
  void control_load_pid();                // PID控制主函数（每时间步调用）
  void init_pid_parameters();             // 初始化PID控制相关参数
  bool is_load_stable(double tolerance_percent = 1.0);  // 载荷稳定性检查
  
  // 自动PID参数计算系统
  void auto_initialize_pid();             // 自动计算最优PID参数
  std::tuple<double, double, double> calculate_physical_pid();    // 基于物理模型计算PID
  std::tuple<double, double, double> calculate_empirical_pid();   // 基于经验公式计算PID
  
  // 自适应控制系统
  void adaptive_pid_adjustment();         // 实时自适应调节PID参数
  PerformanceMetrics calculate_performance_metrics();  // 计算控制性能指标
};

}    // namespace LAMMPS_NS

#endif
#endif