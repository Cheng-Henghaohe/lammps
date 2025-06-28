/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/
   完全自动的载荷控制系统 - 无需手动设置PID参数
   
   功能说明：
   - 实现纳米压痕仿真中的精确载荷控制
   - 使用PID控制算法自动调节压头位置
   - 支持自适应参数调节，无需手动调参
   - 适用于圆柱形和球形压头
------------------------------------------------------------------------- */

// 包含自定义头文件
#include "fix_indent_load_improved.h"

// C++标准库头文件
#include <cmath>      // 数学函数：sqrt, pow, fabs等
#include <cstring>    // 字符串操作：strcmp等
#include <algorithm>  // STL算法：min, max等
#include <random>     // 随机数生成（未使用但保留）

// LAMMPS核心头文件
#include "atom.h"      // 原子数据结构
#include "comm.h"      // MPI通信接口
#include "domain.h"    // 仿真域定义
#include "error.h"     // 错误处理
#include "force.h"     // 力计算接口
#include "input.h"     // 输入解析
#include "math_const.h" // 数学常数（如π）
#include "respa.h"     // RESPA积分算法
#include "update.h"    // 时间步更新
#include "variable.h"  // 变量管理
#include "utils.h"     // LAMMPS实用函数（包含utils::numeric）

// 使用LAMMPS命名空间，避免名称冲突
using namespace LAMMPS_NS;  // LAMMPS主命名空间
using namespace MathConst;   // 数学常数命名空间
using namespace FixConst;    // Fix相关常数命名空间

// 压头类型枚举：定义支持的压头几何形状
enum { NONE, CYLINDER, SPHERE };

/* ---------------------------------------------------------------------- */
/* 构造函数：创建载荷控制对象并初始化所有参数
   参数说明：
   - lmp: LAMMPS主对象指针
   - narg: 命令行参数个数
   - arg: 命令行参数数组
   
   C++知识点：
   - 构造函数初始化列表：在函数体执行前初始化成员变量
   - 继承：调用父类Fix的构造函数
*/

FixIndentLoadImproved::FixIndentLoadImproved(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),  // 调用父类Fix的构造函数
    cyl_radius(0.0), cyl_x(0.0), cyl_y(0.0), cyl_z(0.0)  // 初始化压头几何参数
{
  // 参数验证：确保至少有8个必需参数
  // 命令格式：fix ID group_name indent/load/improved k target_load max_speed overshoot_speed type ...
  if (narg < 8) error->all(FLERR, "Illegal fix indent/load/improved command");

  // 设置LAMMPS输出标志：告诉LAMMPS这个fix能提供什么类型的输出
  scalar_flag = 1;           // 可以输出标量值（能量）
  vector_flag = 1;           // 可以输出矢量值（力的分量）
  size_vector = 3;           // 矢量大小为3（fx, fy, fz）
  global_freq = 1;           // 每个时间步都更新输出
  extscalar = 1;             // 扩展标量输出
  extvector = 1;             // 扩展矢量输出
  respa_level_support = 1;   // 支持RESPA多尺度积分
  ilevel_respa = 0;          // RESPA积分层级

  // 解析必需参数 - arg[0]=fix, arg[1]=ID, arg[2]=group, arg[3]开始是我们的参数
  
  // 参数1：力常数K (单位：GPa) - 决定压头与材料接触时的刚度
  // C++知识点：utils::numeric()将字符串转换为数值
  k = utils::numeric(FLERR, arg[3], false, lmp);
  if (k <= 0.0) error->all(FLERR, "Illegal fix indent/load/improved force constant");

  // 参数2：目标载荷 (单位：nN) - 我们想要达到的载荷值
  target_load = utils::numeric(FLERR, arg[4], false, lmp);
  if (target_load <= 0.0) error->all(FLERR, "Illegal fix indent/load/improved target load");

  // 参数3：最大调整速度 (单位：m/s) - 限制压头移动的最大速度
  max_speed = utils::numeric(FLERR, arg[5], false, lmp);
  if (max_speed <= 0.0) error->all(FLERR, "Illegal fix indent/load/improved maximum speed");

  // 参数4：超调修正速度 (单位：m/s) - 发生超调时的最大修正速度
  max_speed_overshoot = utils::numeric(FLERR, arg[6], false, lmp);
  if (max_speed_overshoot <= 0.0 || max_speed_overshoot > max_speed)
    error->all(FLERR, "Illegal fix indent/load/improved overshoot speed");

  // 参数5：压头类型和几何参数
  int iarg = 7;  // 当前参数索引
  // C++知识点：后缀++运算符，先使用值再自增
  istyle = get_indenter_type(arg[iarg++]);  // 解析压头类型（cylinder或sphere）

  // 根据压头类型解析几何参数
  if (istyle == CYLINDER) {
    // 圆柱形压头：需要x,y,z位置和半径
    if (iarg + 4 > narg) error->all(FLERR, "Illegal fix indent/load/improved command");

    cyl_x = utils::numeric(FLERR, arg[iarg++], false, lmp);      // 圆柱中心x坐标
    cyl_y = utils::numeric(FLERR, arg[iarg++], false, lmp);      // 圆柱中心y坐标
    cyl_z = utils::numeric(FLERR, arg[iarg++], false, lmp);      // 圆柱底面z坐标
    cyl_radius = utils::numeric(FLERR, arg[iarg++], false, lmp); // 圆柱半径

    if (cyl_radius <= 0.0) error->all(FLERR, "Illegal fix indent/load/improved cylinder radius");
  } else if (istyle == SPHERE) {
    // 球形压头：需要x,y,z位置和半径
    if (iarg + 4 > narg) error->all(FLERR, "Illegal fix indent/load/improved command");

    cyl_x = utils::numeric(FLERR, arg[iarg++], false, lmp);      // 球心x坐标
    cyl_y = utils::numeric(FLERR, arg[iarg++], false, lmp);      // 球心y坐标
    cyl_z = utils::numeric(FLERR, arg[iarg++], false, lmp);      // 球心z坐标
    cyl_radius = utils::numeric(FLERR, arg[iarg++], false, lmp); // 球半径

    if (cyl_radius <= 0.0) error->all(FLERR, "Illegal fix indent/load/improved sphere radius");
  }

  // 初始化可选参数的默认值
  sliding_speed = 0.0;           // 滑动速度（默认不滑动）
  sliding_dir = 1;               // 滑动方向（1为正方向）
  sliding_active = 0;            // 滑动功能开关（默认关闭）
  enable_adaptive_pid = true;    // 自适应PID调节（默认开启）
  bool manual_pid = false;       // 手动PID参数标志

  // C++知识点：while循环解析可选参数
  // 使用strcmp()比较C风格字符串
  while (iarg < narg) {
    if (strcmp(arg[iarg], "slide") == 0) {
      // slide选项：启用压头滑动功能
      if (iarg + 1 > narg) error->all(FLERR, "Illegal fix indent/load/improved slide command");
      sliding_speed = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      sliding_active = 1;  // 激活滑动
      iarg += 2;           // 跳过已处理的参数
    } else if (strcmp(arg[iarg], "manual_pid") == 0) {
      // manual_pid选项：手动指定PID参数，不使用自动计算
      if (iarg + 3 > narg) error->all(FLERR, "Illegal fix indent/load/improved manual_pid command");
      manual_pid = true;
      pid_kp = utils::numeric(FLERR, arg[iarg + 1], false, lmp);  // 比例增益
      pid_ki = utils::numeric(FLERR, arg[iarg + 2], false, lmp);  // 积分增益
      pid_kd = utils::numeric(FLERR, arg[iarg + 3], false, lmp);  // 微分增益
      iarg += 4;
    } else if (strcmp(arg[iarg], "adaptive") == 0) {
      // adaptive选项：控制是否启用自适应PID调节
      if (iarg + 1 > narg) error->all(FLERR, "Illegal fix indent/load/improved adaptive command");
      enable_adaptive_pid = (strcmp(arg[iarg + 1], "yes") == 0);
      iarg += 2;
    } else {
      // 遇到未知参数时报错
      error->all(FLERR, "Illegal fix indent/load/improved command: {}", arg[iarg]);
    }
  }

  // 初始化基本变量
  dt = update->dt;               // 获取时间步长（从LAMMPS更新器）
  current_load = 0.0;            // 当前载荷初始值
  
  // 调用PID控制系统初始化函数
  init_pid_parameters();

  // 计算原子平均质量：用于单位转换和物理参数估算
  double total_mass = 0.0;
  int ntypes = atom->ntypes;     // 获取原子类型数量
  
  // 安全检查：防止除零错误
  if (ntypes <= 0) {
    error->all(FLERR, "No atom types defined in fix indent/load/improved");
  }
  
  // C++知识点：for循环累加计算
  for (int i = 1; i <= ntypes; i++) {
    if (atom->mass[i] > 0.0) {  // 只累加有效质量
      total_mass += atom->mass[i];  // 累加各类型质量
    }
  }
  
  // 安全检查：防止除零错误
  if (total_mass <= 0.0) {
    atom_mass_avg = 63.546;  // 使用铜的原子质量作为默认值
    if (comm->me == 0) {
      printf("警告：未找到有效原子质量，使用默认值 %.3f g/mol\n", atom_mass_avg);
    }
  } else {
    atom_mass_avg = total_mass / ntypes;  // 计算平均质量
  }

  // 执行单位转换：将输入单位转换为LAMMPS内部单位
  convert_units_properly();

  // 根据用户选择设置PID参数
  if (manual_pid) {
    // 使用用户提供的PID参数
    // C++知识点：条件编译保护（comm->me == 0确保只在主进程输出）
    if (comm->me == 0) {  // 只在主进程输出信息（MPI并行）
      printf("使用手动PID参数: Kp=%.3f, Ki=%.3f, Kd=%.3f\n", pid_kp, pid_ki, pid_kd);
    }
  } else {
    // 自动计算最优PID参数
    auto_initialize_pid();
  }

  // 初始化压头输出数据
  indenter_flag = 0;             // 压头状态标志
  // C++知识点：数组初始化，一次性设置多个元素为0
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;  // [能量, fx, fy, fz]
}

/* ---------------------------------------------------------------------- */
/* 析构函数：对象销毁时调用，用于清理资源
   C++知识点：
   - 析构函数名前有~符号
   - 对象生命周期结束时自动调用
   - 这里是空函数，因为使用的都是自动管理的资源
*/

FixIndentLoadImproved::~FixIndentLoadImproved() {}

/* ---------------------------------------------------------------------- */
/* PID控制系统初始化函数
   功能：初始化所有PID控制相关的参数和变量
   调用时机：构造函数中调用一次
*/

void FixIndentLoadImproved::init_pid_parameters()
{
  // PID控制器状态变量初始化
  integral_error = 0.0;          // 积分误差累积
  previous_error = 0.0;          // 上一步的误差值
  max_integral_windup = target_load * 0.1;  // 积分饱和限制（防止积分项过大）
  
  // 控制系统参数
  pid_history_length = 20;       // 误差历史记录长度
  load_tolerance = 1.0;          // 载荷稳定判断容差（1%）
  adaptive_factor = 0.02;        // 自适应调节因子
  min_kp = 0.1;                  // 比例增益最小值
  max_kp = 5.0;                  // 比例增益最大值
  
  // 统计和监控变量
  control_step_count = 0;        // 控制步数计数器
  max_overshoot = 0.0;           // 最大超调量记录
  settling_time = 0.0;           // 建立时间
  load_stable = false;           // 载荷稳定状态标志
  
  // C++知识点：STL容器操作
  error_history.clear();         // 清空误差历史向量
  error_history.reserve(pid_history_length);  // 预分配内存空间，提高效率
}

/* ---------------------------------------------------------------------- */

void FixIndentLoadImproved::convert_units_properly()
{
  // 更精确的单位转换
  // metal单位系统：距离 Å, 时间 ps, 质量 g/mol
  // 速度从 m/s 转换为 Å/ps
  
  // 详细的单位转换计算：
  // 1 m/s = 10^10 Å/s = 10^10 Å/s * (1 ps / 10^12 s) = 10^-2 Å/ps
  const double conversion_factor = 1e-2; // m/s -> Å/ps
  
  // 安全检查：确保输入速度合理
  if (max_speed <= 0.0 || max_speed_overshoot <= 0.0) {
    error->all(FLERR, "Invalid speed values in fix indent/load/improved");
  }
  
  max_speed *= conversion_factor;
  max_speed_overshoot *= conversion_factor;
  sliding_speed *= conversion_factor;  // 滑动速度也要转换
  
  // 输出转换结果供调试
  if (comm->me == 0) {
    printf("单位转换完成：最大速度 = %.6f Å/ps, 超调速度 = %.6f Å/ps\n", 
           max_speed, max_speed_overshoot);
  }
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

/* setmask函数：告诉LAMMPS何时调用这个fix
   这是LAMMPS调用机制的核心！
   
   返回值：位掩码，指定在时间步的哪些阶段调用这个fix
   
   C++知识点：
   - 位运算操作 |= （按位或赋值）
   - 枚举常量作为位标志
   - 函数返回值用于系统回调注册
*/

int FixIndentLoadImproved::setmask()
{
  int mask = 0;  // 初始化为0
  
  // 设置调用时机的位标志
  mask |= POST_FORCE;        // 在力计算后调用post_force()
  mask |= POST_FORCE_RESPA;  // 在RESPA积分的力计算后调用post_force_respa()
  mask |= MIN_POST_FORCE;    // 在能量最小化的力计算后调用min_post_force()
  
  return mask;  // 返回位掩码给LAMMPS系统
  
  /* 工作原理解释：
     LAMMPS主循环会检查这个mask，在相应的时机调用对应的函数：
     
     时间步循环：
     1. initial_integrate()     // 初始积分
     2. force->compute()        // 计算原子间力
     3. post_force()           // ← 我们的函数在这里被调用！
     4. final_integrate()      // 最终积分
     
     这样我们的载荷控制就能在每个时间步都执行了！
  */
}

/* ---------------------------------------------------------------------- */
/* init函数：LAMMPS系统初始化时调用
   功能：设置RESPA多尺度积分的参数
   
   RESPA (Reversible reference System Propagation Algorithm)：
   - 一种多尺度时间积分方法
   - 允许不同相互作用使用不同的时间步长
   - 快速力（如键力）使用小时间步
   - 慢速力（如压头力）使用大时间步
   
   调用时机：在仿真开始前调用一次
*/

void FixIndentLoadImproved::init()
{
  // 检查是否使用RESPA积分方法
  if (utils::strmatch(update->integrate_style, "^respa")) {
    // 获取RESPA的最高层级（最慢的时间尺度）
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    // 如果用户指定了层级，使用较小的值
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */
/* setup函数：仿真开始时的设置函数
   功能：根据积分方法执行首次力计算
   
   参数说明：
   - vflag: virial计算标志（用于压力计算）
   
   工作原理：
   - Verlet积分：直接调用post_force()
   - RESPA积分：需要在正确的层级上调用，并复制力数据
   
   调用时机：仿真循环开始前调用一次
*/

void FixIndentLoadImproved::setup(int vflag)
{
  // 根据积分方法选择不同的设置方式
  if (utils::strmatch(update->integrate_style, "^verlet")) {
    // Verlet积分：直接调用力计算
    post_force(vflag);
  } else {
    // RESPA积分：需要在指定层级操作
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */
/* min_setup函数：能量最小化时的设置函数
   功能：为能量最小化过程设置压头力
   
   参数说明：
   - vflag: virial计算标志
   
   使用场景：
   - 在进行能量最小化（如弛豫）时调用
   - 确保压头力在优化过程中被正确计算
   
   调用时机：能量最小化开始前调用
*/

void FixIndentLoadImproved::min_setup(int vflag)
{
  // 直接调用力计算函数
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */
/* 核心函数：每个时间步的力计算后调用
   功能：
   1. 计算压头与原子的相互作用力
   2. 统计总载荷
   3. 执行PID载荷控制
   4. 更新压头位置
   
   调用时机：每个时间步都调用（由LAMMPS主循环自动调用）
   
   C++知识点：
   - vflag参数未使用，用注释避免编译器警告
   - 指针和二维数组的使用
*/

void FixIndentLoadImproved::post_force(int /*vflag*/)
{
  // 重置本时间步的压头力数据
  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;

  // 获取原子数据的指针
  // C++知识点：指向指针的指针（二维数组）
  double **x = atom->x;      // 原子位置数组 x[i][0]=x, x[i][1]=y, x[i][2]=z
  double **f = atom->f;      // 原子受力数组 f[i][0]=fx, f[i][1]=fy, f[i][2]=fz
  int *mask = atom->mask;    // 原子组掩码数组
  int nlocal = atom->nlocal; // 本进程的原子数量

  // 声明力计算中使用的局部变量
  double delx, dely, delz, r, dr, fmag, fx, fy, fz;
  double total_fx = 0.0, total_fy = 0.0, total_fz = 0.0;

  // 根据压头类型计算相互作用力
  if (istyle == CYLINDER) {
    // 圆柱形压头的力计算
    for (int i = 0; i < nlocal; i++) {
      // C++知识点：位运算检查原子是否属于指定组
      if (mask[i] & groupbit) {  // 检查原子是否在目标组中
        // 计算原子到圆柱轴线的距离
        delx = x[i][0] - cyl_x;  // x方向距离
        dely = x[i][1] - cyl_y;  // y方向距离
        r = sqrt(delx * delx + dely * dely);  // 径向距离

        // 如果原子在圆柱内部，计算排斥力
        if (r < cyl_radius) {
          dr = cyl_radius - r;     // 穿透深度
          fmag = k * dr * dr;      // 二次势：F = k * dr²

          // 安全检查：防止除零错误
          if (r > 1e-10) {  // 避免除以接近零的数
            // 计算力的各分量（径向排斥）
            fx = delx * fmag / r;    // x分量
            fy = dely * fmag / r;    // y分量
            fz = 0.0;                // 圆柱形压头z方向无力
          } else {
            // 当r接近零时，使用随机方向或特定方向
            fx = fmag * 0.707;  // 使用默认方向（约45度）
            fy = fmag * 0.707;
            fz = 0.0;
          }

          // 将力加到原子上
          f[i][0] += fx;
          f[i][1] += fy;

          // 累积总力（用于载荷计算）
          total_fx += fx;
          total_fy += fy;

          // 累积势能：U = ∫F·dr = k*dr³/3
          indenter[0] += k * dr * dr * dr / 3.0;
        }
      }
    }
  } else if (istyle == SPHERE) {
    // 球形压头的力计算
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        // 计算原子到球心的距离
        delx = x[i][0] - cyl_x;  // x方向距离
        dely = x[i][1] - cyl_y;  // y方向距离
        delz = x[i][2] - cyl_z;  // z方向距离
        r = sqrt(delx * delx + dely * dely + delz * delz);  // 球心距

        // 如果原子在球内，计算排斥力
        if (r < cyl_radius) {
          dr = cyl_radius - r;     // 穿透深度
          fmag = k * dr * dr;      // 二次势

          // 安全检查：防止除零错误
          if (r > 1e-10) {  // 避免除以接近零的数
            // 计算力的各分量（径向排斥）
            fx = delx * fmag / r;    // x分量
            fy = dely * fmag / r;    // y分量
            fz = delz * fmag / r;    // z分量
          } else {
            // 当r接近零时，使用随机方向
            fx = fmag * 0.577;  // 使用默认方向（约 1/√3）
            fy = fmag * 0.577;
            fz = fmag * 0.577;
          }

          // 将力加到原子上
          f[i][0] += fx;
          f[i][1] += fy;
          f[i][2] += fz;

          // 累积总力
          total_fx += fx;
          total_fy += fy;
          total_fz += fz;

          // 累积势能
          indenter[0] += k * dr * dr * dr / 3.0;
        }
      }
    }
  }

  // MPI并行：汇总所有进程的力
  double total_f[3];
  total_f[0] = total_fx;
  total_f[1] = total_fy;
  total_f[2] = total_fz;
  // C++知识点：MPI_Allreduce进行并行求和
  MPI_Allreduce(MPI_IN_PLACE, total_f, 3, MPI_DOUBLE, MPI_SUM, world);

  // 存储压头受到的总力（牛顿第三定律：作用力与反作用力）
  indenter[1] = -total_f[0];  // 压头受力fx（与原子受力相反）
  indenter[2] = -total_f[1];  // 压头受力fy
  indenter[3] = -total_f[2];  // 压头受力fz

  // 计算总载荷（力的大小）
  // C++知识点：使用sqrt计算向量模长
  current_load = sqrt(indenter[1] * indenter[1] + 
                     indenter[2] * indenter[2] + 
                     indenter[3] * indenter[3]);

  // 执行PID载荷控制（核心功能）
  control_load_pid();

  // 如果启用了滑动功能，更新压头滑动位置
  if (sliding_active) { 
    update_sliding_position(); 
  }
  
  // 更新统计信息和监控数据
  update_statistics();
}

/* ---------------------------------------------------------------------- */
/* PID载荷控制函数：系统的核心控制逻辑
   功能：使用PID算法自动调节压头位置以达到目标载荷
   
   PID控制原理：
   - P（比例）：响应当前误差
   - I（积分）：消除稳态误差  
   - D（微分）：预测未来误差，提高稳定性
   
   调用时机：每个时间步在post_force()中调用
*/

void FixIndentLoadImproved::control_load_pid()
{
  // 1. 计算控制误差
  double current_error = target_load - current_load;  // 误差 = 目标值 - 当前值
  control_step_count++;  // 控制步数计数器递增
  
  // 2. 计算微分项（误差变化率）
  double error_derivative = 0.0;
  if (control_step_count > 1 && dt > 1e-15) {  // 安全检查：防止时间步长过小
    // C++知识点：条件判断避免第一步计算微分
    error_derivative = (current_error - previous_error) / dt;
    
    // 限制微分项的值，防止过大的导数引起不稳定
    const double max_derivative = target_load / dt * 10.0;  // 最大允许的导数
    if (fabs(error_derivative) > max_derivative) {
      error_derivative = (error_derivative > 0) ? max_derivative : -max_derivative;
    }
  }
  
  // 3. 更新积分项（误差累积）
  if (dt > 1e-15) {  // 安全检查
    integral_error += current_error * dt;
  }
  
  // 4. 积分饱和限制（防止积分项过大导致系统不稳定）
  if (integral_error > max_integral_windup) integral_error = max_integral_windup;
  if (integral_error < -max_integral_windup) integral_error = -max_integral_windup;
  
  // 5. PID控制器输出计算
  // C++知识点：PID控制算法的标准公式
  double pid_output = pid_kp * current_error +      // 比例项
                     pid_ki * integral_error +      // 积分项
                     pid_kd * error_derivative;     // 微分项
  
  // 6. 限制最大调整速度（安全保护）
  double max_adjustment = max_speed * dt;
  if (fabs(pid_output) > max_adjustment) {
    // C++知识点：三元运算符 ? :
    pid_output = (pid_output > 0) ? max_adjustment : -max_adjustment;
  }
  
  // 7. 超调检测（载荷穿越目标值）
  bool overshoot_detected = false;
  if (error_history.size() >= 2) {
    double prev_error = error_history[error_history.size() - 1];
    
    // 安全检查：确保除数不为零
    if (fabs(target_load) > 1e-10) {
      // 检测误差符号变化（表示穿越了目标值）
      if ((current_error > 0 && prev_error < 0) || 
          (current_error < 0 && prev_error > 0)) {
        overshoot_detected = true;
        integral_error *= 0.7;  // 减少积分项以抑制振荡
        
        // 记录最大超调量用于性能评估
        double overshoot = fabs(current_error) / target_load * 100.0;
        if (overshoot > max_overshoot) max_overshoot = overshoot;
        
        // 添加超调警告（仅在严重超调时）
        if (overshoot > 20.0 && comm->me == 0) {
          printf("警告：检测到较大超调 %.1f%%，正在调整控制参数...\n", overshoot);
        }
      }
    }
  }
  
  // 8. 保守控制策略（在超调或大误差时降低控制强度）
  if (overshoot_detected || (fabs(target_load) > 1e-10 && fabs(current_error) > 0.1 * target_load)) {
    pid_output *= 0.5;  // 减少控制强度50%
  }
  
  // 9. 数值稳定性检查：防止NaN或无穷大值
  if (!std::isfinite(pid_output) || fabs(pid_output) > 1000.0) {
    pid_output = 0.0;  // 重置为安全值
    if (comm->me == 0) {
      printf("警告：PID输出异常，已重置为0\n");
    }
  }
  
  // 10. 更新压头位置（执行控制动作）
  // 注意：向下为负方向，载荷不足时压头下移增加载荷
  cyl_z -= pid_output;
  
  // 11. 更新历史数据
  previous_error = current_error;
  error_history.push_back(current_error);
  // C++知识点：STL容器大小控制
  if (error_history.size() > static_cast<size_t>(pid_history_length)) {
    error_history.erase(error_history.begin());  // 删除最旧的数据
  }
  
  // 12. 自适应PID参数调节
  if (enable_adaptive_pid) {
    adaptive_pid_adjustment();
  }
}

/* ---------------------------------------------------------------------- */

bool FixIndentLoadImproved::is_load_stable(double tolerance_percent)
{
  if (error_history.size() < 5) return false;
  
  // 安全检查：防止除零错误
  if (fabs(target_load) < 1e-10) {
    return fabs(current_load) < 1e-10;  // 如果目标载荷为0，检查当前载荷是否也为0
  }
  
  double max_recent_error = 0.0;
  size_t start_idx = error_history.size() - 5;
  for (size_t i = start_idx; i < error_history.size(); i++) {
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
    // 安全计算误差百分比
    double current_error_percent = 0.0;
    if (fabs(target_load) > 1e-10) {
      current_error_percent = fabs(target_load - current_load) / target_load * 100.0;
    }
    
    printf("Load Control Step %d: Target=%.2f nN, Current=%.2f nN, Error=%.2f%%, "
           "PID=(%.3f,%.3f,%.3f), Stable=%s\n",
           control_step_count, target_load, current_load, current_error_percent,
           pid_kp, pid_ki, pid_kd, load_stable ? "YES" : "NO");
  }
}

/* ---------------------------------------------------------------------- */
/* post_force_respa函数：RESPA积分中的力计算调用
   功能：在RESPA多尺度积分的指定层级调用力计算
   
   参数说明：
   - vflag: virial计算标志
   - ilevel: 当前RESPA层级
   - iloop: 循环索引（未使用，用注释避免编译器警告）
   
   工作原理：
   - 只在我们指定的RESPA层级(ilevel_respa)调用力计算
   - 这样压头力只在较慢的时间尺度上更新，提高效率
   
   调用时机：RESPA积分的每个时间步都会调用，但只在指定层级执行
*/

void FixIndentLoadImproved::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  // 只在指定的RESPA层级调用力计算
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */
/* min_post_force函数：能量最小化中的力计算
   功能：在能量最小化过程中计算压头力
   
   参数说明：
   - vflag: virial计算标志
   
   使用场景：
   - 在能量最小化的每次迭代中调用
   - 确保优化过程考虑压头的约束力
   
   调用时机：能量最小化的每次迭代后调用
*/

void FixIndentLoadImproved::min_post_force(int vflag)
{
  // 直接调用标准力计算函数
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */
/* compute_scalar函数：返回标量输出值
   功能：提供压头的势能给LAMMPS输出系统
   
   返回值：压头与原子相互作用的总势能
   
   使用方法：
   - 在LAMMPS输入文件中可以用 f_<fix_ID> 访问
   - 例如：thermo_style custom step temp press f_indent
   
   调用时机：每次输出热力学数据时调用
*/

double FixIndentLoadImproved::compute_scalar()
{
  return indenter[0];    // 返回存储的势能值
}

/* ---------------------------------------------------------------------- */
/* compute_vector函数：返回矢量输出值
   功能：提供压头受到的力的各分量给LAMMPS输出系统
   
   参数说明：
   - n: 分量索引（0=fx, 1=fy, 2=fz）
   
   返回值：压头受到的力的指定分量
   
   使用方法：
   - 在LAMMPS输入文件中可以用 f_<fix_ID>[n] 访问
   - 例如：variable fz equal f_indent[3]  # 获取z方向力
   
   调用时机：每次输出热力学数据时调用
*/

double FixIndentLoadImproved::compute_vector(int n)
{
  return indenter[n + 1];    // 返回存储的力分量（索引偏移1，因为[0]是势能）
}

/* ---------------------------------------------------------------------- */
/* get_indenter_type函数：解析压头类型字符串
   功能：将用户输入的字符串转换为内部枚举值
   
   参数说明：
   - arg: 用户输入的压头类型字符串
   
   支持的类型：
   - "cylinder": 圆柱形压头，适用于线接触
   - "sphere": 球形压头，适用于点接触
   
   返回值：对应的枚举值（CYLINDER或SPHERE）
   
   异常处理：如果输入不支持的类型，调用error->all()终止程序
*/

int FixIndentLoadImproved::get_indenter_type(const char *arg)
{
  // C++知识点：strcmp()函数用于比较C风格字符串
  if (strcmp(arg, "cylinder") == 0) return CYLINDER;  // 圆柱形
  if (strcmp(arg, "sphere") == 0) return SPHERE;      // 球形
  
  // 遇到不支持的类型时报错并终止程序
  error->all(FLERR, "Illegal fix indent/load/improved type: {}", arg);
  return NONE;  // 这行不会执行，但保留以避免编译器警告
}

/* ---------------------------------------------------------------------- */
/* update_sliding_position函数：更新压头滑动位置
   功能：实现压头在x方向的滑动运动
   
   应用场景：
   - 模拟划痕过程（横向滑动）
   - 摩擦试验（滑动接触）
   - 复杂加载条件下的纳米压痕
   
   计算公式：
   新位置 = 当前位置 + 方向 × 速度 × 时间步长
   
   参数说明：
   - sliding_dir: 滑动方向（+1或-1）
   - sliding_speed: 滑动速度（已转换为Å/ps）
   - dt: 时间步长（ps）
   
   调用时机：每个时间步在post_force()中调用（当sliding_active=1时）
*/

void FixIndentLoadImproved::update_sliding_position()
{
  // 更新压头在x方向的位置
  // C++知识点：使用复合赋值运算符 +=
  cyl_x += sliding_dir * sliding_speed * dt;
  
  // 注意：这里没有边界检查，如果需要可以添加：
  // if (cyl_x > domain->boxhi[0] || cyl_x < domain->boxlo[0]) sliding_dir *= -1;
}