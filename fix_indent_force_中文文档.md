# LAMMPS 力控制压痕功能详细说明文档

## 目录
1. [功能概述](#功能概述)
2. [命令语法](#命令语法)
3. [参数说明](#参数说明)
4. [使用示例](#使用示例)
5. [算法原理](#算法原理)
6. [高级功能](#高级功能)
7. [变量访问](#变量访问)
8. [力的符号约定](#力的符号约定)
9. [最佳实践](#最佳实践)
10. [故障排除](#故障排除)
11. [性能优化](#性能优化)
12. [完整仿真示例](#完整仿真示例)

---

## 功能概述

`fix indent/force` 是LAMMPS中的**力控制压痕**功能，扩展了标准的 `fix indent` 命令。与传统的位置控制压痕不同，该功能允许用户指定**目标力值**（fx, fy, fz），系统会自动调整压头位置以达到指定的力值。

### 主要特点

- ✅ **力控制模式**：指定目标力而非位置
- ✅ **多方向控制**：独立控制x、y、z三个方向的力
- ✅ **自动收敛**：内置迭代算法确保力值收敛
- ✅ **MPI并行支持**：正确处理多处理器并行计算
- ✅ **实时位置访问**：可获取压头的实时位置信息
- ✅ **高级算法**：自适应步长、振荡阻尼、梯度加速

### 适用场景

1. **载荷控制压痕实验** - 模拟真实实验中的力控制条件
2. **纳米压痕仿真** - 更贴近实验的力反馈机制
3. **接触力学研究** - 维持特定接触力的研究
4. **材料测试仿真** - 力控制加载的材料性能测试
5. **多轴加载研究** - 复杂应力状态下的材料响应

---

## 命令语法

### 基本语法

```lammps
fix ID group-ID indent/force K geometry fx fy fz [keyword value ...]
```

### 参数详解

| 参数 | 说明 | 单位 |
|------|------|------|
| `ID` | Fix标识符 | - |
| `group-ID` | 原子组名称 | - |
| `K` | 力常数 | 力/距离² |
| `geometry` | 几何形状（sphere或plane） | - |
| `fx, fy, fz` | 目标力分量 | 力单位 |

### 支持的几何形状

#### 1. 球形压头 (sphere)
```lammps
fix ID group indent/force K sphere x y z R fx fy fz
```
- `x, y, z`: 球心坐标（距离单位）
- `R`: 球半径（距离单位）

#### 2. 平面压头 (plane)
```lammps
fix ID group indent/force K plane dim pos side fx fy fz
```
- `dim`: 平面法向（x、y或z）
- `pos`: 平面位置（距离单位）
- `side`: 压痕方向（lo或hi）

### 可选关键字

| 关键字 | 值 | 说明 |
|--------|----|----- |
| `side` | `in`/`out` | 内部/外部压痕（仅球形）|
| `units` | `lattice`/`box` | 晶格单位/盒子单位 |

---

## 参数说明

### 力常数 K

力常数K与标准indent fix相同，定义了压头与原子间的相互作用强度：

$$F(r) = -K(r-R)^2$$

其中：
- `r`: 原子到压头表面的距离
- `R`: 特征尺寸（球半径或平面位置）

**建议值**：
- 小系统（< 1000原子）：K = 10-100
- 中等系统（1000-10000原子）：K = 100-1000  
- 大系统（> 10000原子）：K = 1000-10000

### 目标力 fx, fy, fz

目标力分量可以是：
- **数值常量**：`-100.0`（压缩100个力单位）
- **变量**：`v_force_z`（使用变量动态控制）

**力的方向**：
- **正值** = 拉伸（压头拉扯材料）
- **负值** = 压缩（压头推入材料，常用于压痕）

---

## 使用示例

### 1. 基础球形压痕

```lammps
# 系统设置
dimension 3
boundary p p p
atom_style atomic

# 创建原子、设置势能等...
# （标准LAMMPS设置）

# 应用力控制球形压头
# 目标：z方向压缩100个力单位
fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -100.0

# 监控力和位置
variable fx equal f_1[1]    # x方向力
variable fy equal f_1[2]    # y方向力  
variable fz equal f_1[3]    # z方向力
variable ix equal f_1[4]    # 压头x位置
variable iy equal f_1[5]    # 压头y位置
variable iz equal f_1[6]    # 压头z位置

thermo_style custom step temp v_fx v_fy v_fz v_ix v_iy v_iz
thermo 100

run 10000
```

### 2. 基础平面压痕

```lammps
# 从顶部（z方向）推入的平面压头
fix 2 all indent/force 15.0 plane z 20.0 hi 0.0 0.0 -50.0

# 监控相关力分量
variable fz equal f_2[3]    # z方向力
variable pz equal f_2[6]    # 平面位置

thermo_style custom step temp v_fz v_pz
run 5000
```

### 3. 渐进加载（避免系统冲击）

```lammps
# 定义渐进加载函数
variable max_force equal -200.0
variable ramp_steps equal 5000
variable current_force equal "min(v_max_force*step/v_ramp_steps, v_max_force)"

fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 v_current_force

# 输出力的演变
fix output all print 100 "${step} ${current_force}" file force_ramp.dat
thermo_style custom step v_current_force v_fx v_fy v_fz
```

### 4. 循环加载

```lammps
# 正弦波力变化
variable frequency equal 0.001    # 每时间步的循环数
variable amplitude equal 50.0     # 振幅
variable baseline equal -100.0    # 基线
variable cyclic_force equal "v_baseline + v_amplitude*sin(2*PI*v_frequency*step)"

fix 1 all indent/force 8.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 v_cyclic_force

# 跟踪力和位置演变
variable fx equal f_1[1]
variable fy equal f_1[2]
variable fz equal f_1[3]
variable iz equal f_1[6]

fix output all print 50 "${step} ${cyclic_force} ${fz} ${iz}" file cyclic_test.dat
```

### 5. 多方向同时加载

```lammps
# 定义复杂力模式
variable fx_target equal "10.0*sin(step*dt*0.1)"     # x方向正弦力
variable fy_target equal "5.0*cos(step*dt*0.05)"    # y方向余弦力  
variable fz_target equal "-50.0 - 20.0*step/10000.0" # z方向渐进压缩

fix 1 all indent/force 12.0 sphere 5.0 5.0 7.0 1.0 v_fx_target v_fy_target v_fz_target

# 监控所有分量
variable fx equal f_1[1]
variable fy equal f_1[2]
variable fz equal f_1[3]

thermo_style custom step v_fx_target v_fx v_fy_target v_fy v_fz_target v_fz
```

### 6. 移动压头与力控制组合

```lammps
# 移动的压头中心
variable center_x equal "5.0 + 2.0*step/10000.0"
variable center_y equal "5.0"
variable center_z equal "7.0"

# 恒定向下力
fix 1 all indent/force 10.0 sphere v_center_x v_center_y v_center_z 1.0 0.0 0.0 -75.0

# 跟踪压头位置演变
variable ix equal f_1[4]
variable iy equal f_1[5]
variable iz equal f_1[6]

fix output all print 100 "${step} ${center_x} ${ix} ${iy} ${iz}" file moving_indenter.dat
```

### 7. 使用变量的高级示例

```lammps
# 温度相关的力
compute sys_temp all temp
variable temp_force equal "-100.0*(c_sys_temp/300.0)"

# 时间和位置相关的复杂力模式
variable time_factor equal "1.0 + 0.5*sin(step*dt*0.02)"
variable complex_fz equal "v_temp_force*v_time_factor"

fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 v_complex_fz

# 输出详细信息
thermo_style custom step c_sys_temp v_time_factor v_temp_force v_complex_fz v_fz
```

---

## 算法原理

### 核心算法流程

力控制压痕算法在每个时间步执行**内部迭代循环**：

```
1. 力计算 → 使用标准indent物理计算当前压头力
2. MPI求和 → MPI_Allreduce在所有处理器间求和力
3. 误差计算 → 比较实际力 (fx, fy, fz) 与目标力
4. 收敛检查 → 如果误差在容差内则停止
5. 位置调整 → 基于力误差使用自适应算法调整压头位置
6. 重复 → 直到收敛或达到最大迭代次数
```

### 物理原理

压头与原子的相互作用遵循与标准indent相同的力表达式：

$$F(r) = -K(r-R)^2$$

其中：
- **K**: 力常数（用户指定）
- **r**: 原子到压头表面的距离
- **R**: 特征尺寸（球半径或平面位置）

**关键创新**：通过调整R值（压头位置）来控制总力，而不是直接指定位置。

### 收敛判据

算法使用相对误差作为收敛判据：

$$\text{收敛条件}: \frac{|F_{\text{actual}} - F_{\text{target}}|}{|F_{\text{target}}|} < \text{tolerance}$$

- **默认容差**: 0.01 (1%)
- **零力处理**: 当目标力接近零时，使用绝对误差

### 高级算法特性

#### 1. 六级自适应步长控制

根据相对误差大小动态调整位置调整步长：

| 相对误差范围 | 步长因子 | 策略 |
|-------------|----------|------|
| > 200% | 0.002 | 极大误差，保守步长 |
| > 100% | 0.001 | 大误差，小步长 |
| > 50% | 0.0008 | 中大误差，较小步长 |
| > 10% | 0.0005 | 中等误差，保守步长 |
| > 1% | 0.0003 | 小误差，很小步长 |
| ≤ 1% | 0.0001 | 极小误差，微小步长 |

```cpp
double adaptive_factor = get_adaptive_factor(error, target);
double adjustment = adaptive_factor * error;
```

#### 2. 动态位置限制

位置调整被动态限制以防止系统不稳定：

$$\text{限制} = \text{基础限制} \times (1 + 0.5 \times \text{相对误差})$$

- **基础限制**: 0.05（保守）
- **最大限制**: 0.15
- **目的**: 防止大幅位置跳跃导致的系统不稳定

#### 3. 振荡阻尼机制

检测连续调整方向是否交替，如果检测到振荡：

```cpp
if (current_adjustment * previous_adjustment < 0) {
    adjustment *= 0.5;  // 将步长减少50%
}
```

**作用**: 防止在最优位置附近的振荡行为。

#### 4. 基于梯度的加速

使用保守牛顿法加速收敛：

```cpp
double gradient = (current_force - prev_force) / (current_pos - prev_pos);
if (abs(gradient) > 1e-10) {
    double newton_step = error / gradient;
    adjustment = 0.8 * newton_step;  // 80%缩放，保守起见
}
```

**条件**: 仅当牛顿步长小于自适应步长时使用。

#### 5. 历史跟踪

算法维护历史信息以支持高级策略：

- **前一位置**: 用于梯度估计
- **前一力值**: 用于梯度计算  
- **前一调整**: 用于振荡检测
- **首次迭代标志**: 特殊处理初始迭代

### 算法参数

当前算法使用以下硬编码参数：

| 参数 | 值 | 说明 |
|------|----|----- |
| `tolerance` | 0.01 | 1%力容差 |
| `max_iterations` | 500 | 每时间步最大迭代数（增强版）|
| `base_adjustment` | 动态 | 0.0001到0.002变化 |
| `oscillation_damping` | 0.5 | 振荡检测时的阻尼因子 |
| `newton_scaling` | 0.8 | 牛顿法的保守缩放 |

### MPI并行处理

**关键修复**: 原始实现的主要bug是在并行环境中力计算错误。

**解决方案**:
```cpp
// 在每次迭代中进行MPI力求和
double indenter_global[4];
MPI_Allreduce(indenter, indenter_global, 4, MPI_DOUBLE, MPI_SUM, world);

// 使用全局力计算误差
fx_error = target_fx - indenter_global[1];
fy_error = target_fy - indenter_global[2]; 
fz_error = target_fz - indenter_global[3];
```

**重要性**: 这确保了无论处理器数量如何，都能产生一致的结果。

### 初始接触处理

当目标力非零但当前力接近零时（无接触），算法自动：

1. **计算原子组质心**
2. **向质心移动压头** (移动因子 = 0.1)
3. **建立初始接触**

这防止了压头"滞留"在远离材料的位置。

---

## 高级功能

### 1. 静默运行模式

力控制算法默认静默运行，不输出详细的收敛信息。所有必要信息都通过fix变量访问。

**优势**:
- ✅ 更清洁的输出
- ✅ 更容易分析结果  
- ✅ 用户完全控制输出内容

**推荐监控方式**:
```lammps
thermo_style custom step temp v_fx v_fy v_fz v_ix v_iy v_iz
thermo 100
```

### 2. 可扩展几何形状

当前实现支持：
- ✅ **球形压头** (sphere) - 完全实现
- ✅ **平面压头** (plane) - 完全实现
- ❌ **圆柱形压头** (cylinder) - 待实现
- ❌ **圆锥形压头** (cone) - 待实现

### 3. 能量计算支持

支持 `fix_modify energy` 选项：

```lammps
fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -100.0
fix_modify 1 energy yes

# 压头能量将包含在系统总势能中
thermo_style custom step temp pe ke etotal
```

### 4. r-RESPA集成支持

支持多时间步长集成：

```lammps
run_style respa 4 2 inner 1 4.0 middle 2 8.0 outer 1
fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -100.0
fix_modify 1 respa 3  # 在第3层级应用力
```

---

## 变量访问

### 6分量向量访问

`fix indent/force` 提供6个分量的向量访问：

| 索引 | 球形压头 | 平面压头 | 说明 |
|------|----------|----------|------|
| `[1]` | fx | fx | x方向力 |
| `[2]` | fy | fy | y方向力 |
| `[3]` | fz | fz | z方向力 |
| `[4]` | x位置 | x位置(固定) | 压头x坐标 |
| `[5]` | y位置 | y位置(固定) | 压头y坐标 |
| `[6]` | z位置 | 平面位置 | 压头z坐标或平面位置 |

### 典型输出示例

```bash
Step          Temp           v_fx           v_fy           v_fz          v_ix           v_iy           v_iz
     0   0.16999953     0.00039873     0.0027783      99.17299      5.0000         5.0000         7.0000
   100   0.1            0.0018872     -0.0050041      100.01008     5.0001         4.9998         6.9995
   200   0.1            0.0021303     -0.0043898      99.996308     5.0002         4.9999         6.9993
```

### 应用示例

```lammps
# 定义所有变量
variable fx equal f_1[1]    # x方向力
variable fy equal f_1[2]    # y方向力
variable fz equal f_1[3]    # z方向力
variable ix equal f_1[4]    # x位置
variable iy equal f_1[5]    # y位置
variable iz equal f_1[6]    # z位置

# 计算力的大小
variable force_mag equal "sqrt(v_fx*v_fx + v_fy*v_fy + v_fz*v_fz)"

# 计算位置变化
variable pos_change equal "sqrt((v_ix-5.0)^2 + (v_iy-5.0)^2 + (v_iz-7.0)^2)"

# 输出到文件
fix output all print 100 "${step} ${fx} ${fy} ${fz} ${force_mag} ${pos_change}" file analysis.dat
```

---

## 力的符号约定

### 物理意义

力的符号遵循牛顿第三定律 - 输出的力是**作用在压头上的反作用力**：

```lammps
# 压缩（常见压痕）
fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -100.0
#                                                              ^^^
#                                                       负值 = 压缩
# 物理解释：压头对材料向下施力，材料对压头向上反作用力为+100
# 但我们指定-100，表示希望压头受到向下的反作用力（即压头向下压）

# 拉伸（特殊应用）
fix 2 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 +100.0
#                                                              ^^^  
#                                                       正值 = 拉伸
# 物理解释：压头从材料中"拉出"，压头受到向上的反作用力
```

### 符号对照表

| 应用场景 | 力方向 | 符号 | 物理含义 |
|---------|-------|------|---------|
| 标准压痕 | 向下压 | 负值 | 压缩应力 |
| 拉拔试验 | 向上拉 | 正值 | 拉伸应力 |
| 侧向推入 | x/y方向 | ±值 | 剪切应力 |

### 验证符号的方法

```lammps
# 设置压缩测试
fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -50.0

variable fz equal f_1[3]
variable iz equal f_1[6]

# 输出监控
fix check all print 10 "Step: ${step} Target: -50.0 Actual: ${fz} Position: ${iz}" screen yes

# 预期结果：
# - fz应该接近-50.0（负值确认压缩）
# - iz应该减小（压头向下移动）
```

---

## 最佳实践

### 1. 系统准备

**步骤1：能量最小化**
```lammps
minimize 1.0e-6 1.0e-8 1000 10000
```

**步骤2：热平衡**
```lammps
velocity all create 300.0 12345
fix nvt all nvt temp 300.0 300.0 100.0
run 20000
```

**步骤3：切换到力控制**
```lammps
unfix nvt
fix nve all nve
fix indent all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -50.0
```

### 2. 力大小指导

根据系统规模选择合适的力值：

| 系统规模 | 原子数 | 推荐力范围 | 备注 |
|---------|--------|------------|------|
| 小型 | < 1,000 | 1-10 | 快速测试 |
| 中型 | 1,000-10,000 | 10-100 | 常见研究 |
| 大型 | > 10,000 | 100-1,000 | 产业应用 |

**经验法则**: 力值应与系统的特征能量尺度匹配。

### 3. 时间步长建议

```lammps
# 标准平衡阶段
timestep 0.001

# 力控制阶段（推荐减小）
timestep 0.0005   # 或更小，如0.0002

fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -100.0
```

**原因**: 力控制需要多次迭代，较小的时间步长有助于稳定性。

### 4. 温度控制

力控制过程中监控和控制温度：

```lammps
# 方法1：Langevin恒温器（推荐）
compute mobile_temp mobile temp
fix thermostat mobile langevin 300.0 300.0 100.0 48279

# 方法2：速度重标（温和）
fix temp_rescale all temp/rescale 100 300.0 300.0 0.02 1.0

# 方法3：Nosé-Hoover（动力学保守）
fix nvt mobile nvt temp 300.0 300.0 100.0
```

### 5. 渐进加载策略

避免系统冲击的渐进加载：

```lammps
# 线性渐进
variable max_force equal -200.0
variable ramp_time equal 10000
variable current_force equal "min(v_max_force*step/v_ramp_time, v_max_force)"

# 指数渐进（更温和）
variable exp_force equal "v_max_force*(1.0 - exp(-step/v_ramp_time))"

# 三角波渐进
variable tri_force equal "v_max_force*min(step/v_ramp_time, 1.0)"

fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 v_current_force
```

### 6. 监控策略

设置全面的监控：

```lammps
# 基本变量
variable fx equal f_1[1]
variable fy equal f_1[2]
variable fz equal f_1[3]
variable ix equal f_1[4]
variable iy equal f_1[5]
variable iz equal f_1[6]

# 系统诊断
compute mobile_temp mobile temp
compute mobile_press mobile pressure mobile_temp

# 误差分析
variable target_fz equal -100.0
variable error_fz equal "abs(v_fz - v_target_fz)"
variable rel_error_fz equal "v_error_fz/abs(v_target_fz)*100"

# 输出设置
thermo_style custom step c_mobile_temp c_mobile_press v_target_fz v_fz v_error_fz v_rel_error_fz v_iz
thermo 100

fix detailed_output all print 100 "${step} ${fx} ${fy} ${fz} ${ix} ${iy} ${iz} ${error_fz} ${rel_error_fz}" file detailed.dat
```

---

## 故障排除

### 常见问题1：力不收敛

**症状**: 
- 力误差很大且持续振荡
- `thermo`输出显示力值远离目标

**可能原因及解决方案**:

```lammps
# 原因1：目标力过大
# 解决：减小目标力
variable small_force equal -10.0  # 从小力开始
fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 v_small_force

# 原因2：时间步长过大  
# 解决：减小时间步长
timestep 0.0001

# 原因3：系统未充分平衡
# 解决：延长平衡时间
fix nvt all nvt temp 300.0 300.0 100.0
run 50000  # 增加平衡步数
```

### 常见问题2：系统过热

**症状**:
- 温度快速上升
- 原子"飞出"系统
- 能量发散

**解决方案**:
```lammps
# 添加强制恒温器
compute mobile_temp mobile temp
fix emergency_thermo mobile langevin 300.0 300.0 50.0 12345

# 减小力和时间步长
variable conservative_force equal -5.0
timestep 0.00005

fix 1 all indent/force 5.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 v_conservative_force

# 监控温度
thermo_modify temp mobile_temp
```

### 常见问题3：无初始接触

**症状**:
- 力始终为零
- 压头位置不变

**解决方案**:
```lammps
# 将压头移近材料
fix 1 all indent/force 10.0 sphere 5.0 5.0 6.0 1.0 0.0 0.0 -50.0
#                                            ^^^
#                                        减小z坐标

# 或增加压头半径
fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 2.0 0.0 0.0 -50.0
#                                               ^^^
#                                           增大半径

# 算法会自动寻找接触，但初始位置要合理
```

### 常见问题4：性能问题

**症状**:
- 仿真非常慢
- CPU使用率低

**原因**: 算法可能无法收敛，每步都运行满500次迭代。

**诊断**:
```lammps
# 添加误差监控
variable error_fz equal "abs(v_fz - v_target_fz)"
fix error_monitor all print 100 "${step} Error: ${error_fz}" screen yes

# 如果误差始终很大，说明收敛困难
```

**解决方案**:
```lammps
# 1. 减小目标力
# 2. 使用渐进加载
# 3. 改善初始条件
# 4. 检查系统稳定性
```

### 常见问题5：并行结果不一致

**症状**:
- 单核和多核结果不同
- 不同处理器数量结果差异

**解决方案**:
此问题已在当前版本修复。如果仍出现，请检查：

```lammps
# 确保使用相同的随机种子
velocity all create 300.0 12345  # 固定种子

# 确保系统充分平衡
run 20000  # 足够的平衡步数

# 监控不同核心数的结果一致性
```

---

## 性能优化

### 1. 收敛性优化

**监控收敛性能**:
```lammps
# 创建收敛诊断变量
variable conv_fx equal "abs(v_fx/v_target_fx - 1.0)*100"
variable conv_fy equal "abs(v_fy/v_target_fy - 1.0)*100"  
variable conv_fz equal "abs(v_fz/v_target_fz - 1.0)*100"

fix convergence_monitor all print 50 "${step} Convergence: ${conv_fx}% ${conv_fy}% ${conv_fz}%" file convergence.dat

# 良好收敛：误差 < 1%
# 可接受收敛：误差 < 5%
# 收敛困难：误差 > 10%
```

### 2. 系统设置优化

**优化压头参数**:
```lammps
# 合适的力常数（不要太大也不要太小）
variable k_optimal equal "100.0/xlat/xlat"  # 基于晶格尺度

# 合适的压头大小（不要太小）
variable radius_min equal "2.0*xlat"  # 至少2个晶格常数

fix 1 all indent/force ${k_optimal} sphere 5.0 5.0 7.0 ${radius_min} 0.0 0.0 -50.0
```

**优化原子分布**:
```lammps
# 确保压头附近有足够原子密度
# 避免压头在空洞区域开始

# 检查原子分布
compute coord_num all coord/atom cutoff 3.0
dump distribution all custom 1 distribution.dump id x y z c_coord_num

run 0  # 仅输出初始分布
undump distribution
```

### 3. 内存和I/O优化

**减少输出频率**:
```lammps
# 高频监控（调试用）
thermo 10
fix output all print 10 "..." file debug.dat

# 生产运行（优化性能）
thermo 1000  
fix output all print 500 "..." file production.dat
```

**优化dump输出**:
```lammps
# 仅输出必要信息
dump trajectory all atom 1000 trajectory.lammpstrj
dump_modify trajectory every 1000 first yes

# 或使用更紧凑格式
dump custom_out all custom 1000 output.dump id x y z fx fy fz
```

---

## 完整仿真示例

### 纳米压痕完整示例

```lammps
# ===================================================================
# LAMMPS 力控制纳米压痕完整仿真示例
# 作者：[您的姓名]
# 日期：[日期]
# 描述：展示完整的力控制压痕仿真流程
# ===================================================================

# 1. 基本设置
dimension 3
boundary p p p
atom_style atomic

# 2. 创建基底
lattice fcc 3.5
region substrate block 0 30 0 30 0 20
create_box 2 substrate
create_atoms 1 substrate

# 3. 设置质量和势能
mass 1 1.0  # 可动原子
mass 2 1.0  # 固定原子  
pair_style lj/cut 8.0
pair_coeff * * 1.0 3.5 8.0

# 4. 创建固定底层
region bottom block INF INF INF INF 0 3
group bottom region bottom
group mobile subtract all bottom
set group bottom type 2

# 5. 初始条件
velocity mobile create 300.0 12345
fix 1 all nve
fix 2 bottom setforce 0.0 0.0 0.0
fix 3 mobile langevin 300.0 300.0 100.0 54321

# 6. 能量最小化
minimize 1.0e-6 1.0e-8 1000 10000
print "Energy minimization completed."

# 7. 热平衡
timestep 0.001
thermo 2000
run 20000
print "Thermal equilibration completed."

# 8. 准备力控制压痕
unfix 3
fix 3 mobile langevin 300.0 300.0 100.0 54321  # 保持温控
timestep 0.0005  # 减小时间步长

# 9. 定义渐进加载
variable max_force equal -300.0
variable ramp_time equal 30000
variable hold_time equal 10000
variable total_time equal "v_ramp_time + v_hold_time"

# 渐进加载函数
variable current_force equal "step < v_ramp_time ? v_max_force*step/v_ramp_time : v_max_force"

# 10. 应用力控制压头
fix indent all indent/force 60.0 sphere 15.0 15.0 25.0 3.0 0.0 0.0 v_current_force

# 11. 定义监控变量
variable fx equal f_indent[1]
variable fy equal f_indent[2]
variable fz equal f_indent[3]
variable ix equal f_indent[4]
variable iy equal f_indent[5]
variable iz equal f_indent[6]

# 系统诊断变量
compute mobile_temp mobile temp
compute mobile_press mobile pressure mobile_temp
variable force_error equal "abs(v_fz - v_current_force)"
variable rel_error equal "v_force_error/abs(v_current_force + 1e-10)*100"
variable penetration equal "25.0 - v_iz"  # 压入深度

# 12. 输出设置
thermo_style custom step c_mobile_temp c_mobile_press v_current_force v_fx v_fy v_fz v_penetration v_rel_error
thermo 500

# 详细数据输出
fix detailed_log all print 200 "${step} ${current_force} ${fx} ${fy} ${fz} ${ix} ${iy} ${iz} ${penetration} ${force_error} ${rel_error}" file nanoindentation_detailed.dat

# 轨迹输出（较低频率以节省空间）
dump trajectory all atom 2000 nanoindentation_trajectory.lammpstrj
dump_modify trajectory first yes

# 13. 执行仿真
print "Starting force-controlled nanoindentation..."
print "Target force: ${max_force} units"
print "Ramp time: ${ramp_time} steps"
print "Hold time: ${hold_time} steps"

run ${total_time}

# 14. 卸载阶段（可选）
variable unload_time equal 15000
variable unload_force equal "v_max_force*(1.0 - step/v_unload_time)"

# 重新定义力
variable final_force equal "step < v_unload_time ? v_unload_force : 0.0"
fix_modify indent fx 0.0 fy 0.0 fz v_final_force  # 如果支持动态修改

print "Starting unloading phase..."
run ${unload_time}

# 15. 结果汇总
print "============ SIMULATION COMPLETED ============"
print "Final Results:"
print "  Final applied force: ${fz} (target was ${current_force})"
print "  Maximum penetration: ${penetration}"
print "  Final indenter position: (${ix}, ${iy}, ${iz})"
print "  Final force error: ${rel_error}%"
print "  Final system temperature: ${mobile_temp}"

# 16. 清理
undump trajectory
unfix detailed_log
unfix indent

# 17. 自由松弛（可选，观察材料恢复）
print "Starting free relaxation..."
run 10000

print "Nanoindentation simulation completed successfully!"
print "Data saved to: nanoindentation_detailed.dat"
print "Trajectory saved to: nanoindentation_trajectory.lammpstrj"

# ===================================================================
# 分析建议：
# 1. 绘制力-位移曲线：current_force vs penetration
# 2. 分析收敛性：监控 rel_error
# 3. 检查系统稳定性：监控温度和压力
# 4. 可视化：使用OVITO或VMD分析trajectory文件
# ===================================================================
```

### 运行命令

```bash
# 编译LAMMPS（如果需要）
cd lammps/build
make -j4

# 运行仿真
./lmp -in nanoindentation_complete.in -log nanoindentation.log

# 分析结果
python analyze_results.py nanoindentation_detailed.dat
```

### 后处理脚本示例（Python）

```python
import numpy as np
import matplotlib.pyplot as plt

# 读取数据
data = np.loadtxt('nanoindentation_detailed.dat')
step = data[:, 0]
target_force = data[:, 1]
actual_force_z = data[:, 4]  
penetration = data[:, 8]
force_error = data[:, 9]
rel_error = data[:, 10]

# 绘制结果
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))

# 力-时间曲线
ax1.plot(step, target_force, 'b-', label='Target Force')
ax1.plot(step, actual_force_z, 'r--', label='Actual Force')
ax1.set_xlabel('Time Step')
ax1.set_ylabel('Force')
ax1.set_title('Force Control Performance')
ax1.legend()
ax1.grid(True)

# 力-位移曲线
ax2.plot(penetration, -actual_force_z, 'g-', linewidth=2)
ax2.set_xlabel('Penetration Depth')
ax2.set_ylabel('Applied Force (Magnitude)')
ax2.set_title('Force-Displacement Curve')
ax2.grid(True)

# 收敛误差
ax3.semilogy(step, rel_error, 'orange')
ax3.set_xlabel('Time Step')
ax3.set_ylabel('Relative Error (%)')
ax3.set_title('Convergence Performance')
ax3.grid(True)

# 累积统计
ax4.hist(rel_error, bins=50, alpha=0.7, color='purple')
ax4.set_xlabel('Relative Error (%)')
ax4.set_ylabel('Frequency')
ax4.set_title('Error Distribution')
ax4.grid(True)

plt.tight_layout()
plt.savefig('nanoindentation_analysis.png', dpi=300)
plt.show()

# 输出统计信息
print(f"Average relative error: {np.mean(rel_error):.2f}%")
print(f"Maximum penetration: {np.max(penetration):.3f}")
print(f"Final force accuracy: {rel_error[-1]:.2f}%")
```

---

## 总结

本文档详细介绍了LAMMPS中`fix indent/force`功能的各个方面，从基本使用到高级算法原理。该功能通过先进的力控制算法，为分子动力学仿真提供了更贴近实验的载荷控制能力。

### 关键优势

1. **实验相关性** - 模拟真实实验中的力控制条件
2. **算法先进性** - 500次迭代、自适应步长、振荡阻尼等
3. **并行支持** - MPI并行环境中的正确实现
4. **灵活监控** - 完整的力和位置信息访问
5. **稳定性** - 保守参数设计确保系统稳定

### 使用要点

- **充分准备系统** - 能量最小化和热平衡
- **合理选择参数** - 力大小、时间步长、温度控制
- **密切监控** - 收敛性、稳定性、准确性
- **渐进加载** - 避免系统冲击
- **结果验证** - 检查物理合理性

希望本文档能帮助用户成功应用力控制压痕功能，推进相关科学研究和工程应用。

---

## 参考资料

1. LAMMPS 官方文档：https://docs.lammps.org/
2. fix indent 原始功能：https://docs.lammps.org/fix_indent.html
3. 相关研究论文和技术报告
4. 本项目源代码：`src/WANG/fix_indent_force.cpp`

---

*本文档最后更新：2025年1月*  
*版本：1.0*