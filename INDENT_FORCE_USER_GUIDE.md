# Fix indent/force User Guide

## Introduction

The `fix indent/force` command provides force-controlled indentation for LAMMPS, extending the standard position-controlled `fix indent`. This guide provides practical information for using this advanced feature effectively.

## Quick Start

### Basic Spherical Indentation

```lammps
# Setup system first
dimension 3
boundary p p p
atom_style atomic

# Create atoms, set potentials, etc.
# ... (standard LAMMPS setup)

# Apply force-controlled spherical indenter
# Target: 100 units compression in z-direction
fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -100.0

# Monitor forces and positions
variable fx equal f_1[1]
variable fy equal f_1[2] 
variable fz equal f_1[3]
variable ix equal f_1[4]
variable iy equal f_1[5]
variable iz equal f_1[6]

thermo_style custom step temp v_fx v_fy v_fz v_ix v_iy v_iz
thermo 100

run 10000
```

### Basic Planar Indentation

```lammps
# Planar indenter pushing from the top (z direction)
fix 2 all indent/force 15.0 plane z 20.0 hi 0.0 0.0 -50.0

# Monitor only the relevant force component
variable fz equal f_2[3]
variable pz equal f_2[6]  # plane position

thermo_style custom step temp v_fz v_pz
run 5000
```

## Advanced Usage Examples

### 1. Ramped Force Loading

Apply gradually increasing force to avoid system shock:

```lammps
# Define ramping function
variable max_force equal -200.0
variable ramp_steps equal 5000
variable current_force equal "min(v_max_force*step/v_ramp_steps, v_max_force)"

fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 v_current_force

# Output force evolution
fix output all print 100 "${step} ${current_force}" file force_ramp.dat
```

### 2. Cyclic Loading

Apply sinusoidal force variation:

```lammps
variable frequency equal 0.001  # cycles per timestep
variable amplitude equal 50.0
variable baseline equal -100.0
variable cyclic_force equal "v_baseline + v_amplitude*sin(2*PI*v_frequency*step)"

fix 1 all indent/force 8.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 v_cyclic_force

# Track force and position evolution
variable fx equal f_1[1]
variable fy equal f_1[2]
variable fz equal f_1[3]
variable iz equal f_1[6]

fix output all print 50 "${step} ${cyclic_force} ${fz} ${iz}" file cyclic_test.dat
```

### 3. Multi-Directional Loading

Apply forces in multiple directions simultaneously:

```lammps
# Define complex force patterns
variable fx_target equal "10.0*sin(step*dt*0.1)"
variable fy_target equal "5.0*cos(step*dt*0.05)"
variable fz_target equal "-50.0 - 20.0*step/10000.0"

fix 1 all indent/force 12.0 sphere 5.0 5.0 7.0 1.0 v_fx_target v_fy_target v_fz_target

# Monitor all components
variable fx equal f_1[1]
variable fy equal f_1[2]
variable fz equal f_1[3]

thermo_style custom step v_fx_target v_fx v_fy_target v_fy v_fz_target v_fz
```

### 4. Moving Indenter with Force Control

Combine position and force control:

```lammps
# Moving indenter center
variable center_x equal "5.0 + 2.0*step/10000.0"
variable center_y equal "5.0"
variable center_z equal "7.0"

# Constant downward force
fix 1 all indent/force 10.0 sphere v_center_x v_center_y v_center_z 1.0 0.0 0.0 -75.0

# Track indenter position evolution
variable ix equal f_1[4]
variable iy equal f_1[5]
variable iz equal f_1[6]

fix output all print 100 "${step} ${center_x} ${ix} ${iy} ${iz}" file moving_indenter.dat
```

## Monitoring and Diagnostics

### Essential Monitoring Variables

```lammps
# Force components (actual forces on indenter)
variable fx equal f_1[1]  # x-direction force
variable fy equal f_1[2]  # y-direction force  
variable fz equal f_1[3]  # z-direction force

# Indenter position (adjusted automatically)
variable ix equal f_1[4]  # x position
variable iy equal f_1[5]  # y position
variable iz equal f_1[6]  # z position (or plane position for planar indenters)

# System diagnostics
compute mobile_temp mobile temp
compute mobile_press all pressure mobile_temp

thermo_style custom step c_mobile_temp c_mobile_press v_fx v_fy v_fz v_ix v_iy v_iz
```

### Convergence Monitoring

```lammps
# Target forces (for comparison)
variable target_fx equal 0.0
variable target_fy equal 0.0
variable target_fz equal -100.0

# Force errors
variable error_fx equal "abs(v_fx - v_target_fx)"
variable error_fy equal "abs(v_fy - v_target_fy)"
variable error_fz equal "abs(v_fz - v_target_fz)"

# Relative errors (percentage)
variable rel_error_fx equal "v_error_fx/abs(v_target_fx + 1e-10)*100"
variable rel_error_fy equal "v_error_fy/abs(v_target_fy + 1e-10)*100"
variable rel_error_fz equal "v_error_fz/abs(v_target_fz + 1e-10)*100"

fix output all print 100 "${step} ${error_fx} ${error_fy} ${error_fz} ${rel_error_fx} ${rel_error_fy} ${rel_error_fz}" file convergence.dat
```

## Best Practices

### 1. System Preparation

Always equilibrate your system before applying force control:

```lammps
# Step 1: Energy minimization
minimize 1.0e-6 1.0e-8 1000 10000

# Step 2: Thermal equilibration
velocity all create 300.0 12345
fix nvt all nvt temp 300.0 300.0 100.0
run 20000

# Step 3: Switch to force control
unfix nvt
fix nve all nve
fix indent all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -50.0
```

### 2. Force Magnitude Guidelines

Start with conservative force values:

- **Small systems (< 1000 atoms)**: Forces of 1-10 units
- **Medium systems (1000-10000 atoms)**: Forces of 10-100 units  
- **Large systems (> 10000 atoms)**: Forces of 100-1000 units

Scale forces based on your system's characteristic energy scale.

### 3. Timestep Considerations

Use smaller timesteps when force control is active:

```lammps
# Standard timestep for equilibration
timestep 0.001

# Reduced timestep for force control
timestep 0.0005  # or smaller if needed

fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -100.0
```

### 4. Temperature Control

Monitor system temperature carefully:

```lammps
# Add gentle thermostat to prevent overheating
compute mobile_temp mobile temp
fix thermostat mobile langevin 300.0 300.0 100.0 48279

# Or use velocity rescaling
fix temp_rescale all temp/rescale 100 300.0 300.0 0.02 1.0
```

## Troubleshooting

### Issue: Forces Not Converging

**Symptoms:** Large force errors, oscillating values
**Solutions:**
1. Reduce target forces
2. Use smaller timesteps
3. Ensure proper system equilibration
4. Check for system instabilities (temperature spikes)

```lammps
# Example: More conservative parameters
timestep 0.0001
variable small_force equal -10.0  # Start smaller
fix 1 all indent/force 5.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 v_small_force
```

### Issue: System Overheating

**Symptoms:** Temperature spikes, atoms flying apart
**Solutions:**
1. Add appropriate thermostat
2. Reduce force magnitude
3. Use gradual force ramping

```lammps
# Add temperature control
compute mobile_temp mobile temp
fix thermostat mobile nvt temp 300.0 300.0 100.0

# Monitor temperature
thermo_modify temp mobile_temp
```

### Issue: No Initial Contact

**Symptoms:** Forces remain at zero despite non-zero targets
**Solutions:**
1. Move indenter closer to material
2. Increase indenter size
3. Let the algorithm find contact automatically (it should)

```lammps
# Manual positioning closer to material
fix 1 all indent/force 10.0 sphere 5.0 5.0 6.0 1.0 0.0 0.0 -50.0
#                                            ^-- moved closer (was 7.0)
```

### Issue: Performance Problems

**Symptoms:** Very slow simulation
**Solutions:**
1. The algorithm runs 500 iterations per timestep maximum
2. Well-converged forces typically need 50-200 iterations
3. Poor convergence will always use the full 500 iterations

Monitor convergence to identify performance issues:

```lammps
# Add custom output to track performance
fix perf all print 100 "${step} Force_errors: ${error_fx} ${error_fy} ${error_fz}" screen yes
```

## Force Sign Convention Reference

Understanding force signs is crucial:

```lammps
# COMPRESSION (typical indentation)
# Negative z-force = indenter pushes DOWN into material
fix 1 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -100.0

# TENSION (pulling)  
# Positive z-force = indenter pulls UP on material
fix 2 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 100.0

# MIXED LOADING
# Complex stress states
fix 3 all indent/force 10.0 sphere 5.0 5.0 7.0 1.0 50.0 25.0 -100.0
#                                                    ^    ^     ^
#                                                    |    |     compression in z
#                                                    |    tension in y  
#                                                    tension in x
```

## Integration with Other LAMMPS Features

### With Compute Commands

```lammps
# Analyze stress during indentation
compute stress all stress/atom NULL
compute virial all reduce sum c_stress[1] c_stress[2] c_stress[3]

thermo_style custom step v_fx v_fy v_fz c_virial[1] c_virial[2] c_virial[3]
```

### With Dump Commands

```lammps
# Output atomic trajectories and forces
dump trajectory all atom 100 indent_trajectory.lammpstrj
dump forces all custom 100 forces.dump id x y z fx fy fz

# Custom dump with indenter info
variable step_num equal step
fix custom_output all print 100 "${step_num} ${fx} ${fy} ${fz} ${ix} ${iy} ${iz}" file indenter_data.txt
```

### With Regions and Groups

```lammps
# Apply force control only to specific regions
region contact sphere 5.0 5.0 7.0 2.0
group contact_atoms region contact

fix 1 contact_atoms indent/force 10.0 sphere 5.0 5.0 7.0 1.0 0.0 0.0 -50.0
```

## Example: Complete Indentation Simulation

```lammps
# Complete example: Nanoindentation simulation
dimension 3
boundary p p p
atom_style atomic

# Create substrate
lattice fcc 3.5
region substrate block 0 20 0 20 0 15
create_box 2 substrate
create_atoms 1 substrate

# Set masses and potentials
mass 1 1.0
mass 2 1.0
pair_style lj/cut 8.0
pair_coeff * * 1.0 3.5 8.0

# Create fixed bottom layer
region bottom block INF INF INF INF 0 2
group bottom region bottom
group mobile subtract all bottom
set group bottom type 2

# Initial conditions
velocity mobile create 300.0 12345
fix 1 all nve
fix 2 bottom setforce 0.0 0.0 0.0
fix 3 mobile langevin 300.0 300.0 100.0 54321

# Equilibration
timestep 0.001
thermo 1000
run 20000

# Switch to force-controlled indentation
unfix 3
fix 3 mobile langevin 300.0 300.0 100.0 54321
timestep 0.0005

# Gradual force increase
variable max_force equal -200.0
variable ramp_time equal 20000
variable current_force equal "min(v_max_force*step/v_ramp_time, v_max_force)"

fix indent all indent/force 50.0 sphere 10.0 10.0 18.0 2.0 0.0 0.0 v_current_force

# Monitor everything
variable fx equal f_indent[1]
variable fy equal f_indent[2]
variable fz equal f_indent[3]
variable ix equal f_indent[4]
variable iy equal f_indent[5]
variable iz equal f_indent[6]
compute mobile_temp mobile temp

thermo_style custom step c_mobile_temp v_current_force v_fx v_fy v_fz v_ix v_iy v_iz
thermo 500

# Output data
fix output all print 100 "${step} ${current_force} ${fx} ${fy} ${fz} ${ix} ${iy} ${iz}" file nanoindentation.dat
dump trajectory all atom 1000 trajectory.lammpstrj

# Run the simulation
run 50000

print "Nanoindentation simulation completed!"
print "Final force: ${fz} (target was ${current_force})"
print "Final indenter position: ${iz}"
```

This guide should help you get started with force-controlled indentation and troubleshoot common issues. Remember that force control is inherently more complex than position control, so patience and careful monitoring are essential for successful simulations.