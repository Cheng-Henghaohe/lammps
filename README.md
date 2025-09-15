# Force-Controlled Indentation for LAMMPS

This directory contains the implementation of a force-controlled indentation fix for LAMMPS, based on the original position-controlled `fix indent` command.

## Overview

The `fix indent/force` command provides load-controlled indentation instead of position-controlled indentation. Instead of specifying a fixed indenter position, you specify target forces (fx, fy, fz), and the indenter position is automatically adjusted each timestep to achieve those target forces.

## Key Features

- **Force Control**: Specify target forces instead of positions
- **Automatic Convergence**: Internal iteration loop adjusts indenter position until forces converge
- **Multi-directional**: Control forces independently in x, y, and z directions
- **Initial Contact Handling**: Automatically moves indenter to establish contact when forces are initially zero
- **Flexible Tolerance**: Configurable force tolerance and iteration limits

## Files

- `fix_indent_force.h` - Header file
- `fix_indent_force.cpp` - Implementation
- `in.indent_force_test` - Test input script
- `README.md` - This documentation

## Syntax

```
fix ID group-ID indent/force K geometry fx fy fz [keyword value ...]
```

- `K` = force constant (same as original indent)
- `geometry` = sphere, plane, etc. (same geometries as original indent)
- `fx, fy, fz` = target forces in each direction
- Optional keywords: `side`, `units` (same as original indent)

## Examples

### Spherical Indenter with Compression Force
```lammps
# Apply -100 force in z direction (compression)
fix 1 all indent/force 10.0 sphere 20.0 10.0 15.0 5.0 0.0 0.0 -100.0
```

### Planar Indenter
```lammps  
# Apply -50 force with planar indenter
fix 2 all indent/force 10.0 plane z 10.0 lo 0.0 0.0 -50.0
```

### Using Variables for Dynamic Forces
```lammps
variable force_z equal -100.0*step/1000.0  # Gradually increase force
fix 3 all indent/force 10.0 sphere 20.0 10.0 15.0 5.0 0.0 0.0 v_force_z
```

## Algorithm Details

1. **Force Control Loop**: Each timestep contains an internal iteration loop
2. **Force Calculation**: Calculate current indenter forces using standard indent physics
3. **Error Calculation**: Compare actual forces with target forces  
4. **Position Adjustment**: Adjust indenter position based on force errors
5. **Convergence Check**: Continue until forces are within tolerance or max iterations reached

## Control Parameters

- `tolerance = 0.01` (1% of target force)
- `max_iterations = 100` (per timestep)  
- `adjustment_factor = 0.01` (position adjustment step size)

These are currently hard-coded but could be made configurable in future versions.

## Force Sign Convention

- **Positive forces**: Compression (indenter pushing into material)
- **Negative forces**: Tension (indenter pulling material)
- This follows the physics: indenter[1-3] represents reaction forces on the indenter

## Compilation

To compile this fix into LAMMPS:

1. Copy `fix_indent_force.h` and `fix_indent_force.cpp` to your LAMMPS `src/` directory
2. Recompile LAMMPS
3. The fix will be available as `indent/force`

Alternatively, you can compile as a separate package by including these files in a new package directory.

## Limitations

- Currently only implements sphere and plane geometries
- Control parameters are hard-coded (not user-configurable)
- Cylinder and cone geometries need implementation
- No validation for unrealistic force targets

## Future Improvements

- Add user-configurable control parameters
- Implement remaining geometries (cylinder, cone)
- Add convergence diagnostics and warnings
- Implement adaptive adjustment factors
- Add force ramping capabilities

## Testing

Use the provided `in.indent_force_test` script to test the functionality:

```bash
lmp < in.indent_force_test
```

The script tests both spherical and planar indenters with different force targets.