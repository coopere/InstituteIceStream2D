[![DOI](https://zenodo.org/badge/1121945435.svg)](https://doi.org/10.5281/zenodo.18063098)

# Institute Ice Stream 2D Model

This repository contains MATLAB code for simulating two-dimensional cross-sectional ice stream dynamics using convex optimization. The model solves the ice flow problem using a variational formulation that minimizes total power dissipation subject to physical constraints.

## Overview

The code implements a finite element method for solving ice flow equations with:
- **Glen's flow law** (n=3, p=4/3) for ice rheology
- **Spatially variable basal friction** (rock/sediment transitions, channels)
- **Temperature-dependent viscosity** (thermal-mechanical coupling in dimensional version)
- **Unstructured triangular mesh** discretization

## Files

### Executable Scripts

- **`obstacle2D_cvx_dimensional.m`** - Dimensional version with physical units (SI)
  - Includes thermal-mechanical coupling
  - Temperature-dependent ice rheology
  - Frictional heating
  - Comparison with observational data
  
- **`obstacle2D_cvx_nondimensional.m`** - Nondimensional version
  - Useful for understanding scaling relationships
  - Simplified for testing and validation
  - Faster execution

### Data Files

- **`vel_profile_full.mat`** - Observed velocity profiles for comparison
- **`vel_profiles.mat`** - Additional velocity profile data
- **`iceColorMap.mat`** - Colormap for temperature visualization (optional - code will use default colormap if not available)

## Dependencies

### Required MATLAB Toolboxes

1. **CVX Optimization Toolbox** (http://cvxr.com/cvx/)
   - Required for solving the convex optimization problem
   - Free academic license available
   - Install and run `cvx_setup` before use

2. **distmesh2d** (http://persson.berkeley.edu/distmesh/)
   - Unstructured triangular mesh generation
   - Must be in MATLAB path

## Usage

### Running the Dimensional Model

```matlab
% Make sure CVX and distmesh2d are in your path
% Run the script
obstacle2D_cvx_dimensional
```

The script will:
1. Generate an unstructured mesh
2. Initialize temperature field
3. Iteratively solve for velocity and temperature (thermal-mechanical coupling)
4. Display visualization plots

### Running the Nondimensional Model

```matlab
% Make sure CVX and distmesh2d are in your path
% Run the script
obstacle2D_cvx_nondimensional
```

The script will:
1. Generate an unstructured mesh
2. Solve for velocity field
3. Display visualization plots

## Model Description

### Governing Equations

The model solves a variational problem minimizing total power dissipation:

```
J = ∫ (1/p) * a * |∇u|^p dΩ - ∫ f * u dΩ + ∫_Γ τ_c * |u| dΓ
```

where:
- `u` is the velocity field
- `a` is the temperature-dependent viscosity coefficient
- `f` is the driving stress (ρg sin(α))
- `τ_c` is the basal friction law
- `p = 4/3` (Glen's law with n=3)

### Basal Friction Law

The basal friction `τ_c(x,y,u)` is piecewise defined:

1. **Left boundary** (x < sed_trans_l): Very high friction (no slip)
2. **Rock region** (sed_trans_l < x < rock_trans): Power law friction
   ```
   τ = (β_var*(x/L) + β)^2 * |u|^(5/2)
   ```
3. **Sediment region** (rock_trans < x < sed_trans_r): Linear friction with channel
   ```
   τ = (sed_var*(1-x/L) + ch_str*exp(-|x-ch_loc|/ch_decay) + sed_const) * |u|
   ```
4. **Right boundary** (x > sed_trans_r): Very high friction (no slip)

### Temperature-Dependent Rheology

The viscosity coefficient `a` depends on temperature via Arrhenius relation:

```
a = (1/2) * [A * E * exp(-Q_h/R * (1/T_h - 1/T_star))]^(-1/3)
```

where:
- `A` is the preexponential constant
- `E` is the enhancement factor
- `Q_h` is the activation energy (temperature-dependent)
- `T_h` is the pressure-adjusted temperature
- `R` is the gas constant

### Thermal-Mechanical Coupling

The dimensional model includes:
1. **Frictional heating**: `f_therm = 2 * τ_E * ε_E`
2. **Heat conduction**: `∇·(k∇T) = f_therm`
3. **Iterative coupling**: Velocity and temperature fields are solved iteratively until convergence

## Model Parameters

### Physical Constants (Dimensional Version)

- Ice density: `ρ = 900 kg/m³`
- Gravitational acceleration: `g = 9.8 m/s²`
- Surface slope: `α = 2.4×10⁻³ rad`
- Flow law exponent: `p = 4/3` (Glen's law n=3)
- Preexponential constant: `A = 3.5×10⁻²⁵ 1/(s·Pa³)`
- Gas constant: `R = 8.314 J/(K·mol)`
- Activation threshold: `T_star = 263.15 K` (-10°C)

## Output

The scripts generate several visualization plots:

1. **Velocity field**: 2D contour plot of velocity magnitude
2. **Surface velocity profile**: Comparison with observations
3. **Velocity gradient**: Comparison with observations
4. **Temperature field** (dimensional only): 2D temperature distribution
5. **Basal friction parameterization**: Plot of friction coefficients

## Citation 

For more details on this model and application to Institute Ice Stream, please refer to the accompanying paper:
```
Suckale, Jenny and Elsworth, Cooper W. An antiplane strain model for evaluating shear-margin stability (Ortholine v1.0), under review, 2026
```

## Contact

Please contact the authors with any inquiries on this code at cooper.elsworth@gmail.com

## Acknowledgments

- CVX optimization toolbox developers
- distmesh2d developers
- Institute Ice Stream field observations
