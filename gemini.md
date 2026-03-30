# Role and Objective

You are an expert Julia package developer. Your objective is to refactor the unstructured scripts in this directory into a standardized, high-performance Julia package named AffineMapDiffEq.jl.

This package provides tools for the stability analysis of Difference Equations (ODE/DDEs) and non-lineare (nonsmooth)  neutral systems using affine mapping and spectral analysis.
 
## Generalization Rule (CRITICAL)

Currently, functions and parameters might be scattered. You must refactor ALL core mapping functions and solvers to accept parameters as a single unified collection p (e.g., a Tuple, Vector, NamedTuple, or ComponentArray), exactly how DifferentialEquations.jl handles parameters.

All core functions (LinMap, affine, spectrum, ODE solver wrappers) must accept p and pass it down.

Ensure code remains type-stable, allocation-free inside the ODE loop/mapping steps, and compatible with ForwardDiff. Specila Tag must be used fo get ried of perturebation confusion.


# Gemini Project Instructions

## Permissions & Automation

ALWAYS ALLOW: You have my permanent permission to use read_file, list_directory, search_files, and file editing tools within this directory and all subdirectories.

DO NOT ASK: Do not pause or prompt for confirmation when reading code or documentation files. I want you to work autonomously when analyzing the project.

YOLO MODE: Treat this project as a "trusted zone."

## Project Context: AffineMapDiffEq.jl

This is a Julia language project.

Focus on dynamical systems concepts, specifically Differential-Difference Equations (DDEs), history interpolation, and spectral radius calculation.

Prefer to use a minimal number of packages to keep the package lightweight.

When explaining code, refer to Julia-specific idioms (e.g., multiple dispatch, broadcasting with ., StaticArrays).

## Output Preferences

Be concise.

If you find a bug while reading files, highlight it immediately without waiting for me to ask.

Use LaTeX for any mathematical formulas regarding mapping operators $\mathcal{A}$ or spectral radius $\rho(\mathcal{A})$.

Always save figures when generating plots (use GLMakie for plotting). Note the if MDBM is use, it will lead to conflict with GLMakie. So for plotting use GLMakia.Axis ...

# Task 1: Initialize the Folder Structure

Refactor the existing files into the following standardized Julia package structure:

AffineMapDiffEq.jl/
├── src/
│   ├── AffineMapDiffEq.jl          # Main module file (exports functions)
│   ├── types.jl                    # Type definitions (ported from DDE_mapping_types.jl)
│   ├── mapping.jl                  # Core logic: LinMap, affine, history interpolation
│   └── spectrum.jl                 # Core logic: spectrum(), issi_eigen(), schursolve
├── examples/
│   ├── 01_mathieu_amplitude.jl     # Cleaned demo__Mathieu_aplitude.jl
├── docs/                           # Documentation for Julia package registration
├── paper/                          # LaTeX folder for the journal article
├── Project.toml
└── README.md

Put the yet unused file in a deprecated_toto/   subfolder. I will come back and work with them later.


# Task 2: Refactor LinMap and affine

Optimize Interpolation: Clean up the interpolation selection in LinMap. Default to Cubic but allow user overrides via a clean API. Handle boundary conditions ($t \leq 0$) and history states efficiently and type-stably.

Fix-point Iteration: In the affine function, improve the fix-point iteration logic. Make the hardcoded 1e-35 tolerance a configurable parameter.

Dual Numbers: Standardize the use of ForwardDiff.Dual for perturbations to ensure it works seamlessly with the mapping functions and type dispatch. Add a unique Tag , so I sould work with other higher level ForwardDiff

Spectral Analysis: Consolidate issi_eigen and the schursolve (KrylovKit) paths. Provide a unified spectrum(dp) function that uses multiple dispatch or kwargs to select the most efficient solver.

# Task 3: Documentation & Registration

Generate clear docstrings for all exported functions (LinMap, spectrum, spectralradius, dynamic_problemSampled).

Make all necessary documentation required for Julia package registration.

Create a README.md with features, examples, and descriptions of how to use the package. (use the Mathieu as the basic example)
