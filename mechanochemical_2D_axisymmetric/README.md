# mechanochemical_2D_axisymmetric

This repository contains a Julia implementation for mechanochemical surface flow simulations in axisymmetric geometries. It is derived from the `SurfaceBulkViscousFlows` by Eric Neiva codebase and refocused on surface-only dynamics for mechanochemical models, removing the bulk flow contribution and increasing resolution on the cell surface.

## Overview

`mechanochemical_2D_axisymmetric` provides tools to simulate coupled surface shape evolution and chemical/mechanical fields on an axisymmetric surface. The code uses the Gridap ecosystem for finite element assembly and supports numerical examples in the `examples` folder.

Key points:
- Surface-only viscous flow and mechanochemical coupling
- Axisymmetric geometry suited for 2D reduction of rotationally symmetric problems
- Demonstrator examples in `examples/SurfaceViscousFlows`
- Based on the original `SurfaceBulkViscousFlows` repository ( https://github.com/ericneiva/SurfaceBulkViscousFlows ) , adapted for mechanochemical focus

## Requirements

- Julia 1.11 or later
- Recommended packages are listed in `Project.toml`

Required dependencies include:
- `Gridap`
- `GridapEmbedded`
- `GridapPETSc`
- `Plots`
- `PythonPlot`
- `Symbolics`
- `DelimitedFiles`
- `SparseMatricesCSR`

## Installation

1. Clone this repository:

```bash
git clone <repository-url> mechanochemical_2D_axisymmetric
cd mechanochemical_2D_axisymmetric
```

2. Start Julia with the project environment:

```bash
julia --project=.
```

3. Instantiate dependencies once:

```julia
julia> import Pkg; Pkg.instantiate()
```

## Usage

Open a Julia session from the repository root and run the relevant example script.

Example:

```julia
julia> include("examples/SurfaceViscousFlows/SurfaceViscousFlows.jl")
```

Output is written to `output/` and may include VTK files for visualization.

## Repository structure

- `src/`
  - Core implementation files such as `MechanochemicalAxisymmetricVector.jl`, `SurfaceBulkViscousFlows.jl`, `WeakForms.jl`, `SolverFunctions.jl`
  - Helper modules for activity, density initialization, tangent operators, and plotting
- `examples/`
  - `SurfaceViscousFlows/` contains the current mechanochemical and surface flow demonstrator example
  - `Verification/` contains verification examples for baseline tests
- `output/`
  - Generated simulation results and VTK files
- `Project.toml`, `Manifest.toml`
  - Julia package environment and dependency manifest

## Examples

Current example directories:
- `examples/SurfaceViscousFlows/` 

Use `include(...)` on the desired example file from within a Julia session.

## Development notes

This repository is intended for research and demonstration of mechanochemical surface flow models in axisymmetric domains. It is not packaged as a registered Julia package, but it can be run directly from the repository root using Julia's project environment.

## Contact

For questions about this repository, contact the authors or maintainers listed in `Project.toml` and the source files.

## License

This project inherits the license found in `LICENSE`.
