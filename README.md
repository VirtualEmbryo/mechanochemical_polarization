# Mechanochemical Polarization

`mechanochemical_polarization` is a collection of Julia codes for simulating mechanochemical models of cell polarization and migration. This repository is distributed as supplemental material to the paper:

> Henry De Belly*, Andreu F. Gallen, Evelyn Strickland, Dorothy C. Estrada, David Sanchez Godinez, Eric Neiva, Patrick J. Zager, Tamas L. Nagy, Janis K Burkhardt, Hervé Turlier*, Orion D. Weiner*. "Long range mutual activation establishes Rho and Rac polarity during cell migration" **Nature Cell Biology** 2026.

The simulations model the coupled dynamics of Rho GTPase (Rho and Rac) activity and mechanical forces on the cell surface, leading to polarization during migration.

---

## Key Features

- Finite-element based simulations using [Gridap.jl](https://github.com/gridap/Gridap.jl)
- Mechanochemical coupling between biochemical signaling and mechanical deformation
- Support for various scenarios: optogenetic activation, local inhibition, heatmaps of parameter spaces
- Axisymmetric 2D simulations for efficient computation
- Visualization tools for Rho/Rac distributions and surface evolution
- 2D axisymmetric implementation based on the original `SurfaceBulkViscousFlows` repository ( https://github.com/ericneiva/SurfaceBulkViscousFlows ) , adapted for mechanochemical focus

## Repository Structure

- `Mechanochemical_general_code.jl`: Main simulation script for 1D mechanochemical models
- `Plots_RhoRacA.jl`: Plotting utilities for Rho and Rac distributions
- `Examples/`: Collection of example simulations
  - `Heatmap/`: Parameter sweep simulations
  - `Opto_activation/`: Optogenetic perturbation studies
  - `Opto_activation_local_inhibition/`: Combined opto and inhibition
- `mechanochemical_2D_axisymmetric/`: Julia package for 2D axisymmetric surface flow simulations
  - `src/`: Source code modules
  - `examples/`: Example scripts
  - `output/`: Simulation outputs


## Requirements

- **Julia**: Version 1.6 or later (tested with 1.11)
- **Gridap.jl**: Finite element library
- **Plots.jl**: For visualization
- Other dependencies listed in `mechanochemical_2D_axisymmetric/Project.toml`

## Installation

1. Clone this repository:
   ```bash
   git clone <repository-url>
   cd mechanochemical_polarization
   ```

2. Install Julia dependencies:
   ```bash
   julia --project=mechanochemical_2D_axisymmetric
   ```
   ```julia
   julia> import Pkg; Pkg.instantiate()
   ```

3. For the main scripts, ensure Gridap and Plots are installed:
   ```julia
   julia> using Pkg; Pkg.add(["Gridap", "Plots"])
   ```

---

## Usage

### Running Main Simulations

The main simulation script is `Mechanochemical_general_code.jl`. Modify the parameters in the main section and run:

```bash
julia Mechanochemical_general_code.jl
```

Key parameters include:
- `L`: Domain length
- `T`: Simulation time
- `Δt`: Time step
- Mechanical parameters: `vCTE`, `α`, `β`, `σₐ₀`
- Biochemical parameters: reaction rates, diffusion coefficients

### Running Examples

Navigate to the `Examples/` directory and run specific example scripts:

```bash
cd Examples
julia --project=../mechanochemical_2D_axisymmetric
julia> include("mechanochemical_opto_front2back.jl")
```

Available examples:
- `Heatmap_Time2LosePol.jl`: Parameter sweeps for polarization loss
- `localinhibition_opto_front2back.jl`: Local inhibition with opto perturbations
- `mechanochemical_opto_front2back.jl`: Basic mechanochemical with opto

Results are saved in subfolders based on simulation parameters.

### 2D Axisymmetric Simulations

For full 2D surface simulations, use the `mechanochemical_2D_axisymmetric` package:

```bash
cd mechanochemical_2D_axisymmetric
julia --project=.
julia> include("examples/SurfaceViscousFlows/SurfaceViscousFlows.jl")
```

## Theory and Implementation

The model couples Rho/Rac GTPase signaling with surface mechanics. Rho promotes contractility, while Rac drives protrusion. Mechanical feedback creates positive feedback loops leading to polarization.

For detailed mathematical formulation, see the supplemental materials of the associated paper: .

## Output and Visualization

- Simulation data saved as VTU files for ParaView visualization
- PNG plots of concentration profiles over time
- Text files with parameter values and timing data

Use `Plots_RhoRacA.jl` for custom plotting.

---

## Authors and Contributing

This code is provided as supplemental material. For questions or modifications, contact the authors.

- Andreu Fernández Gallén
- [Eric Neiva](https://github.com/ericneiva)
- [Hervé Turlier](https://github.com/VirtualEmbryo)

This work re-used parts of original code from Eric Neiva: [SurfaceBulkViscousFlows](https://github.com/VirtualEmbryo/SurfaceBulkViscousFlows), based on a work by Eric Neiva and Hervé Turlier: [arXiv:2505.05723](https://arxiv.org/abs/2505.05723)

## Funding

This project was funded by the European Union's Horizon 2020 research and innovation programme under the **Marie Skłodowska-Curie** grant agreement No. 101150259 (A.F.G)  and received funding under the **Marie Skłodowska-Curie** grant agreement No. 101105565 (E.N) and the **European Research Council** grant agreement No. 949267 (H.T.).

## Citing

If you use this code in your research, please cite:

> Henry De Belly*, Andreu F. Gallen, Evelyn Strickland, Dorothy C. Estrada, David Sanchez Godinez, Eric Neiva, Patrick J. Zager, Tamas L. Nagy, Janis K Burkhardt, Hervé Turlier*, Orion D. Weiner*. "Long range mutual activation establishes Rho and Rac polarity during cell migration" **Nature Cell Biology** 2026.

<!-- Repository DOI:
> vertAX contributors (2019). vertAX: a differentiable vertex model framework. [doi:10.5281/zenodo.3555620](https://zenodo.org/record/3555620) -->

## License

This code is released under the [Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)](https://creativecommons.org/licenses/by-sa/4.0/) [`License`](LICENSE).

