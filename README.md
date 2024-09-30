# Mechanochemical polarization

`mechanochemical_polarization` is distributed as a supplemental material to the paper:

> Henry De Belly, Andreu Fernández Gallén, Evelyn Strickland, Dorothy C. Estrada, Patrick J. Zager, Janis K Burkhardt, Hervé Turlier,  Orion D. Weiner Title: Long range mutual activation establishes Rho and Rac polarity during cell migration 

The code is based on the [Gridap](https://github.com/gridap/Gridap.jl) for the finite-element implementation.

## Theory and implementation

See the [Theory and implementation](theory_implementation.md) document.

## Requirements

* **`Gridap`** see [installation instructions here](https://github.com/gridap/Gridap.jl?tab=readme-ov-file#readme).

Gridap is a registered package in the official Julia package registry. Thus, the installation of Gridap is straight forward using the Julia's package manager. Open the Julia REPL, type ] to enter package mode, and install as follows

>pkg> add Gridap

* **`Plots`** Julia library for Plotting


## Usage

On the main section in the code there is specified the different constants and length of the simulation. There the values can be changed to the desired ones and then one can run the code as a normal Julia code.

Each example is in the folder `./Examples`, including:
* `./Examples/Heatmap`
* `./Examples/Opto_activation` 
* `./Examples/Opto_activation_local_inhibition` 

To run one of the examples, run it like
```
cd ./Examples
julia
include("./mechanochemical_opto_front2back.jl")

```

Results will be stored in the relevant folder, each example code has a folder name in the "simulation" parameter. The folder name will be for example Heatmap and the subfolders that will create to store the simulation data will depend on the simulation parameters introduced.

The code stores a copy of itself in the output folder, to be able to check after a given time for the exact version of the code and of mechanical parameters that where used.
