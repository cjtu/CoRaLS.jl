# CoRaLS

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rprechelt.github.io/CoRaLS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rprechelt.github.io/CoRaLS.jl/dev)
[![Build Status](https://github.com/rprechelt/CoRaLS.jl/workflows/CI/badge.svg)](https://github.com/rprechelt/CoRaLS.jl/actions) -->

## Installation

Install the latest version of Julia from [julialang.org/downloads](https://julialang.org/downloads/) and [add Julia to your system path](https://julialang.org/downloads/platform/#optional_add_julia_to_path). 

Also ensure python and the `matplotlib` package are installed with [Anaconda](https://www.anaconda.com/)/[Miniconda](https://conda.io/miniconda.html)/[Miniforge](https://github.com/conda-forge/miniforge) (e.g. so that it is visible with `conda list matplotlib`).

Note: This code was originally written for v1.8 but most recently tested on v1.11. Julia takes care to not break backwards compatibility, so the code should run on newer v1.X.Y versions. To test on multiple versions, it is best to use a Julia install helper like [juliaup](https://github.com/JuliaLang/juliaup) or [jill](https://github.com/abelsiqueira/jill).

Clone this repository:

```sh
git clone git@github.com:cjtu/CoRaLS.jl.git
```

Change directory to the cloned repository and start Julia.

```sh
cd CoRaLS.jl
julia
```

In the Julia REPL, enter package mode by pressing `]` (notice `julia>` changes to `>pkg`). Adding and building the PyCall package will give Julia access to the Python installation specified at the path supplied in quotes to `ENV["PYTHON"]` (leave blank `""` to use system default python):

```julia
julia> ENV["PYTHON"]=""
julia> ]
pkg> add PyCall
pkg> build PyCall
```

Then activate the CoRaLS environment and instantiate it to install the required packages and compile CoRaLS.jl. You will need to activate the environment each time you start a new Julia session:

```julia
pkg> activate .
pkg> instantiate
```

Now CoRaLS is ready to use. To exit package mode, press `backspace` (prompt is now `julia>`). Try importing CoRaLS and running a test:

```julia
julia> using CoRaLS
julia> plot_acceptance(acceptance(10000, 20))
```

Or equivalently:

```julia
julia> import CoRaLS
julia> CoRaLS.plot_acceptance(CoRaLS.acceptance(10000, 20))
```

Exit julia with `ctrl+D` or `exit()`:

```julia
julia> exit()
```

## Using the CoRaLS environment

After the first build of PyCall and importing of CoRaLS code, subsequent runs will be faster and proceed as:

```julia
pkg> activate .
```

Backspace to exit `pkg>` mode, then:

```julia
julia> using CoRaLS
```

And run your functions/scripts.

## Description of model and included files

`src/` contains the CoRaLS.jl Monte Carlo model code

- `CoRaLS.jl`: Defines main Julia module / environment
- `acceptance.jl`: Main detection Monte Carlo loop, ntrials across nbins of energy `acceptance(ntrials, nbins; ...)`
  - returns `Acceptance` object containing # of successful events in each energy bin
  - computes polar PSR geometry and # of detectable events
  - samples the Auger spectrum for energies in input energy bins
  - passes Auger energies, altitude, other kwargs to main simulation, `throw_cosmicray()`
- `simulate.jl`: Main simulation of one cosmic ray event propagated to payload
  - `throw_cosmicray()`: Simulate single cosmic ray trail
    - returns structs containing E-field, polarization at payload from direct and reflected (off ice layer) signals
    - randomize detector location and cosmic ray angle
    - pass inputs to `compute_direct()` and `compute_reflected()`
  - `compute_direct()`: Compute RF solution direct from cosmic ray event with scalar geometry
    - returns `Direct` struct containing E-field, polarization at payload (and angles, path lengths)
    - compute path lengths in subsurface and vacuum accounting for refraction
    - pass inputs to `regolith_field()` to get E-field at payload
    - use `surface_transmission()` to get modified Fresnel coeffs correcting for surface roughness and distance to spacecraft (near/far-field correction, default: `MixedFieldDivergence`)
    - compute polarization vector of signal
  - `compute_reflected()`: Compute RF solution reflected off ice layer with scalar geometry (differeces from direct in **bold**)
    - returns `Reflected` struct containing E-field, polarization at payload (and angles, path lengths)
    - compute path lengths in subsurface and vacuum accounting for refraction, **accounting for pre- and post-reflection path length and checking for total internal reflection**
    - pass inputs to `regolith_field()` to get E-field at payload, **note: small extra distance before reflection is not accounted for**
    - **apply `ice_roughness()` to electric field**
    - use `surface_transmission()` to get modified Fresnel coeffs correcting for surface roughness and distance to spacecraft (near/far-field correction, default: `MixedFieldDivergence`)
    - compute polarization vector of signal **from top of ice layer, also compute polarization vector of signal reflected off the bottom of buried ice layer**
  - `efield.jl`: Compute integrated electric field at payload
    - `regolith_field()`: Load electric field at payload from ARW lookup table and attenuate through regolith. Returns frequency ν and electric field E
- Helper files called by `acceptance()` and subfunctions:
  - `constants.jl`: Physical constants
  - `cosmic_ray.jl`: Get cosmic ray grammage with `propagate_to_xmax()`
  - `fresnel.jl`: Compute Fresnel coefficients with & without near/far/mixed divergence models
  - `geometry.jl`: Various geometry helper functions
  - `ice.jl`: Correct for `ice_roughness()` when reflecting off ice layer
  - `moon.jl`: Compute PSR area and whether a point hits a PSR
  - `plots.jl`: Plotting functions
  - `regolith.jl`: Get `regolith_density()` for attenuation calculation
  - `spectrum.jl`: Load Auger spectrum with/without random sampling
  - `surface.jl`: Correct for `surface_transmission()`, with / without roughness
  - `trigger.jl`: Defines the trigger conditions for positive detections at the payload
- Other files **NOT** used in `acceptance()`:
  - `detector.jl` Antenna simulation `create_antenna()`
  - `raytrace.jl`: Raytrace from surface to ice layer and get angles assuming specular reflection
  - `simulate_vector.jl`: Alternate direct/reflected simulation using vector geometry
  - `sky.jl`: Compute `sky_temperature()` and `(extra)galactic_noise()` for antenna simulation
  - `toy.jl`: "Toy" models with different MC sampling strategies
  - `utils.jl`: Helper functions `polarization_angle()` used in `plots.jl` and `retrigger()` used in tests

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
