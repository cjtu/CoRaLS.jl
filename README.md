# CoRaLS

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rprechelt.github.io/CoRaLS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rprechelt.github.io/CoRaLS.jl/dev)
[![Build Status](https://github.com/rprechelt/CoRaLS.jl/workflows/CI/badge.svg)](https://github.com/rprechelt/CoRaLS.jl/actions) -->


## Installation

Install Julia v1.8. Since this is an older release, either go to [old releases](https://julialang.org/downloads/oldreleases/) and [add Julia to your dydtem path](https://julialang.org/downloads/platform/#optional_add_julia_to_path) or use an install helper like [juliaup](https://github.com/JuliaLang/juliaup) or [jill](https://github.com/abelsiqueira/jill).

Also ensure python and the `matplotlib` package are installed with [Anaconda](https://www.anaconda.com/)/[Miniconda](https://conda.io/miniconda.html)/[Miniforge](https://github.com/conda-forge/miniforge) (e.g. via `conda install matplotlib`).

Clone this repository:

```sh
git clone git@github.com:cjtu/CoRaLS.jl.git
```

Change directory to the cloned repository. If using multiple Julia versions on this machine, you may need to specify `juliaup default 1.8.5` or `jill switch 1.8`. Then start Julia via:

```sh
cd CoRaLS.jl
julia
```

In the Julia REPL, enter package mode by pressing `]` (note the `>pkg`). Building the PyCall package will give Julia access to the Python installation:

```julia
pkg> build PyCall
```

Then activate the CoRaLS environment and instantiate it to install the required packages and compile CoRaLS.jl:

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

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
