# CoRaLS

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://rprechelt.github.io/CoRaLS.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://rprechelt.github.io/CoRaLS.jl/dev)
[![Build Status](https://github.com/rprechelt/CoRaLS.jl/workflows/CI/badge.svg)](https://github.com/rprechelt/CoRaLS.jl/actions) -->

## Installation

Install the latest version of Julia from [julialang.org/downloads](https://julialang.org/downloads/) and [add Julia to your system path](https://julialang.org/downloads/platform/#optional_add_julia_to_path). 

Also ensure python and the `matplotlib` package are installed with [Anaconda](https://www.anaconda.com/)/[Miniconda](https://conda.io/miniconda.html)/[Miniforge](https://github.com/conda-forge/miniforge) (e.g. so that it is visible with `conda list matplotlib`).

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

From now on, simply start julia with the environment activated with `--project=/path/to/root/CoRaLS/directory`. If you are in the folder simply run:

```bash
julia --project=.
```

Then you can import CoRaLS and use it as above:

```julia
julia> using CoRaLS
```

## Testing

To test the installation, run the following from the CoRaLS directory:

```bash
julia --project=. test/runtests.jl
```

## Developers

To build documentation, run:

```bash
julia --project=docs/ docs/make.jl
```

To preview docs locally, run, then go to the url in a browser:

```bash
julia --project=docs/ -e 'using LiveServer; serve(dir="docs/build")'
```

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
