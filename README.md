# CoRaLS

The `CoRaLS.jl` (Cosmic Ray Lunar Sounder) Monte Carlo model computes detection rates of Askaryan emission from cosmic ray interactions in lunar regolith.

## Installation

1. [Install Julia](https://julialang.org/install/) 
2. Install Python and the `matplotlib` package, the easiest way is with [Anaconda](https://www.anaconda.com/) (e.g. `conda list matplotlib`).
3. Fork or clone this repository:

```sh
git clone git@github.com:cjtu/CoRaLS.jl.git
```

4. (First time only): Setup path to python/matplotlib. In the Julia REPL, supply the path to your python environment from #2 in quotes to `ENV["PYTHON"]=""` (leave blank to use system default python). Then add and build the `PyCall` package:

```bash
$ julia

julia> ENV["PYTHON"]=""
julia> using Pkg
julia> Pkg.add("PyCall")
julia> Pkg.build("PyCall")
```

5. Run the test suite to test the installation:

```bash
julia --project=/path/to/CoRaLS.jl /path/to/CoRaLS.jl/test/runtests.jl
```

If the tests were successful, `CoRaLS.jl` is compiled and ready to use!

## Calculating rates with CoRaLS

1. Start Julia with the CoRaLS project active using the `--project` flag pointing to the CoRaLS directory (if you are in the directory use `--project=.`). For multithreaded mode, use `-t` to specify the number of threads (default is 1, "auto" chooses for you).

```bash
julia --project=/path/to/CoRaLS.jl -t "auto"
julia> using CoRaLS
```

2. import `CoRaLS` with:

```julia
julia> using CoRaLS
```

3. Compute and plot an acceptance:

```julia
julia> A = acceptance(10000, 20; region=create_region("psr:south"), spacecraft=CircularOrbit(50.0km))
julia> plot_acceptance(A)
```

See full documentation online at... (coming soon)

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

Paper coming soon.
