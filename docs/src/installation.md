# Installation

This guide will walk you through the process of installing and setting up the CoRaLS simulation package.

## Quick Install Steps
1. ** Install Julia v1.8.**: 
    Since this is an older release, either go to [old releases](https://julialang.org/downloads/oldreleases/) and [add Julia to your system path](https://julialang.org/downloads/platform/#optional_add_julia_to_path) or use an install helper like [juliaup](https://github.com/JuliaLang/juliaup) or [jill](https://github.com/abelsiqueira/jill).

    You will also ensure Python and `matplotlib` are installed.

2. **Clone the CoRaLS repository**:
    ```sh
    git clone git@github.com:cjtu/CoRaLS.jl.git
    ```

    
3. ** Install CoraLS.jl**:
    Change directory to the cloned repository. If using multiple Julia versions on this machine, you may need to specify `juliaup default 1.8.5` or `jill switch 1.8`. Then start Julia via:

    ```sh
    cd CoRaLS.jl
    julia
    ```

    In the Julia REPL (the interoreter that starts when you type julie, stands for read-eval-print loop), enter package mode by pressing `]` (note `julia>` changes to `>pkg`). Adding and building the PyCall package will give Julia access to the Python installation specified at the path supplied in quotes to `ENV["PYTHON"]` (leave blank `""` to use system default python):

    ```julia
    julia> ENV["PYTHON"]=""
    julia> ]
    pkg> add PyCall
    pkg> build PyCall
    ```

    Then activate the CoRaLS environment and instantiate it to install the required packages and compile CoRaLS.jl:

    ```julia
    pkg> activate .
    pkg> instantiate
    ```

    Now CoRaLS is ready to use. To exit package mode, press `backspace` (prompt is now `julia>`). Try importing CoRaLS and running a test:

4. **Verify Installation**:
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
## Troubleshooting
If you encounter any issues during the installation process, please refer to the FAQs or contact the simulations team.