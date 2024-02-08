# Getting Started with CoRaLS

After successfully installing CoRaLS, this guide will help you get started with using the CoRaLS environment for your simulations.

## Using the CoRaLS Environment

Once you have completed the initial build of PyCall and imported the CoRaLS code, subsequent uses of CoRaLS will be more straightforward and quicker.

### Activating the CoRaLS Environment

To begin, activate the CoRaLS environment each time you start a new session:

1. Enter the package mode in Julia's REPL by pressing `]`.
2. Activate the CoRaLS environment:
   ```julia
   pkg> activate .
   pkg> instantiate
    ```
    Exit package mode by pressing `backspace`

### Importing CoRaLS
    ```julia
    julia> using CoRaLS
    ```

    Once you have executed this command, CoRaLS is loaded into your Julia session, and you are ready to run functions and scripts specific to CoRaLS.

### Running a Test
    To verify that CoRaLS is working as expected, you can run a simple test:

    ```
      julia> CoRaLS.plot_acceptance(CoRaLS.acceptance(10000, 20))
    ```

    You can also run some of the test scripts in the `test` directory to verify that CoRaLS is working as expected. First, navigate to the `test` directory and run the test scripts:

    ```julia
    julia> include("test/test_acceptance.jl")
    julia> include("test/test_simulation.jl")
    ```
    If the tests pass, you can be confident that CoRaLS is working as expected.
    
