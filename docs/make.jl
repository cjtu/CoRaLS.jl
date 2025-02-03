using Documenter
using CoRaLS

makedocs(
    sitename="CoRaLS Documentation",
    format=Documenter.HTML(),
    pages=[
        "Getting Started" => "index.md",
        "acceptance.jl: Main Monte Carlo Loop" => "acceptance_jl.md",
        "simulate.jl: Main Cosmic Ray Simulation" => "simulate_jl.md",
        "Submodules of acceptance.jl" => "submodules_acceptance.md",
        "Other Modules" => "other_modules.md",
        "Tutorials" => "tutorials.md",
    ],
    modules=[CoRaLS],
    warnonly=true
)

deploydocs(
    repo="github.com/cjtu/CoRaLS.jl.git",
    push_preview=true,
)
