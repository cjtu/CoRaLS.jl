using Documenter
using CoRaLS

makedocs(
    sitename="CoRaLS Documentation",
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Getting Started" => "getting_started.md",
        "Tutorials" => "tutorials.md",
        "The CoRaLS.jl Overview" => "corals_jl_overview.md",
        "acceptance.jl: The Main Monte Carlo Loop" => "acceptance_jl.md",
        "simulate.jl: The Main Cosmic Ray Simulation Script" => "simulate_jl.md",
        "Modules Used by acceptance.jl" => "modules_used_by_acceptance.md",
        "All Remaining Modules" => "all_remaining_modules.md",
        "Introduction to Julia" => "intro_to_julia.md",
        "Additional Details" => "additional_details.md"
    ],
    modules=[CoRaLS],
    warnonly=true
)

deploydocs(
    repo="github.com/yourusername/CoRaLS.jl.git",
    branch="jr-docstrings",
    push_preview=true,
)
