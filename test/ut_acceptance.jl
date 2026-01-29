import Logging; Logging.disable_logging(Logging.Info)
using Test
using CoRaLS
using CoRaLS: ustrip, EeV
using DelimitedFiles

# Helper functions for plaintext I/O
function save_acceptance_txt(A::CoRaLS.Acceptance, filename::String, header::String="")
    # Save energies, rAΩ, dAΩ, and differential spectra to a simple text file
    # energies has nbins+1 elements (bin edges), rAΩ and dAΩ have nbins elements
    # Format: labeled lines with space-separated values, all rounded to 4 decimal places
    r_diff = differential_spectrum(A.energies, A.rAΩ, 1yr)
    d_diff = differential_spectrum(A.energies, A.dAΩ, 1yr)
    
    open(filename, "w") do io
        if !isempty(header)
            # Write header with # prefix for each line
            for line in split(header, '\n')
                println(io, "# ", line)
            end
        end
        println(io, "Energy_bins: ", join(round.(ustrip.(A.energies), digits=4), ' '))
        println(io, "Acceptance_reflected: ", join(round.(ustrip.(A.rAΩ), digits=4), ' '))
        println(io, "Counts_reflected_1yr: ", join(round.(r_diff, digits=4), ' '))
        println(io, "Acceptance_direct: ", join(round.(ustrip.(A.dAΩ), digits=4), ' '))
        println(io, "Counts_direct_1yr: ", join(round.(d_diff, digits=4), ' '))
    end
end

function load_acceptance_txt(filename::String)
    # Load energies, rAΩ, dAΩ, and differential spectra from text file
    # Returns a named tuple with energies, rAΩ, r_diff, dAΩ, d_diff
    # Skip lines starting with #
    lines = filter(line -> !startswith(line, '#'), readlines(filename))
    
    # Parse each line after the label
    energies = parse.(Float64, split(split(lines[1], ": ")[2]))
    rAΩ = parse.(Float64, split(split(lines[2], ": ")[2]))
    r_diff = parse.(Float64, split(split(lines[3], ": ")[2]))
    dAΩ = parse.(Float64, split(split(lines[4], ": ")[2]))
    d_diff = parse.(Float64, split(split(lines[5], ": ")[2]))
    
    return (energies=energies, rAΩ=rAΩ, r_diff=r_diff, dAΩ=dAΩ, d_diff=d_diff)
end

@testset verbose = true "acceptance.jl" begin
    @testset "Check acceptance is consistent after changes" begin
        import Random
        baseline_f = joinpath(@__DIR__, "ut_acceptance_test.txt")

        # Explicitly set parameters for current baseline test
        seed = 1254251
        altitude = 50km
        Random.seed!(seed)
        trigger = LPDA(Nant=8, Ntrig=3, θ0=-48.0, altitude=altitude)
        acc_params = (
            min_energy = 6EeV,
            max_energy = 600.0EeV,
            ν_min = 150MHz,
            ν_max = 1000MHz,
            dν = 30MHz,
            geometrymodel = ScalarGeometry(),
            indexmodel = StrangwayIndex(),
            fieldmodel = ARW(),
            divergencemodel = MixedFieldDivergence(),
            densitymodel = StrangwayDensity(),
            slopemodel = GaussianSlope(7.6),
            roughnessmodel = GaussianRoughness(2.0),
            iceroughness = GaussianIceRoughness(2.0cm),
            Nice = 1.305,
            Nbed = 1.6,
            ice_depth = 5.0m,
            ice_thickness = 1.0m,
            spacecraft = CircularOrbit(altitude),
            region = create_region("psr:south"),
            simple_area = true,
            trigger = trigger
        )
        # Make header of acc params to save to file
        header = """Random seed: $(seed)
            Trigger: LPDA(Nant=8, Ntrig=3, θ0=-48.0, altitude=$(altitude))
            Acceptance parameters (200000 events, 7 bins):"""
        for (key, val) in pairs(acc_params)
            if key != :trigger
                header *= "\n  - $(key): $(val)"
            end
        end

        # Run the acceptance calculation
        Anew = acceptance(200000, 5; acc_params...)
        rAΩnew = round.(ustrip.(Anew.rAΩ), digits=4)
        rAnew_diff = round.(differential_spectrum(Anew.energies, Anew.rAΩ, 1yr), digits=4)

        # Make baseline acceptance file if it doesn't exist
        if !isfile(baseline_f)
            @info "Baseline acceptance file $(baseline_f) not found; creating it."
            save_acceptance_txt(Anew, baseline_f, header)
        end

        # Load baseline acceptance from file
        Aold = load_acceptance_txt(baseline_f)
        rAΩold = Aold.rAΩ
        rAold_diff = Aold.r_diff
        println("(old) Reflected events in 1 yr = ", rAold_diff)
        println("(new) Reflected events in 1 yr = ", rAnew_diff)

        # Check if new acceptance (reflected) deviates from baseline by more than 1%
        # If so, print both, then fail. Delete and commit old file to accept the new baseline.
        if !all(isapprox.(rAΩold, rAΩnew; rtol=0.01))
            println("Acceptance changed more than 1%:")
            println("Aold.rAΩ = ", rAΩold)
            println("Anew.rAΩ = ", rAΩnew)
            println("To confirm, delete $(baseline_f), then rerun tests.")
        end
        @test all(isapprox.(rAΩold, rAΩnew; rtol=0.01))
    end
end
