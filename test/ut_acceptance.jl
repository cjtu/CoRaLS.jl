import Logging; Logging.disable_logging(Logging.Info)
using Test
using CoRaLS
using CoRaLS: ustrip

@testset verbose = true "acceptance.jl" begin
    @testset "Check acceptance is consistent after changes" begin
        import Random
        baseline_f = joinpath(@__DIR__, "ut_acceptance_test.jld2")

        # Make baseline acceptance file if it doesn't exist
        if !isfile(baseline_f)
            @info "Baseline acceptance file $(baseline_f) not found; creating it."
            Random.seed!(1234)
            acceptance(100000, 5; region=WholeMoonRegion(), savefile=baseline_f)
        end
        # Load baseline acceptance from file
        Aold = load_acceptance(baseline_f)
        rAΩold = round.(ustrip.(Aold.rAΩ), digits=2)
        
        # Run current baseline with same random seed as above
        Random.seed!(1234)
        Anew = acceptance(100000, 5; region=WholeMoonRegion())
        rAΩnew = round.(ustrip.(Anew.rAΩ), digits=2)
        println("Aold.rAΩ = ", rAΩold)
        println("Anew.rAΩ = ", rAΩnew)
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
