import Logging; Logging.disable_logging(Logging.Info)

@testset verbose = true "acceptance.jl" begin
    @testset "Save and load events" begin
        A = acceptance(10, 5; region=create_region("whole_moon"), savefile="ut_acceptance.jld2")
        B = load_acceptance("ut_acceptance.jld2")
        @test A.ntrials == B.ntrials
        rm("ut_acceptance.jld2")
    end
end
