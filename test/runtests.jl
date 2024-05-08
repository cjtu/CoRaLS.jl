using Test
# using Random

# Set a random seed for reproducible tests
# SEED = 10  # trunc(Int, time())
# println("# Random seed ($SEED), on VERSION == $VERSION")
# Random.seed!(SEED)

@testset "All Tests" begin
    include("ut_geometry.jl");
    include("ut_fresnel.jl");
    include("ut_ice.jl");    
end
