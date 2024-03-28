using Test
import Unitful: cm, km, m, s, H, eV, K, MHz, V

# include("../src/ice.jl")
using CoRaLS

@testset verbose = true "ice.jl" begin
    # Test Variables
    ν_test = 300MHz
    E_test = 1.0V / m
    θ_i_test = π / 4 # 45 degrees

    @testset "NoIceRoughness Model" begin
        # Expected behavior: No change to the electric field
        model = NoIceRoughness()
        ν_modified, E_modified = ice_roughness(model, [ν_test], [E_test], θ_i_test)
        @test ν_modified == [ν_test]
        @test E_modified == [E_test]
    end

    @testset "GaussianIceRoughness Model" begin
        # Expected behavior: Electric field is modified based on the roughness
        σ_test = 0.01m
        model = GaussianIceRoughness(σ_test)

        n = sqrt(3.0) # Refractive index for the test

        λ_test = c_0 / (ν_test * n) # Calculate wavelength for the test frequency
        g_test = (4π * (σ_test / λ_test) * cos(θ_i_test))^2.0 # Calculate g parameter

        ν_modified, E_modified = ice_roughness(model, [ν_test], [E_test], θ_i_test, n)

        # Verify the electric field is modified according to the Gaussian roughness model
        @test ν_modified == [ν_test]
        @test E_modified ≈ [E_test * exp(-g_test / 2.0)] atol = 1e-5V / m
    end
end
