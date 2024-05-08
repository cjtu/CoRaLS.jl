using Test
using LinearAlgebra

include("../src/fresnel.jl")

@testset verbose = true "fresnel.jl" begin
    theta0 = 0
    theta30 = π / 6
    theta45 = π / 4
    theta60 = π / 3
    n_air = 1.0
    n_water = 1.33
    n_glass = 1.5

    @testset "Parallel Refelection: fresnel_rpar" begin
        # Air to water transition at 45 degrees
        rpar = fresnel_rpar(theta45, n_air, n_water)
        @test 0.0 <= rpar <= 1.0

        rpar = fresnel_rpar(0.9826204, n_air, n_glass)
        @test 0.0 <= rpar <= 0.001
    end

    @testset "Pependicular Reflection: fresnel_rperp" begin
        # Air to water transition at 45 degrees
        rperp = fresnel_rperp(theta45, n_air, n_water)
        @test -1.0 <= rperp <= 0.0
    end

    # # Test for total internal reflection
    # @testset "Total Internal Reflection" begin
    #     ni = 1.5   # Refractive index of glass
    #     nt = 1.0   # Refractive index of air, to simulate the possibility of total internal reflection
    #     thetai = asin(nt / ni) + 0.01  # Angle just over the critical angle for total internal reflection
    #     rperp = fresnel_rperp(thetai, ni, nt)
    #     @test rperp == 1.0  # Expecting total internal reflection, so reflection coefficient should be 1
    # end

    @testset "Parallel Transmission: fresnel_tpar" begin
        # Scenario: Air to water transition at a 30-degree angle of incidence
        @testset "Air to Water at 30 degrees" begin
            tpar = fresnel_tpar(theta30, n_air, n_water)
            @test 0.0 <= tpar <= 1.0  # Transmission coefficient should be between 0 and 1
        end

        # Scenario: Checking conservation of energy
        @testset "Conservation of Energy" begin
            tpar = fresnel_tpar(theta45, n_air, n_glass)
            rpar = fresnel_rpar(theta45, n_air, n_glass)
            @test tpar + rpar ≈ 1.0 atol = 1e-5  # Sum of transmission and reflection coefficients should be 1
        end
    end

    @testset "Perpendicular Transmission: fresnel_tperp" begin
        @testset "Normal Incidence" begin
            tperp = fresnel_tperp(theta0, n_air, n_water)
            @test 0.0 <= tperp <= 1.0  # Expecting coefficient to be in [0,1]
        end

        @testset "Conservation of Energy" begin
            tperp = fresnel_tperp(theta45, n_air, n_glass)
            rperp = fresnel_rperp(theta45, n_air, n_glass)
            @test tperp + rperp ≈ 1.0 atol = 1e-5  # Ensuring energy conservation
        end
    end
end