using Test
using LinearAlgebra

include("../src/fresnel.jl")

@testset verbose = true "fresnel.jl" begin
    # Validate Fresnel reflectance / transmission 
    # note: most sources give the POWER reflection/transmission R and T.
    # we're concerned with the amplitude reflection/transmission r and t.
    # The conversion is below (used to check results from power calculations):
    T(t, θ, n1, n2) = t^2 * (n2 * cos(snell_θt(θ, n1, n2))) / (n1 * cos(θ))
    R(r) = r^2

    theta30 = π / 6
    theta45 = π / 4
    n_air = 1.0
    n_water = 1.33
    n_glass = 1.5

    @testset "Fresnell Air -> Water (30 deg)" begin
        θi = theta30
        n1 = n_air
        n2 = n_water
        rpar, rperp, tpar, tperp = fresnel_coeffs(θi, n1, n2)
        rpar2, rperp2, tpar2, tperp2 = fresnel_coeffs(θi, n1, n2; simple_t = false)

        @testset "Reflected" begin
            @test R(rpar) ≈ 0.0117 atol = 1e-3
            @test R(rperp) ≈ 0.0305 atol = 1e-3
        end

        @testset "Transmitted" begin
            @test T(tpar, θi, n1, n2) ≈ 0.988 atol = 1e-3
            @test T(tperp, θi, n1, n2) ≈ 0.969 atol = 1e-3
        end

        @testset "Conservation of Energy" begin
            @test R(rpar) + T(tpar, θi, n1, n2) ≈ 1.0 atol = 1e-5 
            @test R(rperp) + T(tperp, θi, n1, n2) ≈ 1.0 atol = 1e-5
        end

        @testset "Transmitted (full formulas)" begin
            @test rpar2 ≈ rpar
            @test rperp2 ≈ rperp
            @test tpar2 ≈ tpar
            @test tperp2 ≈ tperp
        end

        @testset "Convservation of Energy (full formulas)" begin
            @test R(rpar2) + T(tpar2, θi, n1, n2) ≈ 1.0 atol = 1e-5 
            @test R(rperp2) + T(tperp2, θi, n1, n2) ≈ 1.0 atol = 1e-5
        end
    end

    @testset "Fresnel Water -> Glass (45 deg)" begin
        θi = theta45
        n1 = n_water
        n2 = n_glass
        rpar, rperp, tpar, tperp = fresnel_coeffs(θi, n1, n2)
        rpar2, rperp2, tpar2, tperp2 = fresnel_coeffs(θi, n1, n2; simple_t = false)

        @testset "Reflected" begin
            @test R(rpar) ≈ 0.000137 atol = 1e-3
            @test R(rperp) ≈ 0.0117 atol = 1e-3
        end

        @testset "Transmitted" begin
            @test T(tpar, θi, n1, n2) ≈ 1.00 atol = 1e-3
            @test T(tperp, θi, n1, n2) ≈ 0.988 atol = 1e-3
        end

        @testset "Conservation of Energy" begin
            @test R(rpar) + T(tpar, θi, n1, n2) ≈ 1.0 atol = 1e-5 
            @test R(rperp) + T(tperp, θi, n1, n2) ≈ 1.0 atol = 1e-5
        end

        @testset "Transmitted (full formulas)" begin
            @test rpar2 ≈ rpar
            @test rperp2 ≈ rperp
            @test tpar2 ≈ tpar
            @test tperp2 ≈ tperp
        end

        @testset "Convservation of Energy (full formulas)" begin
            @test R(rpar2) + T(tpar2, θi, n1, n2) ≈ 1.0 atol = 1e-5 
            @test R(rperp2) + T(tperp2, θi, n1, n2) ≈ 1.0 atol = 1e-5
        end

    end

    @testset "Water -> Air (n2 > n1)" begin
        θi, n1, n2 = theta30, n_water, n_air
        rpar, rperp, tpar, tperp = fresnel_coeffs(θi, n1, n2)
        @testset "Reflected" begin
            @test R(rpar) ≈ 0.00469 atol = 1e-3
            @test R(rperp) ≈ 0.0455 atol = 1e-3
        end

        @testset "Transmitted" begin
            @test T(tpar, θi, n1, n2) ≈ 0.995 atol = 1e-3
            @test T(tperp, θi, n1, n2) ≈ 0.955 atol = 1e-3
        end
    end


    @testset "Normal Incidence" begin
        rpar, rperp, tpar, tperp = fresnel_coeffs(0, n_air, n_water)
        @test rpar == 0.0
        @test rperp == 0.0
        @test tpar == 1.0
        @test tperp == 1.0
    end

    @testset "Total Internal Reflection" begin
        # Choose angle just over the critical angle for TIR
        θi = fresnel_critical(n_air, n_glass) + 0.01  
        rpar, rperp, tpar, tperp = fresnel_coeffs(θi, n_air, n_water)
        @test rpar == 1.0
        @test rperp == 1.0
        @test tpar == 0.0
        @test tperp == 0.0
    end
end