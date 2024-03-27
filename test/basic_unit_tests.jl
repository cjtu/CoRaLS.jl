using Test
using LinearAlgebra
import Unitful: cm, km, m, s, H, eV, K, MHz, V

include("../src/constants.jl")
include("../src/geometry.jl")
include("../src/fresnel.jl")
include("../src/ice.jl")

# Test constants
@testset verbose = true "constants.jl" begin
    @test Rmoon ≈ 1736.0km
    @test c_0 ≈ 2.99792458e8m / s
    @test μ_0 ≈ 1.25663706212e-6H / m
    @test k_b ≈ 8.617333e-5eV / K
end

@testset verbose = true "geometry.jl" begin
    # Test random points on poles
    @testset "Random Points on Poles" begin
        point = random_north_pole_point()
        @test length(point) == 3
        @test norm(point) ≈ Rmoon atol = 1e-5km

        point = random_south_pole_point()
        @test length(point) == 3
        @test norm(point) ≈ Rmoon atol = 1e-5km
    end

    # Test random angles on cap
    @testset "Random Angles on Cap" begin
        theta_max = π / 4
        theta, phi = random_angles_on_cap(theta_max)
        @test 0 ≤ theta ≤ theta_max
        @test 0 ≤ phi ≤ 2π
    end

    # Test random point on cap
    @testset "Random Point on Cap" begin
        theta_max = π / 4
        point = random_point_on_cap(theta_max)
        @test length(point) == 3
        @test norm(point) ≈ Rmoon atol = 1e-5km
    end

    # Test spherical to cartesian conversion and vice versa
    @testset "Spherical and Cartesian Conversions" begin
        theta, phi, r = π / 4, π / 4, Rmoon
        cart = spherical_to_cartesian(theta, phi, r)
        sph = cartesian_to_spherical(cart[1], cart[2], cart[3])
        @test sph[1] ≈ theta atol = 1e-5
        @test sph[2] ≈ phi atol = 1e-5
        @test sph[3] ≈ r atol = 1e-5km
    end

    # Test horizon angle calculation
    @testset "Horizon Angle Calculation" begin
        altitude = 100km
        angle = horizon_angle(altitude)
        @test angle < 0 # Should be a negative angle
    end

    # Test random vector generation
    @testset "Random Vector Generation" begin
        vec = random_vector()
        @test length(vec) == 3
        @test norm(vec) ≈ 1.0 atol = 1e-5
    end

    # Test random direction generation
    @testset "Random Direction Generation" begin
        dir = random_direction()
        @test length(dir) == 3
        @test norm(dir) ≈ 1.0 atol = 1e-5
    end

    # Test spherical cap area calculation
    @testset "Spherical Cap Area Calculation" begin
        theta = π / 4
        area = spherical_cap_area(theta)
        expected_area = 2π * Rmoon^2 * (1.0 - cos(theta))
        @test area ≈ expected_area atol = 1e-5km^2
    end

    # Test intersection with sphere
    @testset "Intersection with Sphere" begin
        start = SVector{3}(100.0, 0.0, 0.0)
        direction = normalize(SVector{3}(-1.0, 0.0, 0.0))
        radius = 1
        intersection = intersect_with_sphere(start, direction, radius)
        @test intersection !== nothing
        @test norm(intersection) ≈ radius atol = 1e-5
        @test (norm(intersection) - radius) / radius < 1e-6
    end
end


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

@testset verbose = true "Ice Surface Roughness Models" begin
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
