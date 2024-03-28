using Test
using LinearAlgebra
import Unitful: km, m, s, H, eV, K, MHz, V

# include("../src/constants.jl")
# include("../src/geometry.jl")
# include("../src/CoRaLS.jl")
using CoRaLS

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