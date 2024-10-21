using Test
using LinearAlgebra
import Unitful: km, m, s, H, eV, K, MHz, V

include("../src/geometry.jl")

@testset verbose = true "geometry.jl" begin

    @testset "Random Points on Poles" begin
        point = random_north_pole_point()
        @test length(point) == 3
        @test norm(point) ≈ Rmoon atol = 1e-5km
        @test cartesian_to_latlon(point)[1] > 80

        point = random_south_pole_point()
        @test length(point) == 3
        @test norm(point) ≈ Rmoon atol = 1e-5km
        @test cartesian_to_latlon(point)[1] < -80
    end

    @testset "Random Angles on Cap" begin
        theta_max = π / 4
        theta, phi = random_angles_on_cap(theta_max)
        @test 0 ≤ theta ≤ theta_max
        @test 0 ≤ phi ≤ 2π
    end

    @testset "Random Point on Cap" begin
        theta_max = π / 4
        point = random_point_on_cap(theta_max)
        @test length(point) == 3
        @test norm(point) ≈ Rmoon atol = 1e-5km
    end

    @testset "Spherical and Cartesian Conversions" begin
        theta, phi, r = π / 4, π / 4, Rmoon
        cart = spherical_to_cartesian(theta, phi, r)
        sph = cartesian_to_spherical(cart)
        @test sph[1] ≈ theta atol = 1e-5
        @test sph[2] ≈ phi atol = 1e-5
        @test sph[3] ≈ r atol = 1e-5km
    end

    # TODO: improve this test
    @testset "Horizon Angle Calculation" begin
        altitude = 100km
        angle = horizon_angle(altitude)
        @test angle < 0 # Should be a negative angle
    end

    @testset "Random Vector Generation" begin
        vec = random_vector()
        @test length(vec) == 3
        @test norm(vec) ≈ 1.0 atol = 1e-5
    end

    @testset "Random Direction Generation" begin
        dir = random_direction()
        @test length(dir) == 3
        @test norm(dir) ≈ 1.0 atol = 1e-5
    end

    @testset "Spherical Cap Area Calculation" begin
        a_moon = 4*pi*Rmoon^2

        # Null case
        @test spherical_cap_area(0) == 0km^2
        
        # Whole Moon
        @test spherical_cap_area(π) ≈ a_moon atol = 1e-5km^2

        # Half of the Moon
        @test spherical_cap_area(π/2) ≈ a_moon/2 atol = 1e-5km^2
    end

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
;