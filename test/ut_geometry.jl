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

    @testset "Region Sampling" begin
        # Circle AOI
        circle = Circle(0.0, 0.0, 100km)
        for _ in 1:100
            p = random_point_in_aoi(circle)
            @test is_in_aoi(p, circle)
            @test norm(p) ≈ Rmoon atol = 1e-5km
        end

        # SphericalCap AOI
        cap = SphericalCap(0.0, 0.0, 10.0)
        for _ in 1:100
            p = random_point_in_aoi(cap)
            @test is_in_aoi(p, cap)
            @test norm(p) ≈ Rmoon atol = 1e-5km
        end

        # Quadrangle AOI
        quad = Quadrangle(-10.0, 10.0, -10.0, 10.0)
        for _ in 1:100
            p = random_point_in_aoi(quad)
            @test is_in_aoi(p, quad)
            @test norm(p) ≈ Rmoon atol = 1e-5km
        end

        # Pole AOIs
        south = SouthPoleAOI(-80.0)
        north = NorthPoleAOI(80.0)
        for _ in 1:100
            ps = random_point_in_aoi(south)
            pn = random_point_in_aoi(north)
            @test is_in_aoi(ps, south)
            @test is_in_aoi(pn, north)
            @test norm(ps) ≈ Rmoon atol = 1e-5km
            @test norm(pn) ≈ Rmoon atol = 1e-5km
        end
    end

    @testset "Region Creation and Containment" begin
        # WholeMoonRegion accepts all points
        p1 = latlon_to_cartesian(0.0, 0.0, Rmoon)
        p2 = latlon_to_cartesian(-180.0, 90.0, Rmoon)
        p3 = latlon_to_cartesian(180.0, -90.0, Rmoon)
        @test is_in_region(p1, WholeMoonRegion())
        @test is_in_region(p2, WholeMoonRegion())
        @test is_in_region(p3, WholeMoonRegion())

        # Region with Circle AOI
        region = create_region("circle:20,60,100,1.0")
        p_in = latlon_to_cartesian(20.5, 60.5, Rmoon)
        p_out = latlon_to_cartesian(0.0, 0.0, Rmoon)
        @test is_in_region(p_in, region)
        @test !is_in_region(p_out, region)

        # Region with SphericalCap AOI
        region = create_region("cap:0,0,10,1.0")
        p_in = latlon_to_cartesian(5.0, 0.0, Rmoon)
        p_out = latlon_to_cartesian(15.0, 0.0, Rmoon)
        @test is_in_region(p_in, region)
        @test !is_in_region(p_out, region)

        # PolarRegion
        region = create_region("polar:south,-80,1.0")
        p_in = latlon_to_cartesian(-85.0, 0.0, Rmoon)
        p_out = latlon_to_cartesian(-75.0, 0.0, Rmoon)
        @test is_in_region(p_in, region)
        @test !is_in_region(p_out, region)

        # Test AOI fraction
        region = create_region("polar:south,-80,0.0")
        @test !is_in_region(p_in, region)
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