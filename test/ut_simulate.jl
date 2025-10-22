using Test
using CoRaLS
using LinearAlgebra
using StaticArrays
using Unitful
using Unitful: m, km, MHz

@testset "Simulate Helper Functions" begin
    
    @testset "refract_angle" begin
        # Test air to glass transition (n=1 to n=1.5)
        θ_in = π/4
        θ_out = CoRaLS.refract_angle(θ_in, 1.0, 1.5)
        @test θ_out ≈ asin(sin(π/4) / 1.5)
        
        # Test critical angle behavior
        n1, n2 = 1.5, 1.0
        θ_critical = asin(n2/n1)
        θ_out_critical = CoRaLS.refract_angle(θ_critical, n1, n2)
        @test θ_out_critical ≈ π/2
        
        # Test normal incidence (0 degrees)
        @test CoRaLS.refract_angle(0.0, 1.0, 1.5) ≈ 0.0
        
        # Test symmetry: if we refract in then refract out, we get back original angle
        θ_in = π/6
        θ_intermediate = CoRaLS.refract_angle(θ_in, 1.0, 1.5)
        θ_back = CoRaLS.refract_angle(θ_intermediate, 1.5, 1.0)
        @test θ_back ≈ θ_in
    end
    
    @testset "calculate_surface_geometry" begin
        # Create simple test geometry
        origin = SA[0.0km, 0.0km, 1737.4km]  # Point on lunar surface
        antenna = SA[0.0km, 0.0km, 1757.4km]  # 20km above origin
        slopemodel = CoRaLS.NoSlope()
        
        normal, surface_normal, obs, view = CoRaLS.calculate_surface_geometry(origin, antenna, slopemodel)
        
        # Normal should be radially outward and normalized
        @test norm(normal) ≈ 1.0
        @test normal ≈ origin / norm(origin)
        
        # With NoSlope, surface_normal should equal normal
        @test surface_normal ≈ normal
        
        # Observation vector should point from surface to antenna
        @test obs ≈ antenna - origin
        
        # View should be normalized obs
        @test norm(view) ≈ 1.0
        @test view ≈ obs / norm(obs)
    end
    
    @testset "calculate_direct_distances" begin
        depth = 10.0m
        θ_reg = π/4  # 45 degrees
        obs = SA[0.0km, 0.0km, 20.0km]
        
        Drego, Dvacuum, x = CoRaLS.calculate_direct_distances(depth, θ_reg, obs)
        
        # Check horizontal distance
        @test x ≈ depth * tan(θ_reg)
        
        # Check regolith distance (Pythagorean theorem)
        @test Drego ≈ sqrt(x^2 + depth^2)
        
        # Check vacuum distance
        @test Dvacuum ≈ norm(obs)
        
        # Test with vertical ray (θ_reg = 0)
        Drego_vert, _, x_vert = CoRaLS.calculate_direct_distances(depth, 0.0, obs)
        @test x_vert ≈ 0.0m
        @test Drego_vert ≈ depth
    end
    
    @testset "calculate_emission_vector" begin
        # Setup basic geometry
        normal = SA[0.0, 0.0, 1.0]
        proj = SA[1.0, 0.0, 0.0]
        depth = 5.0m
        θ_i = π/6
        Nsurf = 1.8
        NXmax = 1.7
        Rmoon = 1737.4km
        
        # Test direct emission
        emit_direct, θ_emit_direct, invariant = CoRaLS.calculate_emission_vector(
            normal, proj, depth, θ_i, Nsurf, NXmax, Rmoon, false)
        
        @test norm(emit_direct) ≈ 1.0
        @test invariant ≈ Nsurf * Rmoon * sin(θ_i)
        
        # Test reflected emission (should have opposite vertical component)
        emit_refl, θ_emit_refl, _ = CoRaLS.calculate_emission_vector(
            normal, proj, depth, θ_i, Nsurf, NXmax, Rmoon, true)
        
        @test norm(emit_refl) ≈ 1.0
        # Vertical component should be opposite sign
        @test emit_direct[3] ≈ -emit_refl[3]
        # Horizontal components should be the same
        @test emit_direct[1] ≈ emit_refl[1]
        @test emit_direct[2] ≈ emit_refl[2]
    end
    
    @testset "calculate_incident_polarization" begin
        # Test with simple geometry
        emit = SA[0.0, 0.0, 1.0]  # Vertical emission
        axis = SA[0.0, 1.0, 0.0]  # Horizontal axis
        
        pol = CoRaLS.calculate_incident_polarization(emit, axis)
        
        # Polarization should be normalized
        @test norm(pol) ≈ 1.0
        
        # Polarization should be perpendicular to emission direction
        @test abs(pol ⋅ emit) < 1e-10
        
        # Test that pol is in the plane perpendicular to emit
        # pol = -emit × (emit × axis) should be parallel to axis projection onto plane perp to emit
        axis_perp_component = axis - (axis ⋅ emit) * emit
        if norm(axis_perp_component) > 1e-10
            axis_perp_normalized = axis_perp_component / norm(axis_perp_component)
            @test abs(abs(pol ⋅ axis_perp_normalized) - 1.0) < 1e-10
        end
    end
    
    @testset "construct_incident_coordinate_system" begin
        normal = SA[0.0, 0.0, 1.0]
        proj = SA[1.0, 0.0, 0.0]
        θ_i = π/4
        
        incident, niperp, nipar = CoRaLS.construct_incident_coordinate_system(normal, proj, θ_i)
        
        # All vectors should be normalized
        @test norm(incident) ≈ 1.0
        @test norm(niperp) ≈ 1.0
        @test norm(nipar) ≈ 1.0
        
        # niperp and nipar should be perpendicular
        @test abs(niperp ⋅ nipar) < 1e-10
        
        # niperp should be perpendicular to normal
        @test abs(niperp ⋅ normal) < 1e-10
        
        # nipar should be perpendicular to incident
        @test abs(nipar ⋅ incident) < 1e-10
        
        # Form right-handed coordinate system
        @test nipar ≈ niperp × incident
    end
    
    @testset "construct_transmission_coordinate_system" begin
        view = SA[0.0, 0.0, 1.0]
        normal = SA[0.0, 1.0, 0.0]
        
        ntperp, ntpar = CoRaLS.construct_transmission_coordinate_system(view, normal)
        
        # Both should be normalized
        @test norm(ntperp) ≈ 1.0
        @test norm(ntpar) ≈ 1.0
        
        # Should be perpendicular to each other
        @test abs(ntperp ⋅ ntpar) < 1e-10
        
        # ntperp should be perpendicular to normal
        @test abs(ntperp ⋅ normal) < 1e-10
        
        # ntpar should be perpendicular to view
        @test abs(ntpar ⋅ view) < 1e-10
        
        # Form right-handed system
        @test ntpar ≈ ntperp × view
    end
    
    @testset "calculate_transmitted_polarization" begin
        # Setup coordinate systems
        pol = SA[1.0, 0.0, 0.0]
        niperp = SA[1.0, 0.0, 0.0]
        nipar = SA[0.0, 1.0, 0.0]
        ntperp = SA[1.0, 0.0, 0.0]
        ntpar = SA[0.0, 1.0, 0.0]
        tperp = 0.8
        tpar = 0.6
        
        poltr = CoRaLS.calculate_transmitted_polarization(pol, niperp, nipar, ntperp, ntpar, tperp, tpar)
        
        # Should be in the ntperp-ntpar plane
        z_component = poltr[3]
        @test abs(z_component) < 1e-10
        
        # Test that coefficients are applied correctly
        @test poltr ≈ tperp * (pol ⋅ niperp) * ntperp + tpar * (pol ⋅ nipar) * ntpar
    end
    
    @testset "calculate_polarization_angle" begin
        # Test pure perpendicular polarization
        ntperp = SA[1.0, 0.0, 0.0]
        ntpar = SA[0.0, 1.0, 0.0]
        poltr_perp = SA[1.0, 0.0, 0.0]
        
        θpol_perp = CoRaLS.calculate_polarization_angle(poltr_perp, ntpar, ntperp)
        @test θpol_perp ≈ 0.0 atol=1e-10
        
        # Test pure parallel polarization
        poltr_par = SA[0.0, 1.0, 0.0]
        θpol_par = CoRaLS.calculate_polarization_angle(poltr_par, ntpar, ntperp)
        @test θpol_par ≈ π/2 atol=1e-10
        
        # Test 45 degree polarization
        poltr_45 = (ntperp + ntpar) / sqrt(2)
        θpol_45 = CoRaLS.calculate_polarization_angle(poltr_45, ntpar, ntperp)
        @test θpol_45 ≈ π/4 atol=1e-10
    end
    
    @testset "calculate_reflected_distances" begin
        depth = 5.0m
        ice_depth = 10.0m
        θ_reg = π/6
        obs = SA[0.0km, 0.0km, 20.0km]
        origin = SA[0.0km, 0.0km, 1737.4km]
        normal = SA[0.0, 0.0, 1.0]
        proj = SA[1.0, 0.0, 0.0]
        
        Drego, Dvacuum, x1, x2, source = CoRaLS.calculate_reflected_distances(
            depth, ice_depth, θ_reg, obs, origin, normal, proj)
        
        # Check x1 and x2 calculations
        @test x1 ≈ (ice_depth - depth) * tan(θ_reg)
        @test x2 ≈ ice_depth * tan(θ_reg)
        
        # Check total regolith distance
        expected_Drego = sqrt(x1^2 + (ice_depth - depth)^2) + sqrt(x2^2 + ice_depth^2)
        @test Drego ≈ expected_Drego
        
        # Check vacuum distance
        @test Dvacuum ≈ norm(obs)
        
        # Check source calculation
        expected_source = origin - (2.0 * ice_depth - depth) * normal - (x1 + x2) * proj
        @test source ≈ expected_source
    end
    
    @testset "calculate_observation_angles" begin
        origin = SA[0.0km, 0.0km, 1737.4km]
        axis = SA[0.0, 1.0, 0.0]
        normal = origin / norm(origin)
        view = SA[0.0, 0.0, 1.0]
        antenna = SA[0.0km, 0.0km, 1757.4km]
        
        zenith, θ, ϕ, el = CoRaLS.calculate_observation_angles(origin, axis, normal, view, antenna)
        
        # Zenith should be in valid range [0, π]
        @test 0.0 <= zenith <= π
        
        # θ should be in valid range [0, π]
        @test 0.0 <= θ <= π
        
        # ϕ should be in valid range [-π, π] or [0, 2π]
        @test -π <= ϕ <= 2π
        
        # Elevation angle should be in reasonable range
        @test -π/2 <= el <= π/2
    end
    
end
