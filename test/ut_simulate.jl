using Test
using CoRaLS
using LinearAlgebra
using StaticArrays
using Unitful
using Unitful: m, km, MHz, EeV

@testset "Simulate Helper Functions" begin
    
    @testset "init_event_geometry" begin
        origin = SA[0.0km, 0.0km, 1737.4km]
        antenna = SA[0.0km, 0.0km, 1757.4km]
        Xmax = SA[0.0km, 0.0km, 1737.395km]  # 5m deep
        indexmodel = CoRaLS.ConstantIndex()
        slopemodel = CoRaLS.NoSlope()
        
        normal, obs, view, depth, Nsurf, NXmax, proj = 
            CoRaLS.init_event_geometry(origin, antenna, Xmax, indexmodel, slopemodel)
        
        # Test depth calculation
        @test depth ≈ 5.0m
        
        # Test normalization of vectors
        @test norm(normal) ≈ 1.0
        @test norm(view) ≈ 1.0
        @test norm(proj) ≈ 1.0
        
        # Test observation vector calculation
        @test obs ≈ antenna - origin
        
        # Test projection is perpendicular to origin normal
        origin_normal = origin / norm(origin)
        @test abs(proj ⋅ origin_normal) < 1e-10
        
        # Test refractive indices are reasonable (>1 for regolith)
        @test Nsurf > 1.0
        @test NXmax > 1.0
        
        # Test with sloped surface - normal should differ from geometric normal
        slopemodel_with_slope = CoRaLS.GaussianSlope(10.0)
        normal_sloped, _, _, _, _, _, _ = 
            CoRaLS.init_event_geometry(origin, antenna, Xmax, indexmodel, slopemodel_with_slope)
        # With slope, normal may differ from origin/norm(origin)
        @test norm(normal_sloped) ≈ 1.0  # Still normalized
    end
    
    @testset "calculate_direct_distances" begin
        depth = 10.0m
        θ_reg = π/4
        obs = SA[0.0km, 0.0km, 20.0km]
        
        Drego, Dvacuum, x = CoRaLS.calculate_direct_distances(depth, θ_reg, obs)
        
        # Test Pythagorean theorem
        @test Drego^2 ≈ x^2 + depth^2
        @test Dvacuum ≈ 20.0km
        
        # Test vertical ray (θ_reg = 0) - critical edge case
        Drego_vert, _, x_vert = CoRaLS.calculate_direct_distances(depth, 0.0, obs)
        @test x_vert ≈ 0.0m
        @test Drego_vert ≈ depth
        
        # Test shallow angle (θ_reg near π/2) - should have large x
        Drego_shallow, _, x_shallow = CoRaLS.calculate_direct_distances(depth, π/2 - 0.1, obs)
        @test x_shallow > depth * 5  # Should be much larger than depth
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
        
        # Test path geometry
        @test x1 ≈ (ice_depth - depth) * tan(θ_reg)
        @test x2 ≈ ice_depth * tan(θ_reg)
        
        # Test total regolith distance is sum of two segments
        d1 = sqrt(x1^2 + (ice_depth - depth)^2)
        d2 = sqrt(x2^2 + ice_depth^2)
        @test Drego ≈ d1 + d2
        
        # Test edge case: depth == ice_depth (xmax at ice boundary)
        Drego_edge, _, x1_edge, x2_edge, _ = CoRaLS.calculate_reflected_distances(
            ice_depth, ice_depth, θ_reg, obs, origin, normal, proj)
        @test x1_edge ≈ 0.0m  # No distance before ice
    end
    
    @testset "calculate_emission_vector" begin
        normal = SA[0.0, 0.0, 1.0]
        proj = SA[1.0, 0.0, 0.0]
        depth = 5.0m
        θ_i = π/6
        Nsurf = 1.8
        NXmax = 1.7
        Rmoon = 1737.4km
        
        # Test direct vs reflected have opposite vertical components
        emit_direct, θ_emit_direct = CoRaLS.calculate_emission_vector(
            normal, proj, depth, θ_i, Nsurf, NXmax, Rmoon, false)
        emit_refl, θ_emit_refl = CoRaLS.calculate_emission_vector(
            normal, proj, depth, θ_i, Nsurf, NXmax, Rmoon, true)
        
        @test norm(emit_direct) ≈ 1.0
        @test norm(emit_refl) ≈ 1.0
        @test emit_direct[3] ≈ -emit_refl[3]  # Opposite vertical components
        @test emit_direct[1] ≈ emit_refl[1]   # Same horizontal
        
        # Test Snell's law is applied (angles should match)
        @test θ_emit_direct ≈ θ_emit_refl
    end
    
    @testset "calculate_incident_polarization" begin
        # Test perpendicularity (main requirement)
        emit = normalize(SA[1.0, 1.0, 1.0])
        axis = SA[0.0, 1.0, 0.0]
        
        pol = CoRaLS.calculate_incident_polarization(emit, axis)
        
        @test norm(pol) ≈ 1.0
        @test abs(pol ⋅ emit) < 1e-10  # Perpendicular to emission
        
        # Test degenerate case: emit parallel to axis
        emit_parallel = SA[0.0, 1.0, 0.0]
        pol_parallel = CoRaLS.calculate_incident_polarization(emit_parallel, axis)
        @test norm(pol_parallel) < 1e-10 || norm(pol_parallel) ≈ 1.0  # Either zero or normalized
    end
    
    @testset "construct_incident_coordinate_system" begin
        normal = SA[0.0, 0.0, 1.0]
        proj = SA[1.0, 0.0, 0.0]
        θ_i = π/4
        
        incident, niperp, nipar = CoRaLS.construct_incident_coordinate_system(normal, proj, θ_i)
        
        # Test orthogonality (critical property)
        @test abs(niperp ⋅ nipar) < 1e-10
        @test abs(niperp ⋅ incident) < 1e-10
        @test abs(nipar ⋅ incident) < 1e-10
        
        # Test right-handedness
        @test nipar ≈ niperp × incident
        
        # Test normalization
        @test norm(incident) ≈ 1.0
        @test norm(niperp) ≈ 1.0
        @test norm(nipar) ≈ 1.0
    end
    
    @testset "construct_transmission_coordinate_system" begin
        view = normalize(SA[1.0, 1.0, 2.0])
        normal = SA[0.0, 0.0, 1.0]
        
        ntperp, ntpar = CoRaLS.construct_transmission_coordinate_system(view, normal)
        
        # Test orthogonality and right-handedness
        @test abs(ntperp ⋅ ntpar) < 1e-10
        @test abs(ntperp ⋅ view) < 1e-10
        @test abs(ntpar ⋅ view) < 1e-10
        @test ntpar ≈ ntperp × view
        
        # Test normalization
        @test norm(ntperp) ≈ 1.0
        @test norm(ntpar) ≈ 1.0
    end
    
    @testset "calculate_transmitted_polarization" begin
        # Test that Fresnel coefficients are correctly applied
        pol = SA[1.0, 0.0, 0.0]
        niperp = SA[1.0, 0.0, 0.0]
        nipar = SA[0.0, 1.0, 0.0]
        ntperp = SA[1.0, 0.0, 0.0]
        ntpar = SA[0.0, 1.0, 0.0]
        tperp = 0.8
        tpar = 0.6
        
        poltr = CoRaLS.calculate_transmitted_polarization(pol, niperp, nipar, ntperp, ntpar, tperp, tpar)
        
        # Should be tperp * ntperp since pol is aligned with niperp
        @test poltr ≈ SA[0.8, 0.0, 0.0]
        
        # Test with mixed polarization
        pol_mixed = normalize(SA[1.0, 1.0, 0.0])
        poltr_mixed = CoRaLS.calculate_transmitted_polarization(pol_mixed, niperp, nipar, ntperp, ntpar, tperp, tpar)
        expected = tperp * (pol_mixed ⋅ niperp) * ntperp + tpar * (pol_mixed ⋅ nipar) * ntpar
        @test poltr_mixed ≈ expected
    end
    
    @testset "calculate_polarization_angle" begin
        ntperp = SA[1.0, 0.0, 0.0]
        ntpar = SA[0.0, 1.0, 0.0]
        
        # Test cardinal angles
        @test CoRaLS.calculate_polarization_angle(ntperp, ntpar, ntperp) ≈ 0.0 atol=1e-10
        @test CoRaLS.calculate_polarization_angle(ntpar, ntpar, ntperp) ≈ π/2 atol=1e-10
        @test CoRaLS.calculate_polarization_angle(-ntperp, ntpar, ntperp) ≈ π atol=1e-10
        
        # Test 45° angle
        poltr_45 = normalize(ntperp + ntpar)
        @test CoRaLS.calculate_polarization_angle(poltr_45, ntpar, ntperp) ≈ π/4 atol=1e-10
    end
    
    @testset "calculate_observation_angles" begin
        origin = SA[0.0km, 0.0km, 1737.4km]
        axis = SA[0.0, 1.0, 0.0]
        normal = origin / norm(origin)
        view = SA[0.0, 0.0, 1.0]
        antenna = SA[0.0km, 0.0km, 1757.4km]
        
        zenith, θ, ϕ, el = CoRaLS.calculate_observation_angles(origin, axis, normal, view, antenna)
        
        # Test valid ranges
        @test 0.0 <= zenith <= π
        @test 0.0 <= θ <= π
        @test -π <= ϕ <= 2π
        @test -π/2 <= el <= π/2
        
        # Test specific geometry: origin at north pole
        origin_pole = SA[0.0km, 0.0km, 1737.4km]
        zenith_pole, θ_pole, _, _ = CoRaLS.calculate_observation_angles(
            origin_pole, -SA[0.0, 0.0, 1.0], SA[0.0, 0.0, 1.0], view, antenna)
        @test θ_pole ≈ 0.0 atol=1e-6  # At north pole
    end
    
    @testset "compute_direct integration" begin
        # Test that compute_direct returns Direct struct for valid geometry
        Ecr = 1e18 * 1.602e-10 * EeV  # 1 EeV
        origin = SA[0.0km, 0.0km, 1737.4km]
        axis = normalize(SA[0.0, 1.0, -1.0])
        antenna = SA[0.0km, 0.0km, 1757.4km]
        Xmax = SA[0.0km, 0.0km, 1737.395km]
        
        result = CoRaLS.compute_direct(CoRaLS.ScalarGeometry(), Ecr, origin, axis, antenna, Xmax)
        
        # Test returns Direct struct (not TIR or other failure)
        @test result isa CoRaLS.Direct
        @test result.Ecr ≈ Ecr
        @test result.depth > 0.0m
        @test result.Drego > 0.0m
        @test result.Dvacuum > 0.0m
        @test result.triggered == false  # Default value
        
        # Test TIR case (extreme angle)
        antenna_tir = SA[0.0km, 2000.0km, 1737.4km]  # Far to the side
        result_tir = CoRaLS.compute_direct(CoRaLS.ScalarGeometry(), Ecr, origin, axis, antenna_tir, Xmax)
        @test result_tir == CoRaLS.TIR
    end
    
    @testset "compute_reflected integration" begin
        # Test that compute_reflected returns Reflected struct for valid geometry
        Ecr = 1e18 * 1.602e-10 * EeV
        origin = SA[0.0km, 0.0km, 1737.4km]
        axis = normalize(SA[0.0, 1.0, -1.0])
        antenna = SA[0.0km, 0.0km, 1757.4km]
        Xmax = SA[0.0km, 0.0km, 1737.397km]  # 3m deep
        ice_depth = 10.0m
        
        result = CoRaLS.compute_reflected(CoRaLS.ScalarGeometry(), Ecr, origin, axis, antenna, Xmax, ice_depth)
        
        # Test returns Reflected struct
        @test result isa CoRaLS.Reflected
        @test result.Ecr ≈ Ecr
        @test result.depth > 0.0m
        @test result.depth < ice_depth  # Xmax before ice
        @test result.Drego > 0.0m
        @test result.θ_ice >= 0.0  # Ice angle calculated
        
        # Test subsurface reflection fields exist
        @test result.rpar != 0.0
        @test result.rperp != 0.0
    end
    
end