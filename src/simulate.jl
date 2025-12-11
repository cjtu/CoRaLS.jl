using Optim
using Distributions
using StructArrays
using LinearAlgebra
using NumericalIntegration
using Unitful
using Unitful: m, km, μV, sr
"""
This script implements various functions and types for simulating cosmic ray events and their corresponding radio frequency (RF) signals. Below is a summary of the main components of the script:

1. Enums:
   - `TrialFailed`: Enumerates possible reasons why a cosmic ray trial failed.
   - `EventGeometry`: Enumerates possible event geometries (direct or reflected).

2. Abstract Types:
   - `GeometryImplementation`: Abstract type for event simulation implementation.
   - `AbstractSignal`: Abstract type for simulated signals, both direct and reflected.

3. Concrete Types:
   - `VectorGeometry`: 3D vector geometry implementation without approximations.
   - `ScalarGeometry`: Approximate scalar geometry based on Andres' MC.

4. Signal Structs:
   - `Direct`: Information stored for a direct RF detection.
   - `Reflected`: Information stored for a reflected RF detection.

5. Functions:
   - `throw_cosmicray`: Simulates a single cosmic ray trial with a given energy, calculating the field at the payload for direct and reflected emissions.
   - `compute_direct`: Computes the 'direct' RF solution using the scalar geometry.
   - `compute_reflected`: Computes the 'reflected' RF solution using the scalar geometry.

6. Miscellaneous:
   - Overloaded `show` function for custom printing of `Direct` and `Reflected` structs.
"""

"""
An enum for possible reasons why the a cosmic ray trial failed.
"""
@enum TrialFailed TIR = 1 XmaxAfterIce = 2 NoXmax = 3 Upgoing = 4 NotVisible = 5 NotInRegion = 6

"""
An enum for possible event geometries (i.e. direct or reflected)
"""
@enum EventGeometry DirectEvent = 0 ReflectedEvent = 1

"""
Abstract type for event simulation implementation.
"""
abstract type GeometryImplementation end

"""
3D vector geometry implementation.

This is the original method of the UH MC - there are no
approximations used in this solution.
"""
struct VectorGeometry <: GeometryImplementation end

"""
*Approximate* scalar geometry based on Andres' MC.

This makes the approximation that the CR impact point is
the same as the refraction point - this is an error of
order <10m out of the 30km distance to the payload.
"""
struct ScalarGeometry <: GeometryImplementation end

"""
Abstract type for simulated signals, both direct and reflected.
"""
abstract type AbstractSignal end

"""
The information stored for a direct RF detection.
"""
mutable struct Direct <: AbstractSignal
    Ecr::typeof(0.0EeV) # the energy of the cosmic ray
    θ::Float64 # the Lunar-centric angle of the cosmic ray (deg)
    ϕ::Float64 # the Lunar-centric azimuthal angle (deg)
    θ_z::Float64 # the CR zenith angle (deg) at θ
    pol::typeof(SA[0.0, 0.0, 0.0]) # the magnitude & polarization vector at the payload
    ν_min::typeof(1.0MHz) # the minimum frequency of this e-field
    dν::typeof(1.0MHz) # the frequency bin spacing of this e-field
    ν_max::typeof(1.0MHz) # the maximum frequency of this e-field
    Ef # the total electric field spectral density at the payload
    θpol::Float64 # the polarization angle (deg) at the payload, 0 for H-Pol, 90 for V-Pol
    θ_el::Float64 # the observed elevation angle below the payload (deg)
    ψ::Float64 # the off-axis angle (deg)
    depth::typeof(0.0km) # the depth below the surface of Xmax
    Drego::typeof(0.0km) # the distance in the regolith
    Dvacuum::typeof(0.0km) # the distance in vacum
    θ_i::Float64 # the incident angle at the surface
    θ_emit::Float64 # the zenith angle of the emitted ray (deg)
    tpar::Float64 # the parallel Fresnel transmission coefficient
    tperp::Float64 # the perpendicular Fresnel transmission coefficient
    triggered::Bool # did this event trigger
end

"""
The information stored for a reflected RF detection.
"""
mutable struct Reflected <: AbstractSignal
    Ecr::typeof(0.0EeV) # the energy of the cosmic ray
    θ::Float64 # the Lunar-centric angle of the cosmic ray (deg)
    ϕ::Float64 # the Lunar-centric azimuthal angle (deg)
    θ_z::Float64 # the zenith angle (deg) at θ
    pol::typeof(SA[0.0, 0.0, 0.0]) # the magnitude & polarization vector
    ν_min::typeof(1.0MHz) # the minimum frequency of this e-field
    dν::typeof(1.0MHz) # the frequency bin spacing of this e-field
    ν_max::typeof(1.0MHz) # the maximum frequency of this e-field
    Ef # the total electric field spectral density at the payload
    polsub::typeof(SA[0.0, 0.0, 0.0]) # the polarization vector from the subsurface
    θpol::Float64 # the polarization angle (deg) at the payload, 0 for H-Pol, 90 for V-Pol
    θpolsub::Float64 # the polarization angle (deg) at the payload, 0 for H-Pol, 90 for V-Pol
    θ_el::Float64 # the observed elevation angle below the payload (deg)
    ψ::Float64 # the off-axis angle (deg)
    depth::typeof(0.0km) # the depth below the surface of Xmax
    Drego::typeof(0.0km) # the distance in the regolith
    Dvacuum::typeof(0.0km) # the distance in vacum
    θ_i::Float64 # the incident angle at the surface (deg)
    θ_emit::Float64 # the zenith angle of the emitted ray (deg)
    θ_ice::Float64 # the incident angle at the ice (deg)
    tpar::Float64 # the parallel Fresnel transmission coefficient
    tperp::Float64 # the perpendicular Fresnel transmission coefficient
    rpar::Float64 # the parallel Fresnel reflection coefficient
    rperp::Float64 # the perpendicular Fresnel reflection coefficient
    subpar::Float64 # the parallel Fresnel coefficient from the subsurface refl.
    subperp::Float64 # the perpendicular Fresnel coefficient from the subsurface refl.
    triggered::Bool # did this event trigger
end

"""
Simplifiy the printing of Signal structs.
"""
Base.show(io::IO, event::Direct) = print(io, "DirectEvent($(event.Ecr), $(event.θ_el))")
Base.show(io::IO, event::Reflected) = print(io, "ReflectedEvent($(event.Ecr), $(event.θ_el))")


"""
Simulate a single cosmic ray trial with a given energy.

This calculates the field at the payload for the direct and reflected emission,
(if they exist) and returns whether the event triggered.
"""
function throw_cosmicray(Ecr, trigger, region, spacecraft; simple_area=false, kwargs...)
    SC = get_position(spacecraft)
    if simple_area && region isa AOIRegion
        surface = random_point_in_aoi(region.aoi)
        if !is_visible(surface, SC)
            return (NotVisible, NotVisible)
        end
    else
        surface = random_point_on_sphere(Rmoon)
        if !is_visible(surface, SC)
            return (NotVisible, NotVisible)
        elseif !is_in_region(surface, region)
            return (NotInRegion, NotInRegion)
        end
    end
    return propagate_cosmicray(Ecr, surface, SC, trigger; kwargs...)
end


"""
Simulate a single cosmic ray trial with a given energy.

This calculates the field at the payload for the direct and reflected emission,
(if they exist) and returns whether the event triggered.
"""
function throw_cosmicray(Ecr, trigger; altitude=20.0km, kwargs...)

    # get the maximum angle that we sample CR impact points from
    θmax = -horizon_angle(altitude)

    # and assume that the polar regions extend out to 10 degrees from the poles
    θpole = deg2rad(10.0)

    # sample a point on the surface of North or South Poles
    # surface = rand() > 0.5 ? random_north_pole_point() : random_south_pole_point()

    # check if this point lies within a PSR
    # within_psr = point_impacts_psr(surface)

    # if we didn't hit a PSR, then we throw away this trial
    # within_psr == false && return NotInRegion, NotInRegion

    # throw for a random origin point on a single hemisphere
    # surface = random_point_on_cap(θmax; r=Rmoon)
    surface = random_point_on_cap(θpole; r=Rmoon)
    # surface = random_point_on_cap(π; r=Rmoon)

    # calculate the sampled value of θ for the surface point
    θ, _, _ = cartesian_to_spherical(surface)

    # check if the point lies outside the polar cap
    # ((θ > θpole) && θ < (π - θpole)) && return (NotInRegion, NotInRegion)
    θ > θpole && return (NotInRegion, NotInRegion)

    # draw a random polar angle for the spacecraft location
    # the edges of the polar cap are visible from twice its
    # opening angle away from the pole.
    # This is drawing an angle from a *circular orbit* so its
    # uniform in θ *not* uniform in cos(θ) like an area would be.
    # we assume this is symmetric around the poles so we still
    # only simulate one of the poles.
    θsc = rand(Uniform(0, θpole + θmax))

    # if the SC is not in view of the lunar pols
    # ((θsc > 2.0*θpole) && θsc < (π - 2.0*θpole)) && (NotVisible, NotVisible)
    # θsc > 2.0*θpole && (NotVisible, NotVisible)

    # and finally calculate the 3D position of the spacecraft
    # SC = spherical_to_cartesian(θsc, 0., Rmoon+altitude)
    SC = spherical_to_cartesian(θsc, rand(Uniform(0, 2π)), Rmoon + altitude)

    # calculate the total central angle between the SC and the event
    Δσ = atan(norm(SC × surface), (SC ⋅ surface))

    # if this point is not within the horizon of the SC, then we can't see it.
    # abs(θ - θsc) > θmax && return (NotVisible, NotVisible)
    abs(Δσ) > θmax && return (NotVisible, NotVisible)
    return propagate_cosmicray(Ecr, surface, SC, trigger; kwargs...)
end

"""
Propagate cosmic ray into surface and compute direct and reflected RF.
"""
function propagate_cosmicray(Ecr, surface, SC, trigger;
    ice_depth=5m,
    geometrymodel=ScalarGeometry(),
    indexmodel=StrangwayIndex(), fieldmodel=ARW(),
    divergencemodel=MixedFieldDivergence(),
    densitymodel=StrangwayDensity(),
    slopemodel=GaussianSlope(7.6),
    roughnessmodel=GaussianRoughness(2.0),
    iceroughness=GaussianIceRoughness(2.0cm),
    kwargs...)
    # throw for a random incident cosmic ray direction in cos^2
    direction = random_direction(surface / norm(surface))

    # if it's an upgoing cosmic ray, reject the trial
    if surface ⋅ direction > 0km
        @debug "Got an upgoing cosmic ray trial! This should never happen."
        return Upgoing, Upgoing
    end

    # we now step the cosmic ray into the regolith until Xmax
    Xmax = propagate_to_Xmax(surface, direction, Ecr, densitymodel)

    # check that we didn't "leave" the regolith due to a highly
    # inclined shower being "too" long.
    if norm(Xmax) > Rmoon
        return NoXmax, NoXmax
    end

    # find the depth of the cosmic ray
    depth = norm(surface) - norm(Xmax)

    # check that Xmax occured before the ice layer
    if depth > ice_depth
        return XmaxAfterIce, XmaxAfterIce
    end

    # compute the solution for the "direct" RF
    direct = compute_direct(geometrymodel, Ecr, surface, direction, SC, Xmax;
        indexmodel=indexmodel,
        fieldmodel=fieldmodel,
        slopemodel=slopemodel,
        roughnessmodel=roughnessmodel,
        densitymodel=densitymodel,
        divergencemodel=divergencemodel,
        kwargs...)

    # compute the solution for the "reflected" RF
    reflected = compute_reflected(geometrymodel, Ecr, surface, direction, SC, Xmax, ice_depth;
        indexmodel=indexmodel,
        fieldmodel=fieldmodel,
        slopemodel=slopemodel,
        roughnessmodel=roughnessmodel,
        iceroughness=iceroughness,
        densitymodel=densitymodel,
        divergencemodel=divergencemodel,
        kwargs...)

    return direct, reflected
end


"""
Calculate distances in regolith and vacuum for direct ray path.

Returns: (Drego, Dvacuum, x) where x is horizontal distance to Xmax.
"""
function calculate_direct_distances(depth, θ_reg, obs)
    x = depth * tan(θ_reg)
    Drego = sqrt(x * x + depth * depth)
    Dvacuum = norm(obs)
    return Drego, Dvacuum, x
end

"""
Calculate emission vector and angle at Xmax.

Returns: (emit, θ_emit) where emit is normalized emission direction.
"""
function calculate_emission_vector(normal, proj, depth, θ_i, Nsurf, NXmax, Rmoon, is_reflected=false)
    θ_emit = snell_θt(θ_i, Nsurf, NXmax)
    
    # For reflected, emission goes opposite direction (down then up) so vert component is negative
    down = is_reflected ? -1 : 1
    emit = down * cos(θ_emit) * normal + sin(θ_emit) * proj
    emit /= norm(emit)
    return emit, θ_emit
end

"""
Calculate incident polarization vector for direct or reflected emission.
"""
function calculate_incident_polarization(emit, axis)
    pol = -emit × (emit × axis)
    pol /= norm(pol)
    return pol
end

"""
Construct coordinate system for incident ray at surface.

Returns: (incident, niperp, nipar) where incident is the ray direction
and niperp/nipar define the perpendicular and parallel basis vectors.
"""
function construct_incident_coordinate_system(normal, proj, θ_i)
    incident = cos(θ_i) * normal + sin(θ_i) * proj
    incident /= norm(incident)
    
    niperp = incident × normal
    niperp /= norm(niperp)
    
    nipar = niperp × incident
    nipar /= norm(nipar)
    
    return incident, niperp, nipar
end

"""
Construct coordinate system for transmitted ray (at spacecraft).

Returns: (ntperp, ntpar) basis vectors for transmitted wave.
"""
function construct_transmission_coordinate_system(view, normal)
    ntperp = view × normal
    ntperp /= norm(ntperp)
    
    ntpar = ntperp × view
    ntpar /= norm(ntpar)
    
    return ntperp, ntpar
end

"""
Calculate transmitted polarization vector after Fresnel transmission.

Projects incident polarization onto incident basis, applies Fresnel coefficients,
then projects onto transmission basis.
"""
function calculate_transmitted_polarization(pol, niperp, nipar, ntperp, ntpar, tperp, tpar)
    poltr = tperp * (pol ⋅ niperp) * ntperp + tpar * (pol ⋅ nipar) * ntpar
    return poltr
end

"""
Calculate polarization angle from transmitted polarization vector.

Returns angle in radians where 0 is perpendicular (H-Pol), π/2 is parallel (V-Pol).
"""
function calculate_polarization_angle(poltr, ntpar, ntperp)
    return atan(poltr ⋅ ntpar, poltr ⋅ ntperp)
end

"""
Calculate distances in regolith and vacuum for reflected ray path.

Returns: (Drego, Dvacuum, x1, x2, source) where:
- Drego: total distance in regolith (pre + post reflection)
- Dvacuum: distance in vacuum to spacecraft
- x1: horizontal distance pre-reflection
- x2: horizontal distance post-reflection
- source: approximate location of Xmax in coordinate system
"""
function calculate_reflected_distances(depth, ice_depth, θ_reg, obs, origin, normal, proj)
    x1 = (ice_depth - depth) * tan(θ_reg)  # pre-reflected
    x2 = ice_depth * tan(θ_reg)  # post-reflection
    
    Drego = sqrt(x1 * x1 + (ice_depth - depth) * (ice_depth - depth)) +
            sqrt(x2 * x2 + ice_depth * ice_depth)
    
    Dvacuum = norm(obs)
    
    # Calculate approximate location of Xmax
    source = origin - (2.0 * ice_depth - depth) * normal - (x1 + x2) * proj
    
    return Drego, Dvacuum, x1, x2, source
end

"""
Calculate observation angles (cosmic ray zenith angle, lunar coordinates, and elevation angle).

Returns: (zenith, θ, ϕ, el) all in radians except θ, ϕ from spherical conversion.
"""
function calculate_observation_angles(origin, axis, normal, view, antenna)
    zenith = acos(normal ⋅ (-axis))
    θ, ϕ, _ = cartesian_to_spherical(origin)
    el = -(acos(-view ⋅ (antenna / norm(antenna))) - π / 2.0)
    return zenith, θ, ϕ, el
end

"""
Initialize common event geometry shared by direct and reflected computations.

Extracts surface geometry, depth, refractive indices, and projection calculations
that are identical in both direct and reflected signal paths.

Returns: (normal, surface_normal, obs, view, depth, Nsurf, NXmax, proj)
"""
function init_event_geometry(origin, antenna, Xmax, indexmodel, slopemodel)
    # Calculate basic surface geometry
    normal = origin / norm(origin)

    # Get new normal relative to random sloped surface (due to topography)
    normal_sloped = random_surface_normal(slopemodel, normal)
    obs = antenna - origin
    view = obs / norm(obs)
    
    # Calculate depth and refractive indices
    depth = norm(origin) - norm(Xmax)
    Nsurf = regolith_index(indexmodel, 0.0m)
    NXmax = regolith_index(indexmodel, depth)
    
    # Project view onto local horizontal for Xmax shift calculations
    proj = view - (view ⋅ normal) * normal
    proj /= norm(proj)
    
    return normal_sloped, obs, view, depth, Nsurf, NXmax, proj
end

"""
Compute the 'direct' RF solution using the scalar geometry.
"""
function compute_direct(::ScalarGeometry,
    Ecr, origin, axis, antenna, Xmax;
    indexmodel=SurfaceDeepIndex(),
    densitymodel=StrangwayDensity(),
    fieldmodel=ARW(),
    divergencemodel=MixedFieldDivergence(),
    slopemodel=NoSlope(),
    roughnessmodel=NoRoughness(),
    ν_min=150MHz, ν_max=600MHz,
    kwargs...)

    # Initialize common event geometry
    normal, obs, view, depth, Nsurf, NXmax, proj = 
        init_event_geometry(origin, antenna, Xmax, indexmodel, slopemodel)
    
    # Calculate exit angle and check for TIR
    θ_r = acos(view ⋅ normal)
    θ_r > π / 2.0 && return TIR

    # Calculate incident and regolith angles using Snell's law
    θ_i = snell_θt(θ_r, 1.0, Nsurf)
    θ_reg = snell_θt(θ_i, Nsurf, NXmax)

    # Calculate distances in regolith and vacuum
    Drego, Dvacuum, x = calculate_direct_distances(depth, θ_reg, obs)

    # Calculate emission vector and angles
    emit, θ_emit = calculate_emission_vector(normal, proj, depth, θ_i, Nsurf, NXmax, Rmoon, false)
    ψ = acos(emit ⋅ axis)

    # Calculate the electric field at the surface
    ν, E = regolith_field(fieldmodel, Ecr, ψ, Drego, Dvacuum;
        n=NXmax,
        density=regolith_density(densitymodel, depth),
        ν_min=ν_min, ν_max=ν_max, kwargs...)

    # Calculate incident polarization vector
    pol = calculate_incident_polarization(emit, axis)

    # Construct coordinate systems for incident and transmitted rays
    incident, niperp, nipar = construct_incident_coordinate_system(normal, proj, θ_i)
    ntperp, ntpar = construct_transmission_coordinate_system(view, normal)

    # Calculate Fresnel transmission coefficients
    tpar, tperp = surface_transmission(roughnessmodel, divergencemodel,
        θ_i, Nsurf, Drego, Dvacuum)

    # Calculate transmitted polarization vector and angle
    poltr = calculate_transmitted_polarization(pol, niperp, nipar, ntperp, ntpar, tperp, tpar)
    θpol = calculate_polarization_angle(poltr, ntpar, ntperp)

    # Calculate observation angles
    zenith, θ, ϕ, el = calculate_observation_angles(origin, axis, normal, view, antenna)

    # and construct and return the signal
    return Direct(Ecr, rad2deg(θ), rad2deg(ϕ), rad2deg(zenith),
        poltr, ν_min, 10MHz, ν_max, E .|> (μV / m / MHz), rad2deg(θpol), rad2deg(el), rad2deg(ψ),
        depth, Drego, Dvacuum, rad2deg(θ_i), rad2deg(θ_emit), tpar, tperp, false)

end


"""
Compute the 'reflected' RF solution using the scalar solution.
"""
function compute_reflected(::ScalarGeometry,
    Ecr, origin, axis, antenna, Xmax, ice_depth;
    Nice=1.305,
    Nbed=1.6, # 1.6 for regolith, 2.8 for bedrock
    ice_thickness=1.0m,
    indexmodel=SurfaceDeepIndex(),
    densitymodel=StrangwayDensity(),
    fieldmodel=ARW(),
    divergencemodel=MixedFieldDivergence(),
    slopemodel=NoSlope(),
    roughnessmodel=NoRoughness(),
    iceroughness=NoIceRoughness(),
    ν_min=150MHz, ν_max=600MHz,
    kwargs...)

    # Initialize common event geometry
    normal, obs, view, depth, Nsurf, NXmax, proj = 
        init_event_geometry(origin, antenna, Xmax, indexmodel, slopemodel)
    
    # Calculate exit angle (for geometry) and actual surface angle (for TIR check)
    θ_r = acos(view ⋅ normal)  # Geometric angle
    θ_r_true = acos(view ⋅ normal)  # Actual surface angle
    θ_r_true > π / 2.0 && return TIR

    # Calculate additional refractive index at ice depth
    Nrego_at_ice = regolith_index(indexmodel, ice_depth)

    # Calculate incident and regolith angles using Snell's law
    θ_i = snell_θt(θ_r, 1.0, Nsurf)
    θ_i_true = snell_θt(θ_r_true, 1.0, Nsurf)
    θ_reg = snell_θt(θ_i_true, Nsurf, NXmax)

    # Calculate distances for reflected path
    Drego, Dvacuum, x1, x2, source = calculate_reflected_distances(
        depth, ice_depth, θ_reg, obs, origin, normal, proj)

    # Check for TIR in the reflected geometry
    (x1 + x2) > (2.0 * ice_depth - depth) / sqrt(regolith_index(indexmodel, 0.0m)^2.0 - 1) && return TIR

    # Calculate emission vector and angles (reflected case)
    emit, θ_emit = calculate_emission_vector(normal, proj, depth, θ_i_true, Nsurf, NXmax, Rmoon, true)
    ψ = acos(emit ⋅ axis)

    # Calculate the electric field at the spacecraft
    ν, E = regolith_field(fieldmodel, Ecr, ψ, Drego, Dvacuum;
        n=NXmax,
        density=regolith_density(densitymodel, depth),
        ν_min=ν_min, ν_max=ν_max, kwargs...)

    # Calculate incident polarization vector
    pol = calculate_incident_polarization(emit, axis)

    # Calculate angle at ice layer using Snell's law
    θ_ice = snell_θt(θ_i_true, Nsurf, Nrego_at_ice)

    # Construct coordinate systems for incident and transmitted rays
    incident, niperp, nipar = construct_incident_coordinate_system(normal, proj, θ_i)
    ntperp, ntpar = construct_transmission_coordinate_system(view, normal)

    # Apply ice roughness to the electric field
    ν, E = ice_roughness(iceroughness, ν, E, θ_ice, NXmax)

    # Calculate Fresnel reflection coefficients at the ice
    rpar, rperp = fresnel_coeffs(θ_ice, Nrego_at_ice, Nice)[1:2]

    # Calculate Fresnel transmission coefficients at surface
    tpar, tperp = surface_transmission(roughnessmodel, divergencemodel,
        θ_i_true, Nsurf, Drego, Dvacuum)

    # Calculate transmitted polarization vector and angle
    poltr = rperp * tperp * (pol ⋅ niperp) * ntperp + rpar * tpar * (pol ⋅ nipar) * ntpar
    θpol = calculate_polarization_angle(poltr, ntpar, ntperp)

    # Construct coordinate system and calculate subsurface reflection if needed
    if ice_thickness > 0.0m
        # Get transmission from regolith into ice
        sub_tpar, sub_tperp = fresnel_coeffs(θ_ice, Nrego_at_ice, Nice)[3:4]

        # Get refracted angle in ice layer using Spherical Snell's law
        θ_bed = snell_θt(θ_ice, Nice, Nbed)

        # Get reflection coefficient at ice->regolith or ice->bedrock interface
        bed_rpar, bed_rperp = fresnel_coeffs(θ_bed, Nice, Nbed)[1:2]

        # Get transmission from ice back into regolith
        ice_tpar, ice_tperp = fresnel_coeffs(θ_bed, Nice, Nrego_at_ice)[3:4]

        # Stack all coefficients for subsurface reflection
        subpar = sub_tpar * bed_rpar * ice_tpar
        subperp = sub_tperp * bed_rperp * ice_tperp

        # Construct polarization vector for subsurface reflection
        polsub = subperp * tperp * (pol ⋅ niperp) * ntperp + subpar * tpar * (pol ⋅ nipar) * ntpar
        θpolsub = calculate_polarization_angle(polsub, ntpar, ntperp)
    else
        # No subsurface reflection
        polsub = SA[0.0, 0.0, 0.0]
        θpolsub = 0.0
        subpar = subperp = 0.0
    end

    # Calculate observation angles
    zenith, θ, ϕ, el = calculate_observation_angles(origin, axis, normal, view, antenna)

    # Construct and return the signal
    return Reflected(Ecr, rad2deg(θ), rad2deg(ϕ), rad2deg(zenith),
        poltr,
        ν_min, 10MHz, ν_max,
        E .|> (μV / m / MHz), polsub, rad2deg(θpol), rad2deg(θpolsub), rad2deg(el), rad2deg(ψ),
        depth, Drego, Dvacuum, rad2deg(θ_i), rad2deg(θ_emit), rad2deg(θ_ice),
        tpar, tperp, rpar, rperp, subpar, subperp, false)

end

