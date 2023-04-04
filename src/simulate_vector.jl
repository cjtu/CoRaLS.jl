"""
    Compute the 'direct' RF solution using the vector geometry.
"""
function compute_direct(::VectorGeometry,
                        Ecr, origin, axis, antenna, Xmax;
                        indexmodel=SurfaceDeepIndex(),
                        densitymodel=StrangwayDensity(),
                        fieldmodel=ARW(),
                        divergencemodel=MixedFieldDivergence(),
                        slopemodel=NoSlope(),
                        roughnessmodel=NoRoughness(),
                        kwargs...)

     @assert (slopemodel isa NoSlope) "Vector geometry does not currently support a slope model."

    # we use a minimization routine to find the launch vector such that
    # after we refract at the surface, the RF emission hits the spacecraft.

    # our initial guess is the straight vector from Xmax to the antenna
    initial = (antenna - Xmax) / norm(antenna - Xmax)

    # use a Newton-Rhapson minimizer to find the correct launch angle
    final_scale = optimize(scale -> propagate_and_refract(Xmax, initial, antenna, scale, indexmodel)[end],
                           [0.05], LBFGS(),
                           Optim.Options(x_tol=deg2rad(1e-3))).minimizer

    # calculate the emission vector and the surface
    emit, surface, θ_i, θ_r, residuals = propagate_and_refract(Xmax, initial, antenna, final_scale, indexmodel)

    # the distance in the regolith is between Xmax and the surface
    Drego = norm(Xmax - surface)

    # and the total distance in vacuum is from surface to the antenna
    Dvacuum = norm(surface - antenna)

    # the off-axis angle is between axis and emit - both are normalized
    ψ = acos( emit ⋅ axis )

    # calculate the depth of Xmax
    depth = norm(origin) - norm(Xmax)

    # calculate the electric field at the surface
    ν, E = regolith_field(fieldmodel, Ecr, ψ, Drego, Dvacuum;
                          n=regolith_index(indexmodel, depth),
                          density=regolith_density(densitymodel, depth), kwargs...)

    # we need the vector from the surface to the payload
    view = (antenna - surface) / norm(antenna - surface)

    # calculate the incident polarization vector for the direct emission
    pol = (axis × emit) × emit
    pol /= norm(pol) # fixes some tiny rounding errors

    # this vector defines "perpendicular" to the plane of incidence
    # this is "slappy" or H-pol
    nperp = view × normal
    nperp /= norm(nperp)

    # this vector defines "parallel" to the plane of incident i.e. "pokey"
    # again assuming that "origin" === "surface"
    npar = normal

    # calculate the modified Fresnel coefficients for transmission
    # use the refractive index at the surface
    tpar = divergence_tpar(divergencemodel, θ_i, θ_r, regolith_index(indexmodel, 0.0m), Drego, Dvacuum)
    tperp = divergence_tperp(divergencemodel, θ_i, θ_r, regolith_index(indexmodel, 0.0m), Drego, Dvacuum)

    # project `pol` onto `npar` and `nperp`, apply Fresnel,
    # and then recombine into the transmitted polarization vector
    poltr = tperp*(pol ⋅ nperp)*nperp + tpar*(pol ⋅ npar)*npar

    # and construct the integrated electric field vector
    Ef = integrate(ν, E) * poltr

    # calculate the zenith angle of the cosmic ray
    # flip the axis so we get the complement of the angle
    zenith = acos( (origin / norm(origin)) ⋅ (-axis))

    # get the Lunar-centric angle of the cosmic ray impact point
    θ, ϕ, _ = cartesian_to_spherical(origin...)

    # calculate the below horizon angle at the payload
    el = -(acos( -view ⋅ (antenna / norm(antenna)) ) - pi/2.)

    # and construct and return the signal
    return Direct(Ecr, rad2deg(θ), rad2deg(ϕ), rad2deg(zenith),
                  Ef .|> (μV/m), rad2deg(el), rad2deg(ψ),
                  depth, Drego, Dvacuum, rad2deg(θ_i), tpar, tperp)

end

"""
Compute the 'reflected' RF solution using the vector solution.
"""
function compute_reflected(::VectorGeometry,
                           Ecr, origin, axis, antenna, Xmax, ice_depth;
                           indexmodel=SurfaceDeepIndex(),
                           densitymodel=StrangwayDensity(),
                           fieldmodel=ARW(),
                           divergencemodel=MixedFieldDivergence(), Nice=1.375,
                           slopemodel=NoSlope(),
                           roughnessmodel=NoRoughness(),
                           iceroughness=NoIceRoughness(),
                           kwargs...)

    @assert (slopemodel isa NoSlope) "Vector geometry does not currently support a slope model."

    # we use a minimization routine to find the launch vector such that
    # after we refract at the surface, the RF emission hits the spacecraft.

    # we need the normalized vector pointing towards Xmax
    normal = Xmax / norm(Xmax)

    # we reflected Xmax below the ice-layer - this is the apparent starting point
    Xmax_refl = 2.0*(Rmoon - ice_depth)*normal - Xmax

    # check that the reflection condition is satisfied
    # @toggled_assert (Rmoon - ice_depth) - norm(Xmax_refl) == norm(Xmax) - ice_depth

    # our initial guess is the straight vector from Xmax to the antenna
    initial = (antenna - Xmax_refl) / norm(antenna - Xmax_refl)

    # use a Newton-Rhapson minimizer to find the correct launch angle
    final_scale = optimize(scale -> propagate_and_refract(Xmax_refl, initial, antenna, scale, indexmodel)[end],
                           [0.05], LBFGS(),
                           Optim.Options(x_tol=deg2rad(1e-3))).minimizer

    # calculate the emission vector and the surface
    emit, surface, θ_i, θ_r, residuals = propagate_and_refract(Xmax_refl, initial, antenna, final_scale, indexmodel)

    # a quick check that Snell's law is satisfied
    # @toggled_assert regolith_index(indexmodel, 0.0m)*sin(θ_i) ≈ sin(θ_r)

    # the distance in the regolith is between Xmax and the surface
    Drego = norm(Xmax_refl - surface)

    # and the total distance in vacuum is from surface to the antenna
    Dvacuum = norm(surface - antenna)

    # calculate the reflected cosmic ray axis
    axis_refl = axis - 2.0(axis ⋅ normal)*normal

    # the off-axis angle is between axis and emit - both are normalized
    ψ = acos( emit ⋅ axis_refl )

    # calculate the depth of Xmax
    depth = norm(origin) - norm(Xmax)

    # calculate the electric field at the surface
    ν, E = regolith_field(fieldmodel, Ecr, ψ, Drego, Dvacuum;
                          n=regolith_index(indexmodel, depth),
                          density=regolith_density(densitymodel, depth), kwargs...)

    # apply roughness to the simulated electric field
    ν, E = ice_roughness(iceroughness, ν, E, θ_i)

    # we need the vector from the surface to the payload
    view = (antenna - surface) / norm(antenna - surface)

    # calculate the incident polarization vector for the direct emission
    pol = (axis_refl × emit) × emit
    pol /= norm(pol) # fixes some tiny rounding errors

    # this vector defines "perpendicular" to the plane of incidence
    # this is "slappy" or H-pol
    nperp = view × normal
    nperp /= norm(nperp)

    # this vector defines "parallel" to the plane of incident i.e. "pokey"
    # again assuming that "origin" === "surface"
    npar = normal

    # calculate the modified Fresnel coefficients for transmission
    # use the refractive index at the surface
    tpar = divergence_tpar(divergencemodel, θ_i, θ_r, regolith_index(indexmodel, 0.0m), Drego, Dvacuum)
    tperp = divergence_tperp(divergencemodel, θ_i, θ_r, regolith_index(indexmodel, 0.0m), Drego, Dvacuum)

    # to be 100% correct, propagate along `emit` to find the incident angle
    # of the RF at the ice surface
    # this correction is *tiny* compared to just using θ_i - on the scale
    # of a thousand-th of a degree. We can probably drop this to save some CPU time.
    ice_hit = intersect_with_sphere(Xmax_refl, emit, Rmoon - ice_depth)

    # calculate the incident angle at the ice-surface
    θ_ice = acos( (ice_hit / norm(ice_hit)) ⋅ emit )

    # calculate the Fresnel reflection coefficients at the ice
    rpar = fresnel_rpar(θ_ice, regolith_index(indexmodel, ice_depth), Nice)
    rperp = fresnel_rperp(θ_ice, regolith_index(indexmodel, ice_depth), Nice)

    # project `pol` onto `npar` and `nperp`, apply Fresnel,
    # and then recombine into the transmitted polarization vector
    poltr = rperp*tperp*(pol ⋅ nperp)*nperp + rpar*tpar*(pol ⋅ npar)*npar

    # and construct the integrated electric field vector
    Ef = integrate(ν, E) * poltr

    # calculate the zenith angle of the cosmic ray
    # flip the axis so we get the complement of the angle
    zenith = acos( (origin / norm(origin)) ⋅ (-axis))

    # get the Lunar-centric angle of the cosmic ray impact point
    θ, ϕ, _ = cartesian_to_spherical(origin...)

    # calculate the below horizon angle at the payload
    el = -(acos( -view ⋅ (antenna / norm(antenna)) ) - pi/2.)

    # and construct and return the signal
    return Reflected(Ecr, rad2deg(θ), rad2deg(ϕ), rad2deg(zenith),
                     Ef .|> (μV/m), rad2deg(el), rad2deg(ψ),
                     depth, Drego, Dvacuum, rad2deg(θ_i), tpar, tperp, rpar, rperp)

end
