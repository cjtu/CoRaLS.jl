
"""
    raytrace(θ_i, Xmax, ice_depth, indexmodel)

    Raytrace a ray from the surface with zenith angle (in the regolith)
    of `θ_i` and calculate the zenith angle at the depth of `Xmax` and
    at `ice_depth` using `indexmodel`.

# Arguments
- `θ_i`: Initial zenith angle of the ray at the lunar surface (in the regolith).
- `Xmax`: Depth at which the cosmic ray reaches its maximum development (Xmax).
- `ice_depth`: Depth of the ice layer beneath the lunar surface.
- `indexmodel`: Model of the regolith index used for ray tracing.

# Returns
- `Drego`: Total path length of the ray within the regolith.
- `θXmax`: Zenith angle of the ray at the Xmax depth.
- `θice`: Zenith angle of the ray at the ice layer depth.

# Description
This function performs ray tracing for a cosmic ray incident on the lunar surface, considering refraction effects as the ray passes through the regolith and approaches the ice layer. The function assumes a perfectly specular reflection at the ice layer, meaning the incident and reflection angles are equal. It computes the path length and zenith angles at crucial depths (Xmax and ice layer) and returns these values along with the total path length within the regolith.
"""
function raytrace(θ_i, Xmax, ice_depth,
    indexmodel::RegolithIndex)

    # make sure Xmax and ice_depth are both measured "positive" below the surface
    @toggled_assert Xmax > 0.0cm
    @toggled_assert ice_depth > 0.0cm && ice_depth > Xmax

    # by default, θXmax and θice are assumed to be θ_i
    # which is the solution without refraction
    θXmax = θ_i
    θice = θ_i

    # this bools make sure that we reach both Xmax and the ice
    hitXmax = false
    hitIce = false

    # we also need to save the total path length up until Xmax
    DXmax = 0.0cm

    # the current depth our ray is at
    depth = 0.0cm

    # the total path length of this ray
    Drego = 0.0cm

    # this is the default radial step that we use for propagation
    dR = 0.1cm

    # the current zenith angle of the ray that is updated in the loop
    θ = θ_i

    # propagate until we reach the depth of the ice
    while depth < ice_depth

        # use different step-sizes based on depth
        # hopefully we never hit the last few branches
        if depth <= 10.0cm
            dR = 0.1cm
        elseif depth <= 20cm
            dR = 0.2cm
        elseif depth <= 50cm
            dR = 0.5cm
        elseif depth <= 100cm
            dR = 1.0cm
        else
            dR = 5.0cm
        end

        # calculate the change in radius corresponding
        # to the current zenith angle of the ray
        dl = dR / cos(θ)

        # get the refractive index at the old and new points
        n_old = regolith_index(indexmodel, depth)
        n_new = regolith_index(indexmodel, depth + dR) # since it is radially symmetric

        # we take advantage of the optical invariant in a
        # radially symmetric refractive index profile (which we assume).
        # n * r * sin(θ) === n' * r' + sin(θ') for all points along the path
        #
        # use the optical invariant to calculate the new zenith angle
        θ = asin((n_old * (Rmoon - depth) * sin(θ)) / ((Rmoon - depth - dR) * n_new))

        # update all the variables
        depth += dR
        Drego += dl

        # check if we are at the location of Xmax
        if abs(depth - Xmax) <= dR
            θXmax = θ
            DXmax = Drego
            hitXmax = true
        end

        # check if we hit the ice layer
        if abs(depth - ice_depth) <= dR
            θice = θ
            hitIce = true
            break
        end

    end

    println("Final Depth: $(depth)")
    println("Final Path: $(Drego)")

    # check that we hit both Xmax and ice
    @assert hitIce && hitXmax

    # println((regolith_index(indexmodel, Xmax) * (Rmoon - Xmax) * sin(θXmax)) |> km)
    # println((regolith_index(indexmodel, 0.0cm) * (Rmoon) * sin(θ_i)) |> km)
    # println((regolith_index(indexmodel, ice_depth) * (Rmoon - ice_depth) * sin(θice)) |> km)

    # and verify that the optical invariant was preserved at Xmax
    # and at the ice layer
    # @assert isapprox(regolith_index(indexmodel, Xmax) * (Rmoon - Xmax) * sin(θXmax),
    #                  regolith_index(indexmodel, 0.0cm) * (Rmoon) * sin(θ_i), rtol=1e-3)
    # @assert isapprox(regolith_index(indexmodel, ice_depth) * (Rmoon - ice_depth) * sin(θice),
    #     regolith_index(indexmodel, 0.0cm) * (Rmoon) * sin(θ_i), rtol=1e-3)

    # and now finally account for the extra distance
    # "back up to" to Xmax after reflection
    Drego += (Drego - DXmax)

    # return the total length of the ray and
    # the zenith angle of the ray at Xmax and at ice_depth
    return Drego, θXmax, θice

end
