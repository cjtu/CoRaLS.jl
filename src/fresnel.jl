using ToggleableAsserts

"""
Abstract type for a divergence model for surface refraction.
"""
abstract type DivergenceModel end

"""
Far-field (spherical) divergence.

This is taken from the FORTE paper.
"""
struct FarFieldDivergence <: DivergenceModel end

"""
Near-field (cylindrical) divergence.

This is currently approximated by an additional `n` factor
out front of the far-field (spherical) divergence.
"""
struct NearFieldDivergence <: DivergenceModel end

"""
Somewhere between near- and far- field divergence.

This is currently approximated by an additional `sqrt(n)` factor
out front of the far-field (spherical) divergence.
"""
struct MixedFieldDivergence <: DivergenceModel end

"""
    fresnel_rpar(thetai, ni, nt)

Calculate the parallel Fresnel reflection coefficient.
"""
function fresnel_rpar(thetai, ni, nt)
    thetat = asin( (ni / nt) * sin(thetai) ) # transmitted angle
    rpar = tan(thetai - thetat) / tan(thetai + thetat)

    # check if valid
    isnan(rpar) && return 0.0

    @toggled_assert abs(rpar) <= 1.0 "Rpar > 1.0", rpar
    return rpar

end

"""
    fresnel_rperp(thetai, ni, nt)

Calculate the perpendicular Fresnel reflection coefficient.
"""
function fresnel_rperp(thetai, ni, nt)
    thetat = asin( (ni / nt) * sin(thetai) ) # transmitted angle
    rperp = -sin(thetai - thetat) / sin(thetai + thetat)

    # check if valid
    isnan(rperp) && return 0.0

    @toggled_assert abs(rperp) <= 1.0 "Rperp > 1.0"
    return rperp

end

"""
    fresnel_tpar(thetai, ni, nt)

Calculate the parallel Fresnel transmission coefficient.
"""
function fresnel_tpar(thetai, ni, nt)
    return 1.0 - fresnel_rpar(thetai, ni, nt)
    # otherwise, go ahead with the Fresnel calculation
    thetat = asin( (ni / nt) * sin(thetai) ) # transmitted angle
    tpar = 2.0 * sin(thetat) * cos(thetai) /
        (sin(thetai + thetat) * cos(thetai - thetat))
    return tpar

end

"""
    fresnel_tperp(thetai, ni, nt)

Calculate the perpendicular Fresnel transmission coefficient.
"""
function fresnel_tperp(thetai, ni, nt)
    return 1.0 - fresnel_rperp(thetai, ni, nt)
    # otherwise, go ahead with the Fresnel calculation
    thetat = asin( (ni / nt) * sin(thetai) ) # transmitted angle
    tperp = 2.0 * sin(thetat) * cos(thetai) / sin(thetai + thetat)
    return tperp
end

"""
    divergence_tperp(θ_i, n, Drego, Dvacuum)

Modified Fresnel transmission coefficients, for the perpendicular
polarization, to include the divergence of spherical waves at the surface.

Note: We construct this so that it can be directly multiplied to an electric field,
so have to back out the total distance, before doing the scaling.

See the FORTE paper, Appendix 2, for more details. arXiV:0309656
"""
function divergence_tperp(::FarFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin( n * sin(θ_i) ) # the refracted angle
    return ((Drego + Dvacuum) / Dvacuum) * ( 2.0cos(θ_r) / (n*cos(θ_i) + cos(θ_r)) ) |> NoUnits
end

"""
    divergence_tpar(θ_i, n, Drego, Dvacuum)

Modified Fresnel transmission coefficients, for the parallel
polarization, to include the divergence of spherical waves at the surface.

Note: We construct this so that it can be directly multiplied to an electric field,
so have to back out the total distance, before doing the scaling.

See the FORTE paper, Appendix 2, for more details. arXiV:0309656
"""
function divergence_tpar(::FarFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin( n * sin(θ_i) ) # the refracted angle
    return ((Drego + Dvacuum) / Dvacuum) * ( 2.0cos(θ_r) / (n*cos(θ_r) + cos(θ_i)) ) |> NoUnits
end


"""
    divergence_tperp(::NearFieldDivergence, θ_i,  n, Drego, Dvacuum)

Modified Fresnel transmission coefficients, for the perpendicular
polarization, to include the divergence of spherical waves at the surface.

Note: We construct this so that it can be directly multiplied to an electric field,
so have to back out the total distance, before doing the scaling.

Note: Since the refraction occurs in the near field (Fresnel region), the standard
far-field (i.e. spherical) divergence assumption in the FORTE paper doesn't hold.
Our current model for this is an additional factor of `n` (which boosts the transmission).

See the FORTE paper, Appendix 2, for more details. arXiV:0309656
"""
function divergence_tperp(::NearFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin( n * sin(θ_i) ) # the refracted angle
    return n * ((Drego + Dvacuum) / Dvacuum) * ( 2.0cos(θ_r) / (n*cos(θ_i) + cos(θ_r)) ) |> NoUnits
end

"""
    divergence_tpar(::NearFieldDivergence, θ_i, n, Drego, Dvacuum)

Modified Fresnel transmission coefficients, for the parallel
polarization, to include the divergence of spherical waves at the surface.

Note: We construct this so that it can be directly multiplied to an electric field,
so have to back out the total distance, before doing the scaling.

Note: Since the refraction occurs in the near field (Fresnel region), the standard
far-field (i.e. spherical) divergence assumption in the FORTE paper doesn't hold.
Our current model for this is an additional factor of `n` (which boosts the transmission).

See the FORTE paper, Appendix 2, for more details. arXiV:0309656
"""
function divergence_tpar(::NearFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin( n * sin(θ_i) ) # the refracted angle
    return n * ((Drego + Dvacuum) / Dvacuum) * ( 2.0cos(θ_r) / (n*cos(θ_r) + cos(θ_i)) ) |> NoUnits
end

"""
    divergence_tperp(::MixedFieldDivergence, θ_i, n, Drego, Dvacuum)

Modified Fresnel transmission coefficients, for the perpendicular
polarization, to include the divergence of spherical waves at the surface.

Note: We construct this so that it can be directly multiplied to an electric field,
so have to back out the total distance, before doing the scaling.

Note: Since the refraction occurs in the near field (Fresnel region), the standard
far-field (i.e. spherical) divergence assumption in the FORTE paper doesn't hold.
Our current model for this is an additional factor of `n` (which boosts the transmission).

See the FORTE paper, Appendix 2, for more details. arXiV:0309656
"""
function divergence_tperp(::MixedFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin( n * sin(θ_i) ) # the refracted angle
    return sqrt(n) * ((Drego + Dvacuum) / Dvacuum) * ( 2.0cos(θ_r) / (n*cos(θ_i) + cos(θ_r)) ) |> NoUnits
end

"""
    divergence_tpar(::MixedFieldDivergence, θ_i, n, Drego, Dvacuum)

Modified Fresnel transmission coefficients, for the parallel
polarization, to include the divergence of the RF at the surface.

See the comments on the other methods for more details.
"""
function divergence_tpar(::MixedFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin( n * sin(θ_i) ) # the refracted angle
    return sqrt(n) * ((Drego + Dvacuum) / Dvacuum) * ( 2.0cos(θ_r) / (n*cos(θ_r) + cos(θ_i)) ) |> NoUnits
end
