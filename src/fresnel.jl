"""
    fresnel.jl

 Fresnel reflections and transmissions related to cosmic ray detection. It includes abstract types for divergence models and specific implementations for different divergence scenarios. It also contains functions for calculating Fresnel reflection and transmission coefficients and their modifications for spherical wave divergence at the surface.

## Main Components
- `DivergenceModel`: An abstract type representing a general model for wave divergence.

- `FarFieldDivergence`, `NearFieldDivergence`, `MixedFieldDivergence`: Specific types representing different divergence models based on field assumptions.

- `fresnel_coeffs`: Calculates Fresnel amplitude coefficients (reflected and transmitted, parallel and perpendicular).

- `divergence_tperp`, `divergence_tpar`: Functions that modify Fresnel transmission coefficients to include the effects of wave divergence at the surface.
"""

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
    snell_θt(θi, ni, nt)

    Calculate the refracted angle from Snell's Law: n_i * sin(θ_i) = n_t * sin(θ_t)
"""
snell_θt(θi, ni, nt) = asin((ni / nt) * sin(θi))

"""
    fresnel_critical(ni, nt)

Calculate the critical angle for total internal reflection.

# Arguments
- `ni`: Refractive index of the initial medium.
- `nt`: Refractive index of the transmitting medium.

# Returns
- Critical angle for total internal reflection.
"""
function fresnel_critical(ni, nt)
    ni > nt && return asin(nt / ni)
    return π / 2  # no total internal reflection
end

"""
    fresnel_coeffs(θi, ni, nt; simple_t = true)

Calculate the parallel and perpendicular Fresnel reflection and transmission 
coefficients. 

Each coefficient is the ratio of reflected/transmitted to incident amplitude
(denoted r and t, not to be confused with fresnel coefficients for power, R and T)

If simple_t is true (default) then the transmission coefficients
are calculated as follows, otherwise they use the full Fresnel equations.
    t_parallel = ni * (r_parallel + 1) / nt
    t_perpendicular = r_perpendicular + 1


# Arguments
- `θi`: Incident angle.
- `ni`: Refractive index of the initial medium.
- `nt`: Refractive index of the transmitting medium.

# Returns
- `(rpar, rperp, tpar, tperp)`: Fresnel coefficients.
"""
function fresnel_coeffs(θi, ni, nt; simple_t::Bool = true)
    if θi >= fresnel_critical(ni, nt)  # Total internal reflection
        return 1.0, 1.0, 0.0, 0.0
    elseif isapprox(θi, 0.0; atol = 1e-6)  # Prevent normal incidence => NaN
        return 0.0, 0.0, 1.0, 1.0
    end
    θt = snell_θt(θi, ni, nt)

    # Refelcted coefficients (Fresnel's sine and tangent law)
    rpar = tan(θi - θt) / tan(θi + θt)
    rperp = -sin(θi - θt) / sin(θi + θt)

    # Transmitted coefficients
    tpar = ni * (rpar + 1) / nt
    tperp = rperp + 1
    if !simple_t
        tpar = 2.0 * ni * cos(θi) / (nt * cos(θi) + ni * cos(θt))
        tperp = 2.0 * ni * cos(θi) / (ni * cos(θi) + nt * cos(θt))
    end

    return rpar, rperp, tpar, tperp
end

"""
    divergence_tperp(::FarFieldDivergence, θ_i, n, Drego, Dvacuum)

Calculate modified Fresnel transmission coefficients for the perpendicular polarization, considering far-field (spherical) divergence of waves at the surface.

# Arguments
- `θ_i`: Incident angle.
- `n`: Refractive index of the medium.
- `Drego`: Distance through regolith.
- `Dvacuum`: Distance through vacuum.

# Returns
- Modified Fresnel transmission coefficient for perpendicular polarization.

Note: We construct this so that it can be directly multiplied to an electric field,
so have to back out the total distance, before doing the scaling.

Note: Since the refraction occurs in the near field (Fresnel region), the standard
far-field (i.e. spherical) divergence assumption in the FORTE paper doesn't hold.
Our current model for this is an additional factor of `n` (which boosts the transmission).

See the FORTE paper, Appendix 2, for more details. arXiV:0309656
"""
function divergence_tperp(::FarFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin(n * sin(θ_i)) # the refracted angle
    return ((Drego + Dvacuum) / Dvacuum) * (2.0cos(θ_r) / (n * cos(θ_i) + cos(θ_r))) |> NoUnits
end

function divergence_tperp(::NearFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin(n * sin(θ_i)) # the refracted angle
    return n * ((Drego + Dvacuum) / Dvacuum) * (2.0cos(θ_r) / (n * cos(θ_i) + cos(θ_r))) |> NoUnits
end

function divergence_tperp(::MixedFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin(n * sin(θ_i)) # the refracted angle
    return sqrt(n) * ((Drego + Dvacuum) / Dvacuum) * (2.0cos(θ_r) / (n * cos(θ_i) + cos(θ_r))) |> NoUnits
end

"""
    divergence_tpar(θ_i, n, Drego, Dvacuum)

Calculate modified Fresnel transmission coefficients for the parallel polarization, considering far-field (spherical) divergence of waves at the surface.

        # Arguments
        - `θ_i`: Incident angle.
        - `n`: Refractive index of the medium.
        - `Drego`: Distance through regolith.
        - `Dvacuum`: Distance through vacuum.

Note: We construct this so that it can be directly multiplied to an electric field,
so have to back out the total distance, before doing the scaling.

Note: Since the refraction occurs in the near field (Fresnel region), the standard
far-field (i.e. spherical) divergence assumption in the FORTE paper doesn't hold.
Our current model for this is an additional factor of `n` (which boosts the transmission).

See the FORTE paper, Appendix 2, for more details. arXiV:0309656
"""
function divergence_tpar(::FarFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin(n * sin(θ_i)) # the refracted angle
    return ((Drego + Dvacuum) / Dvacuum) * (2.0cos(θ_r) / (n * cos(θ_r) + cos(θ_i))) |> NoUnits
end

function divergence_tpar(::NearFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin(n * sin(θ_i)) # the refracted angle
    return n * ((Drego + Dvacuum) / Dvacuum) * (2.0cos(θ_r) / (n * cos(θ_r) + cos(θ_i))) |> NoUnits
end

function divergence_tpar(::MixedFieldDivergence, θ_i, n, Drego, Dvacuum)
    θ_r = asin(n * sin(θ_i)) # the refracted angle
    return sqrt(n) * ((Drego + Dvacuum) / Dvacuum) * (2.0cos(θ_r) / (n * cos(θ_r) + cos(θ_i))) |> NoUnits
end
