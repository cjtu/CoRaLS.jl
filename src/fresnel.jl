using ToggleableAsserts

"""
    fresnel.jl

 Fresnel reflections and transmissions related to cosmic ray detection. It includes abstract types for divergence models and specific implementations for different divergence scenarios. It also contains functions for calculating Fresnel reflection and transmission coefficients and their modifications for spherical wave divergence at the surface.

## Main Components
- `DivergenceModel`: An abstract type representing a general model for wave divergence.

- `FarFieldDivergence`, `NearFieldDivergence`, `MixedFieldDivergence`: Specific types representing different divergence models based on field assumptions.

- `fresnel_rpar`, `fresnel_rperp`: Functions to calculate parallel and perpendicular Fresnel reflection coefficients.

- `fresnel_tpar`, `fresnel_tperp`: Functions to calculate parallel and perpendicular Fresnel transmission coefficients.

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
    fresnel_rpar(thetai, ni, nt)

Calculate the parallel Fresnel reflection coefficient. This coefficient represents the ratio of reflected to incident electric field amplitudes for parallel (p-polarized) light at the interface between two media.

# Arguments
- `thetai`: Incident angle.
- `ni`: Refractive index of the initial medium.
- `nt`: Refractive index of the transmitting medium.

# Returns
- Parallel Fresnel reflection coefficient.
"""
function fresnel_rpar(thetai, ni, nt)
    thetat = asin((ni / nt) * sin(thetai)) # transmitted angle
    rpar = tan(thetai - thetat) / tan(thetai + thetat)

    # check if valid
    isnan(rpar) && return 0.0

    @toggled_assert abs(rpar) <= 1.0 "Rpar > 1.0", rpar
    return rpar

end

"""
    fresnel_rperp(thetai, ni, nt)

Calculate the perpendicular Fresnel reflection coefficient. Similar to `fresnel_rpar`, but for perpendicular (s-polarized) light.

# Arguments
- `thetai`: Incident angle.
- `ni`: Refractive index of the initial medium.
- `nt`: Refractive index of the transmitting medium.

# Returns
- Perpendicular Fresnel reflection coefficient.
"""
function fresnel_rperp(thetai, ni, nt)
    thetat = asin((ni / nt) * sin(thetai)) # transmitted angle
    rperp = -sin(thetai - thetat) / sin(thetai + thetat)

    # check if valid
    isnan(rperp) && return 0.0

    @toggled_assert abs(rperp) <= 1.0 "Rperp > 1.0"
    return rperp

end

"""
    fresnel_tpar(thetai, ni, nt)

Calculate the parallel Fresnel transmission coefficient. It quantifies the amplitude of the transmitted electric field compared to the incident field for parallel polarization.

# Arguments
- `thetai`: Incident angle.
- `ni`: Refractive index of the initial medium.
- `nt`: Refractive index of the transmitting medium.

# Returns
- Parallel Fresnel transmission coefficient.
"""
function fresnel_tpar(thetai, ni, nt)
    return 1.0 - fresnel_rpar(thetai, ni, nt)
    # otherwise, go ahead with the Fresnel calculation
    thetat = asin((ni / nt) * sin(thetai)) # transmitted angle
    tpar = 2.0 * sin(thetat) * cos(thetai) /
           (sin(thetai + thetat) * cos(thetai - thetat))
    return tpar

end

"""
    fresnel_tperp(thetai, ni, nt)

Calculate the perpendicular Fresnel transmission coefficient for s-polarized light at an interface.

# Arguments
- `thetai`: Incident angle.
- `ni`: Refractive index of the initial medium.
- `nt`: Refractive index of the transmitting medium.

# Returns
- Perpendicular Fresnel transmission coefficient.
"""
function fresnel_tperp(thetai, ni, nt)
    return 1.0 - fresnel_rperp(thetai, ni, nt)
    # otherwise, go ahead with the Fresnel calculation
    thetat = asin((ni / nt) * sin(thetai)) # transmitted angle
    tperp = 2.0 * sin(thetat) * cos(thetai) / sin(thetai + thetat)
    return tperp
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
