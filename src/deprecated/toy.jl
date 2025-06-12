""" 
Note (2025/02/02): This is no longer applicable since shifting to new orbital 
and regional sampling strategy (see geometry.jl).
"""

# This file implements several *toy* MCs to evaluate the geometric
# acceptance of various different sampling strategies for CoRaLS
# i.e. the original "non-orbiting" one against newer "orbital" strategies

"""
    nonorbiting(altitude)

Calculate the geometric acceptance using the original "fixed" location algorithm proposed by Andres and used by Remy until Sep. 22 2021.

# Arguments
- `altitude`: Altitude of observation or detection.

# Returns
Geometric acceptance for the non-orbiting model.

# Description
This function calculates the geometric acceptance using the initial fixed location algorithm. It computes the acceptance based on a fixed spherical cap and accounts for the total visible area from a given altitude.
"""
function nonorbiting(altitude)

    # the maximum angle visible from a given altitude
    θmax = -horizon_angle(altitude)

    # we use a fixed spherical cap to sample CRs from
    A = spherical_cap_area(θmax, Rmoon)

    # the corresponding acceptance
    AΩ = π * sr * A

    # the total area of the orbit is a "stripe" that is θmax wide
    Atotal = 2.0 * (spherical_cap_area(π / 2.0) -
                    spherical_cap_area(π / 2.0 + horizon_angle(altitude)))

    # for all these comparisons, we fix the PSR area at 30,000km^2
    AΩ *= (30_000km^2 / Atotal)

    return AΩ

end

"""
    orbiting(altitude; ntrials=2e7)

Calculate the geometric acceptance for the "orbiting" sampling algorithm proposed by Remy in late September 2021.

# Arguments
- `altitude`: Altitude of observation or detection.
- `ntrials`: Number of trials for the Monte Carlo simulation (default: 2e7).

# Returns
Geometric acceptance for the orbiting model.

# Description
This function calculates the geometric acceptance using a newer "orbiting" sampling algorithm proposed by Remy. It factors in the spacecraft's orbital movement and visibility constraints relative to lunar poles. The function conducts a Monte Carlo simulation to estimate the acceptance, considering various trial scenarios.
"""
function orbiting(altitude, ntrials=2e7)

    θmax = -horizon_angle(altitude)
    θpole = deg2rad(10.0)

    A = spherical_cap_area(π / 2.0, Rmoon)
    # A = spherical_cap_area(π, Rmoon)

    # 0.5 to account for the SC being in the *other* hemisphere
    # for half of the time
    AΩ = 0.5 * π * sr * A
    # AΩ = π * sr * A

    # only allowing points within the top spherical cap
    AΩ *= (30_000km^2 / (1.0 * spherical_cap_area(θpole)))

    # points can come from both spherical caps
    # AΩ *= (30_000km^2  / (2.0*spherical_cap_area(θpole)))

    passed = 0

    for i = 1:ntrials

        # # choose point on top hemisphere
        surface = random_point_on_cap(π / 2.0; r=Rmoon)

        # choose a point on the moon
        # surface = random_point_on_cap(π; r=Rmoon)

        # calculate the sampled value of θ for the surface point
        θ, _, _ = CoRaLS.cartesian_to_spherical(surface)

        # check if the point lies outside the polar cap
        ((θ > θpole) && θ < (π - θpole)) && continue

        # choose a SC location on the moon
        # θsc = rand(Uniform(0, π))

        # choose a SC location in the top hemisphere
        θsc = rand(Uniform(0, π / 2.0))

        # if the SC is not in view of the lunar pols
        ((θsc > 2.0 * θpole) && θsc < (π - 2.0 * θpole)) && continue

        # calculate the vector location of the SC
        SC = CoRaLS.spherical_to_cartesian(θsc, rand(Uniform(0, 2π)), Rmoon + altitude)

        # and the great circle angle between the SC and the surface point
        Δσ = atan(norm((SC / norm(SC)) × surface), ((SC / norm(SC)) ⋅ surface))

        # if the point is beyond the visible horizon, then we can't ever see it!
        abs(Δσ) >= θmax && continue

        # if we get here, this is a valid trial
        passed += 1

    end

    return AΩ * (passed / ntrials)

end
