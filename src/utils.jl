"""
    rician(ν, σ, N = 1)

Draw one or more Rician random variables.

Arguments:
- `ν`: The noncentrality parameter, can be a scalar or an array.
- `σ`: The scale parameter, can be a scalar or an array matching the size of `ν`.
- `N = 1`: The number of samples to draw. Used only if `ν` is a scalar.

# Returns
- If `ν` and `σ` are single values, it returns a single Rician random variable.
- If `ν` and `σ` are arrays, it returns an array of Rician random variables with each element corresponding to elements in `ν` and `σ`.
- If `ν` is a single value and `N > 1`, it returns an array of `N` Rician random variables, all generated with the same `ν` and `σ`.

"""
function rician(ν, σ, N=1)
    if N > 1 || length(ν) > 1 || length(σ) > 1
        # Ensure ν and σ are arrays of matching size
        ν_arr = (length(ν) == 1) ? fill(ν, N) : ν
        σ_arr = (length(σ) == 1) ? fill(σ, length(ν_arr)) : σ
        return sqrt.(rand.(Normal.(ν_arr, σ_arr)) .^ 2 .+ rand.(Normal.(zeros(length(ν_arr)), σ_arr)) .^ 2)
    else
        # Single scalar case
        return sqrt(rand(Normal(ν, σ))^2 + rand(Normal(0.0, σ))^2)
    end
end

"""
    polarization_angle(events)

Calculate the polarization angle for a given set of events.

# Arguments
- `events`: A collection of event data, each with a polarization vector.

# Returns
An array of polarization angles in degrees.

# Description
This function computes the polarization angle for each event in the provided collection. The polarization angle is calculated relative to the horizontal plane of the antenna array, with 0 degrees corresponding to horizontal polarization and 90 degrees to vertical polarization.
"""
function polarization_angle(events)
    # construct the Poynting vector for each trial
    S = -spherical_to_cartesian.(π / 2.0 .+ abs.(deg2rad.(events.θ_el)),
        deg2rad.(events.ϕ), 1)

    # calculate the total project of the polarization vector
    # onto the plane perpendicular to S - should `proj`
    # always have the same magnitude as `pol` even in the case
    # of surface roughness? Right now, most events have norm(pol) == norm(proj)
    # but there are rare events that can go down to 0.5.
    proj = events.pol .- (events.pol .⋅ S) .* S

    # define our positive x-axis by pointing to the "right"
    # in the antennas frame; this is H-Pol!
    Haxis = [[0, 0, 1] × S[i] for i = 1:length(S)]
    Haxis ./= norm.(Haxis) # make sure it's normalized

    # and now construct the V-Pol or y-axis by crossing
    # S against with Haxis - should point "up"
    Vaxis = S .× Haxis
    Vaxis ./= norm.(Vaxis) # make sure it's normalized

    # and now project the polarization vector onto both
    # of these basis vectors
    Hpol = events.pol .⋅ Haxis
    Vpol = events.pol .⋅ Vaxis

    # compute the polarization angle - 0 == H-Pol, 90 = V-Pol
    χ = atan.(Vpol, Hpol)

    return rad2deg.(χ)
end

"""
    retrigger(AΩ, trig)

Recalculate the effective acceptance for direct and reflected events.

# Arguments
- `AΩ`: A struct holding initial acceptance calculation results.
- `trig`: A trigger function used to reassess event acceptance.

# Returns
A tuple containing the recalculated effective acceptances for direct and reflected events.

# Description
This function applies a new trigger to a set of cosmic ray events, recalculating the effective acceptance based on which events would pass this new trigger. It's useful for assessing the impact of different trigger criteria on event detection rates.
"""
function retrigger(AΩ, trig)

    # check that both collection of events are not-empty
    @assert length(AΩ.reflected) > 0 "AΩ must have been calculated with save_events=true"
    @assert length(AΩ.direct) > 0 "AΩ must have been calculated with save_events=true"

    # an array to store our trigger pass rates
    Dpassed = zeros(length(AΩ.energies) - 1)
    Rpassed = zeros(length(AΩ.energies) - 1)

    # loop over each bin
    for i = 1:(length(AΩ.energies)-1)

        # find all of the events with energies in this bin
        Dinbin = (AΩ.direct.Ecr .>= AΩ.energies[i]) .& (AΩ.direct.Ecr .<= AΩ.energies[i+1])
        Rinbin = (AΩ.reflected.Ecr .>= AΩ.energies[i]) .& (AΩ.reflected.Ecr .<= AΩ.energies[i+1])

        # apply the trigger to all these events to calculate those that passed
        if sum(Dinbin) > 0
            Dpassed[i] = sum(trig.(AΩ.direct[Dinbin]))
        end
        if sum(Rinbin) > 0 
            Rpassed[i] = sum(trig.(AΩ.reflected[Rinbin]))
        end

    end

    # convert the trigger efficieny into acceptance
    dAΩ = AΩ.gAΩ .* (Dpassed ./ AΩ.ntrials)
    rAΩ = AΩ.gAΩ .* (Rpassed ./ AΩ.ntrials)

    AΩnew = deepcopy(AΩ)
    AΩnew.dcount = Dpassed
    AΩnew.rcount = Rpassed
    AΩnew.dAΩ = dAΩ
    AΩnew.rAΩ = rAΩ

    return AΩnew

end
