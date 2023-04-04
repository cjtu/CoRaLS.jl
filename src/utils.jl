
"""
    rician(ν, σ)

Draw a Rician random variable.
"""
rician(ν::Real, σ::Real) = sqrt( rand(Normal(ν, σ))^2 + rand(Normal(0., σ))^2 )

"""
    rician(ν, σ, N)

Draw `N` Rician random variables where ν and σ are constants.
"""
rician(ν::Real, σ::Real, N) = sqrt.( rand(Normal(ν, σ), N).^2 .+ rand(Normal(0., σ), N).^2 )

"""
    rician(ν, σ, N)

Draw Rician random variables where ν and σ are arrays.
"""
rician(ν::AbstractVector{<:Real},
       σ::AbstractVector{<:Real}) = sqrt.( rand.(Normal.(ν, σ)).^2. .+ rand.(Normal.(zeros(length(ν)), σ)).^2. )

"""
    polarization_angle(events)

Calculate the polarization, in degrees, for a given set of events
"""
function polarization_angle(events)
    # construct the Poynting vector for each trial
    S = -spherical_to_cartesian.(π/2.0 .+ abs.(deg2rad.(events.θ_el)),
                                 deg2rad.(events.ϕ), 1)

    # calculate the total project of the polarization vector
    # onto the plane perpendicular to S - should `proj`
    # always have the same magnitude as `pol` even in the case
    # of surface roughness? Right now, most events have norm(pol) == norm(proj)
    # but there are rare events that can go down to 0.5.
    proj = events.pol .- (events.pol .⋅ S).*S

    # define our positive x-axis by pointing to the "right"
    # in the antennas frame; this is H-Pol!
    Haxis = [[0, 0, 1] × S[i] for i=1:length(S)]
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

Recalculate the effective acceptance for direct and reflected events
given a set of events and a trigger and return the new effective
acceptances.
"""
function retrigger(AΩ, trig)

    # check that both collection of events are not-empty
    @assert length(AΩ.reflected)>0 "AΩ must have been calculated with save_events=true"
    @assert length(AΩ.direct)>0 "AΩ must have been calculated with save_events=true"

    # an array to store our trigger pass rates
    Dpassed = zeros(length(AΩ.energies)-1)
    Rpassed = zeros(length(AΩ.energies)-1)

    # loop over each bin
    Threads.@threads for i = 1:(length(AΩ.energies)-1)

        # find all of the events with energies in this bin
        Dinbin = (AΩ.direct.Ecr .>= AΩ.energies[i]) .& (AΩ.direct.Ecr .<= AΩ.energies[i+1])
        Rinbin = (AΩ.reflected.Ecr .>= AΩ.energies[i]) .& (AΩ.reflected.Ecr .<= AΩ.energies[i+1])

        # apply the trigger to all these events to calculate those that passed
        Dpassed[i] = sum(trig.(AΩ.direct[Dinbin]))
        Rpassed[i] = sum(trig.(AΩ.reflected[Rinbin]))

    end

    # convert the trigger efficieny into acceptance
    dAΩ = AΩ.gAΩ .* (Dpassed ./ AΩ.ntrials)
    rAΩ = AΩ.gAΩ .* (Rpassed ./ AΩ.ntrials)

    return dAΩ, rAΩ

end
