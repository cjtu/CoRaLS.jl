using ProgressMeter
"""
    acceptance.jl

This file is part of the `CoRaLS` module and contains functions and structures for calculating the acceptance of the system for sub-surface Ultra-High-Energy Cosmic Ray (UHECR) reflections from ice. It also includes utilities for analyzing trials and computing differential spectra.

## Main Components

- `struct Acceptance`: A structure to hold results of acceptance calculations, including trials, altitude, energies, and different acceptance types (direct and reflected).

- `function acceptance(...)`: Calculates the acceptance of CoRaLS for sub-surface UHECR reflections. It involves simulations of cosmic ray interactions, triggering conditions, and aggregating results across multiple trials and energy bins. 

- `function trials_passed(...)`: Computes the number of successful trials for both direct and reflected acceptances. 

- `function differential_spectrum(...)`: Calculates the differential spectrum of detected UHECR events based on acceptance calculations and observation time.
"""

"""
A structure to hold results of acceptance calculations, including trials, altitude, energies, and different acceptance types (direct and reflected).
"""
struct Acceptance
    ntrials::Int64 # per bin
    altitude
    gAΩ
    energies
    dAΩ
    direct
    rAΩ
    reflected
    failed_direct
    failed_reflected
end

"""
acceptance(ntrials, nbins; ...)

Calculates the acceptance of CoRaLS for sub-surface UHECR reflections. It involves simulations of cosmic ray interactions, triggering conditions, and aggregating results across multiple trials and energy bins.
"""
function acceptance(ntrials, nbins; min_energy=0.1EeV,
    max_energy=600.0EeV,
    altitude=20.0km,
    trigger=magnitude_trigger(100μV / m),
    save_events=true,
    orbit=false,
    target=AllPSR(),
    kwargs...)

    # create a StructArray for out Detection type
    # this is a vector with the length of the number of threads
    # we concat them together at the end to make a single collection of events
    nthreads = Threads.nthreads()
    devents = fill(StructArray{Direct}(undef, 0), nthreads)
    revents = fill(StructArray{Reflected}(undef, 0), nthreads)
    fdevents = fill(CoRaLS.TrialFailed[], nthreads)
    frevents = fill(CoRaLS.TrialFailed[], nthreads)

    @info "Calculating acceptance using $(ntrials) trials across $(nbins) bins..."

    if orbit
        # the maximum total acceptance is πA km^2 steradians (TODO: see paper X, equation Y)
        A = Moon().area  # [km^2] sampling events over whole Moon
        AΩ = π * sr * A * ones(nbins)
    else
    # we define the lunar poles as everything within 10 degrees
        θpole = deg2rad(10.0)

        # the distance to the horizon from the SC's altitude
        θmax = -horizon_angle(altitude)

        # we sample events only in the lunar poles
        A = spherical_cap_area(θpole, Rmoon)

        # the maximum total acceptance is πA km^2 steradians
        # however the spacecraft is in the *opposite* hemisphere from the
        # cosmic ray trial 50% of the time, so we pick up another factor of 0.5
        # since we are only simulating one "combined" pole (to be more efficient).
        AΩ = 0.5 * π * sr * A * ones(nbins)

        # fold in the fraction of the polar caps that have PSRs
        # assuming 30_000 km^2 of total PSR regions
        AΩ *= 30_000km^2 / spherical_cap_area(θpole, Rmoon)  # TODO: 30000 is both poles, divinding by only 1 sph cap 

        # and account for the fraction of a single hemisphere of the CRs
        # orbit that is visible from the poles, since we only simulate these
        # "visible" SC angles
        AΩ *= (θpole + θmax) / (π / 2.0)  # TODO:rederive this - max area visible by payload
    end
    # construct the bins we sample from in log-space
    energies = (10.0 .^ range(log10(min_energy / 1.0EeV), log10(max_energy / 1.0EeV), length=nbins + 1))EeV

    # construct arrays containing the effective area for direct
    # and reflected events at each energy
    dAΩ = AΩ .* ones(nbins)
    rAΩ = AΩ .* ones(nbins)

    # try the each number of event that was requested
    @showprogress 1 "Simulating..." for bin = 1:nbins

        # the number of trials that have passed for direct and reflected
        dpassed = Threads.Atomic{Int}(0.0)
        rpassed = Threads.Atomic{Int}(0.0)

        # sample the Auger spectrum for UHECR energies within this bin
        Ecr = [sample_auger(energies[bin], energies[bin+1]) for i = 1:ntrials]

        # sample orbital positions from file
        if orbit
            SCs = sample_orbit(parse_orbit(), ntrials)
        end

        # loop over the number of trials in this bin
        Threads.@threads for i = 1:ntrials

            # throw a random cosmic ray trial and get the signal at the payload
            if orbit && altitude <= 0km
                direct, reflected = throw_cosmicray(Ecr[i], SCs[i, :]; target=target, kwargs...)
            elseif orbit && altitude > 0km
                direct, reflected = throw_cosmicray(Ecr[i], altitude; target=target, kwargs...)
            else
                direct, reflected = throw_cosmicray(Ecr[i]; altitude=altitude, target=target, kwargs...)
            end

            # check for a direct trigger
            if !(direct isa TrialFailed)
                if trigger(direct)
                    direct.triggered = true
                    Threads.atomic_add!(dpassed, 1)
                end
                save_events && push!(devents[Threads.threadid()], direct)
            else
                save_events && push!(fdevents[Threads.threadid()], direct)
            end

            # check for a reflected trigger
            if !(reflected isa TrialFailed)
                if trigger(reflected)
                    reflected.triggered = true
                    Threads.atomic_add!(rpassed, 1)
                end
                save_events && push!(revents[Threads.threadid()], reflected)
            else
                save_events && push!(frevents[Threads.threadid()], reflected)
            end

        end # end loop over trials

        # calculate the acceptance at this energy
        dAΩ[bin] *= dpassed[] / ntrials
        rAΩ[bin] *= rpassed[] / ntrials

    end # end loop over events

    # and combine the event structs across all threads
    devents = vcat(devents...)
    revents = vcat(revents...)
    fdevents = vcat(fdevents...)
    frevents = vcat(frevents...)

    return Acceptance(ntrials, altitude, AΩ, energies, dAΩ, devents, rAΩ, revents, fdevents, frevents)

end


"""
    trials_passed

Calculate the number of successful trials in the direct and reflected acceptances.
"""
trials_passed(AΩ) = (AΩ.ntrials * (AΩ.dAΩ ./ AΩ.gAΩ), AΩ.ntrials * (AΩ.rAΩ ./ AΩ.gAΩ))


"""
    differential_spectrum(energies, AΩ, T)

Calculate the differential spectrum of detected UHECR events based on acceptance calculations and T, the total observation time of the experiment.

Calculates the differential spectrum of detected UHECR events .
"""
function differential_spectrum(energies, AΩ, T)

    # array to store the spectrum in each bin
    spectrum = zeros(length(AΩ))

    # loop over each bin
    for bin = 1:(length(energies)-1)

        # calculate the average of the auger spectrum in this bin
        avg = mean(auger_spectrum.(range(energies[bin], stop=energies[bin+1], length=100)))

        # calculate dN/DE using this average value of the spectrum
        # spectrum[bin] = (T * AΩ[bin] * avg) |> (EeV^(-1))
        spectrum[bin] = (T * AΩ[bin] * avg * (energies[bin+1] - energies[bin])) |> NoUnits

    end

    # and we are done
    return spectrum
end
