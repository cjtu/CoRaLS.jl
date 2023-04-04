using ProgressMeter

"""
A struct to hold the results of an acceptance calculation
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
end

"""
     acceptance(ntrials, nbins; ...)

Calculate the acceptance of CoRaLS to sub-surface UHECR reflections from ice.
"""
function acceptance(ntrials, nbins; min_energy=0.1EeV,
                    max_energy=600.0EeV,
                    altitude=20.0km,
                    trigger=magnitude_trigger(100μV/m),
                    save_events=true,
                    kwargs...)

    # create a StructArray for out Detection type
    # this is a vector with the length of the number of threads
    # we concat them together at the end to make a single collection of events
    devents = Vector{StructArray{Direct}}()
    revents = Vector{StructArray{Reflected}}()

    # setup the struct arrays for each thread
    for i in 1:Threads.nthreads()
        push!(devents, StructArray{Direct}(undef, 0))
        push!(revents, StructArray{Reflected}(undef, 0))
    end

    @info "Calculating acceptance using $(ntrials) trials across $(nbins) bins..."

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
    AΩ *= 30_000km^2  / spherical_cap_area(θpole, Rmoon)

    # and account for the fraction of a single hemisphere of the CRs
    # orbit that is visible from the poles, since we only simulate these
    # "visible" SC angles
    AΩ *= (θpole + θmax) / (π/2.0)

    # construct the bins we sample from in log-space
    energies = (10.0 .^ range(log10(min_energy / 1.0EeV), log10(max_energy / 1.0EeV), length=nbins+1))EeV

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
        Ecr = [sample_auger(energies[bin], energies[bin+1]) for i=1:ntrials]

        # loop over the number of trials in this bin
        Threads.@threads for i = 1:ntrials

            # throw a random cosmic ray trial and get the signal at the payload
            direct, reflected = throw_cosmicray(Ecr[i]; altitude=altitude, kwargs...)

            # check for a direct trigger
            if !(direct isa TrialFailed)
                if trigger(direct)
                    direct.triggered = true
                    Threads.atomic_add!(dpassed, 1)
                end
                save_events && push!(devents[Threads.threadid()], direct)
            end

            # check for a reflected trigger
            if !(reflected isa TrialFailed)
                if trigger(reflected)
                    reflected.triggered = true
                    Threads.atomic_add!(rpassed, 1)
                end
                save_events && push!(revents[Threads.threadid()], reflected)
            end

        end # end loop over trials

        # calculate the acceptance at this energy
        dAΩ[bin] *= dpassed[] / ntrials
        rAΩ[bin] *= rpassed[] / ntrials

    end # end loop over events

    # and combine the event structs across all threads
    devents = vcat(devents...)
    revents = vcat(revents...)

    return Acceptance(ntrials, altitude, AΩ, energies, dAΩ, devents, rAΩ, revents)

end


"""
    trials_passed

Calculate the number of successful trials in the direct and reflected acceptances.
"""
trials_passed(AΩ) = (AΩ.ntrials * (AΩ.dAΩ ./ AΩ.gAΩ), AΩ.ntrials * (AΩ.rAΩ ./ AΩ.gAΩ))


"""
    differential_spectrum(energies, AΩ, T)

Calculate the differential spectrum of detected UHECR events.

T is the total observation time of the experiment.
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
