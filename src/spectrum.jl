using DelimitedFiles
using Interpolations
using Unitful: yr

# @unit yr "year" Year 360.2425*24*60*60s false

"""
    sample_power_law(gamma; min_value, max_value)

Draw random samples from a power law with exponent `gammma`.

If provided, `min_value` and `max_value` set the range over which
the power-law is sampled (defaults to 1 EeV and 1000 EeV)

This uses the prescription from:
    https://mathworld.wolfram.com/RandomNumber.html
to convert a uniform deviate into a power law with a given exponent.

"""
function sample_power_law(gamma; min_value=1.0EeV, max_value=1000.0EeV)

    # get a uniform random sample on [0, 1]
    y = rand()

    # we need gamma+1. multiple times
    g1 = gamma+1

    return ((max_value^g1 - min_value^g1)*y + min_value^g1)^(1. / g1)

end

"""
    sample_power_law(gamma, N; min_value, max_value)

Draw `N` random samples from a power law with exponent `gammma`.

If provided, `min_value` and `max_value` set the range over which
the power-law is sampled (defaults to 1 EeV and 1000 EeV)

This uses the prescription from:
    https://mathworld.wolfram.com/RandomNumber.html
to convert a uniform deviate into a power law with a given exponent.

"""
function sample_power_law(gamma, N; min_value=1.0EeV, max_value=1000.0EeV)

    # get a uniform random sample on [0, 1]
    y = rand(N)

    # we need gamma+1. multiple times
    g1 = gamma+1

    return ((max_value^g1 - min_value^g1).*y .+ min_value^g1).^(1. / g1)

end

"""
Evaluate the Auger UHECR spectrum in eV^-1 km^-2 sr^-1 yr^-1.

This uses the parameterization from:
    https://academic.oup.com/ptep/article/2017/12/12A103/4665686

"""
function auger_spectrum(E)

    # J0 = 3.30 \pm 0.15 \pm 0.20 x 10^-19 eV^-1 km^-2 sr^-1 yr^-1
    J0 = 3.30e-19/eV/km^2/sr/yr

    # \gamma_1 = 3.29 \pm 0.02 \pm 0.05
    γ1 = 3.29
    # \gamma_2 = 2.60 \pm 0.02 \pm 0.1
    γ2 = 2.60
    # \Delta\gamma = 3.1 \pm 0.2 \pm 0.4
    Δγ = 3.1

    # the ankle energy, 4.82 \pm 0.07 \pm 0.8
    Eankle = 4.82EeV

    # the cutoff parameter, 42.1 \pm 1.7 \pm 7.6
    Es = 42.1EeV

    # if we are less than the ankle, it's just a power law
    if E <= Eankle
        return J0 * (E / Eankle)^(-γ1)
    else
        powerlaw = (E / Eankle)^(-γ2)
        boost = 1.0 + ( Eankle / Es)^(Δγ)
        cutoff = (1.0 + (E / Es)^(Δγ))^(-1.0)
        return J0 * powerlaw * boost * cutoff
    end
end

"""
Randomly evaluate the Auger UHECR spectrum in eV^-1 km^-2 sr^-1 yr^-1
within experimental uncertainties.

This method differs from the one argument version of auger_spectrum()
in that it randomly samples the parameters within experimental uncertainties
so that you can capture the variation in the cosmic ray spectrum.

This uses the parameterization from:
    https://academic.oup.com/ptep/article/2017/12/12A103/4665686

"""
function auger_spectrum(E, sample)

    # a quick function to generate random normal uncertainties
    uncert(σ) = rand(Normal(0, σ))

    # J0 = 3.30 \pm 0.15 \pm 0.20 x 10^-19 eV^-1 km^-2 sr^-1 yr^-1
    J0 = 3.30e-19/eV/km^2/sr/yr
    # J0 = (3.30 + uncert(0.15) + uncert(0.20))*1e-19/eV/km^2/sr/yr

    # \gamma_1 = 3.29 \pm 0.02 \pm 0.05
    γ1 = 3.29 + uncert(0.02) + uncert(0.05)
    # \gamma_2 = 2.60 \pm 0.02 \pm 0.1
    γ2 = 2.60 + uncert(0.02) + uncert(0.1)
    # \Delta\gamma = 3.1 \pm 0.2 \pm 0.4
    Δγ = 3.1 + uncert(0.2) + uncert(0.4)

    # the ankle energy, 4.82 \pm 0.07 \pm 0.8
    Eankle = (4.82 + uncert(0.07) + uncert(0.8))EeV
    Eankle = 4.82EeV

    # the cutoff parameter, 42.1 \pm 1.7 \pm 7.6
    Es = (42.1 + uncert(1.7) + uncert(7.6))EeV

    # construct the array to store the energies
    spectrum = zeros(length(E)) ./ km ./ km ./ sr ./ yr ./ eV

    # if we are less than the ankle, it's just a power law
    spectrum[E .<= Eankle] = J0 .* (E[E .<= Eankle] ./ Eankle).^(-γ1)

    # otherwise, build up the components and add them to our array
    powerlaw = (E[E .> Eankle] ./ Eankle).^(-γ2)
    boost = 1.0 .+ ( Eankle ./ Es).^(Δγ)
    cutoff = (1.0 .+ (E[E .> Eankle] ./ Es).^(Δγ)).^(-1.0)
    spectrum[E .> Eankle] = J0 .* powerlaw .* boost .* cutoff

    return spectrum
end

let ldata=log10.(readdlm("$(@__DIR__)/../data/CR_FLUX_1e9_1e18.csv", ',')),
    interp = interpolate((ldata[:, 1],), ldata[:, 2], Gridded(Linear()))
    """
    Evaluate the cosmic ray spectrum at lower energies.
    """
    global function cr_spectrum(E)

        # calculate the log10(E) as we interpolate in log-space in eV
        lE = log10.((E ./ eV) .|> NoUnits)

        # interpolate into the file to get a log-flux in log(m^2 sr^1 s^1 GeV^1)-1
        lflux = interp.(lE)

        # and return the flux in sensible units
        return (10.0 .^ lflux) ./ m^2 ./ s ./ sr ./ eV

    end
end

"""
    sample_auger(min_energy, max_energy)

Draw a random sample from the Auger UHECR flux between
`min_energy` and `max_energy`.
"""
function sample_auger(min_energy, max_energy)

    # get the min and max flux values for the rejection sampling
    # assumes a falling spectrum
    min_flux = auger_spectrum(max_energy)
    max_flux = auger_spectrum(min_energy)

    # draw a random sample in EeV
    E = rand(Uniform(min_energy ./ EeV, max_energy ./ EeV))EeV

    # loop until we get a valid sample
    while (max_flux - min_flux)*rand() + min_flux > auger_spectrum(E)
        E = rand(Uniform(min_energy ./ EeV, max_energy ./ EeV))EeV
    end

    return E

end
