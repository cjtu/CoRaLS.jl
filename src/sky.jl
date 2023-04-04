
"""
   sky_temperature(ν)

Calculate the brightness temperature near the galactic poles in K.
"""
function sky_temperature(ν)

    # calculate the combined noise
    I_ν = extragalactic_noise(ν) + galactic_noise(ν)

    # and convert this to a temperature
    return ((I_ν / k_b) .* (c_0*c_0 ./ (2 .* ν .* ν))) .|> K

end

"""
    extragalactic_noise(ν)

Calculate the extragalactic contribution to RF background noise (in
W m^-2 Hz^-1 sr^-1) as function of frequency (in MHz) near the galactic poles.

These parameterization is taken from (Dulk, 2001):
https://www.aanda.org/articles/aa/full/2001/02/aads1858/aads1858.right.html
"""
function extragalactic_noise(ν)

    # convert ν into MHz and remove the units
    freqs = (ν ./ MHz) .|> NoUnits

    # the scaling factor for the galactic contribution.
    Ieg = 1.06e-20 * W / m^2 / Hz / sr

    # and the tau-factors
    tau = 5.0(freqs .^ -2.1)

    # and evaluate the galactic form
    return Ieg .* (freqs .^ -0.8) .* exp.(-tau)
end

"""
    galactic_noise(ν)

Calculate the galactic contribution to RF background noise (in
W m^-2 Hz^-1 sr^-1) as function of frequency (in MHz) near the galactic poles.

These parameterization is taken from (Dulk, 2001):
https://www.aanda.org/articles/aa/full/2001/02/aads1858/aads1858.right.html
"""
function galactic_noise(ν)

    # convert ν into MHz and remove the units
    freqs = (ν ./ MHz) .|> NoUnits

    # the scaling factor for the galactic contribution.
    Ig = 2.48e-20 * W / m^2 / Hz / sr

    # and the tau-factors
    tau = 5.0(freqs .^ -2.1)

    # and evaluate the galactic form
    return Ig .* (freqs .^ -0.52) .* ((1.0 .- exp.(-tau)) ./ tau)
end
