#!/usr/bin/env julia


import Pkg
# Activate the CoRaLS.jl project
Pkg.activate("CoRaLS.jl")    # adjust if your script lives elsewhere

using CoRaLS
using Unitful: km, m, sr, EeV, MHz
using Dates

#–– Parse command-line arguments ––#
if length(ARGS) != 3
    println(stderr, "Usage: julia acceptance.jl <altitude_km> <ice_depth_m>")
    exit(1)
end


# convert to Unitful quantities
altitude = parse(Float64, ARGS[1])km
ice_depth = parse(Float64, ARGS[2])m
index = parse(Float64, ARGS[3])
ENERGY1 = 0.5 * index * EeV
ENERGY2 = 0.5 * (index + 1) * EeV

#–– Set up your run ––#
region   = create_region("polar:south,-80,0.0557")
sc       = CircularOrbit(altitude)
trigger  = LPDA(Nant=4, Ntrig=4, θ0=-90.0, altitude=altitude, skyfrac=0)
kws      = Dict(
    :min_energy => ENERGY1, #30.0EeV,
    :max_energy => ENERGY2,#31.0EeV,
    :ν_min      => 300MHz,
    :ν_max      => 1000MHz,
    :dν         => 30MHz,
    :min_count  => 50,
    :max_tries  => 100,
    :simple_area=> true,
)
ntrials = 10000000
nbins   = 1


#–– Run acceptance ––#
A = acceptance(ntrials, nbins;
    region=region,
    spacecraft=sc,
    trigger=trigger,
    ice_depth=ice_depth,
    kws...,
)

println("FINISHED RUNNING THE ACCEPTANCE")

#–– Define MCSE ––#
function mcse(count::Int, n::Int)
    p = count / n
    return sqrt(p * (1 - p) / n)
end

#–– Compute spectra + errors ––#
r_spectra = differential_spectrum(A.energies, A.rAΩ, 1yr)
r_error   = differential_spectrum(A.energies, mcse.(A.rcount, A.ntrials) .* A.gAΩ, 1yr)

d_spectra = differential_spectrum(A.energies, A.dAΩ, 1yr)
d_error   = differential_spectrum(A.energies, mcse.(A.dcount, A.ntrials) .* A.gAΩ, 1yr)

#–– Print results ––#
timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
println("# Run on $timestamp")
println("# altitude = $(altitude), ice_depth = $(ice_depth)")
println("r_spectra = ", r_spectra)
println("r_error   = ", r_error)
println("d_spectra = ", d_spectra)
println("d_error   = ", d_error)

