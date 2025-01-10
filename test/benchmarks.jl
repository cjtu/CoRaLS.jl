"""
Timing the acceptance model.
"""

using BenchmarkTools
using CoRaLS
using Unitful: km, sr, EeV, MHz, m, μV
 

ntrials = 10000
nbins = 20
altitude = 50km
trigger=gaussian_trigger(sqrt(1000 / 300.)*0.5*67μV/m)
kws = Dict(:ice_depth=>10m,:trigger=>trigger,:min_energy=>10.0EeV, :max_energy=>600.0EeV,:ν_max=>1000MHz,:dν=>30MHz)

# [ Info: Calculating acceptance using 10000 trials across 20 bins...
# > 182.022 ms (1211287 allocations: 37.00 MiB)
@btime A = acceptance(ntrials, nbins, region=AllPSR(), spacecraft=CircularOrbit(altitude); kws...)

# [ Info: Calculating acceptance using 500 trials across 20 bins...
# > 1.860 s (12828195 allocations: 628.07 MiB)
@btime AO = old_acceptance(ntrials÷nbins, nbins; altitude=altitude, save_events=false, kws...)

@info "Finished benchmark."