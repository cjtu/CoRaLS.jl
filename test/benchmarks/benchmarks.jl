"""
Timing the acceptance model.
"""

using Logging
using BenchmarkTools
using CoRaLS
using Unitful: km, sr, EeV, MHz, m, μV

disable_logging(Logging.Info)  # Hide acceptance runtimes

# We modify Ntrials for each run to get similar number of accepted events (benchmarking precision, not total just # events thrown)
ntrials = 10000
nbins = 8
altitude = 50.0km
trigger = LPDA(Nant=8, Ntrig=3, θ0=-85.0, altitude=altitude, skyfrac=0)
kws = Dict(:ice_depth=>6.0m,:trigger=>trigger,:min_energy=>1.0EeV,:max_energy=>100.0EeV,
           :ν_min=>150MHz,:ν_max=>1000MHz,:dν=>30MHz)
region = create_region("psr:south")
sc = FixedPlatform(-90, 0, altitude)

println("Starting benchmark with ", Threads.nthreads(), " threads...")
trial1 = () -> acceptance(ntrials*100, nbins, region=region, spacecraft=sc; kws...)
@btime trial1()

trial2 = () -> acceptance(ntrials, nbins, region=region, spacecraft=sc; simple_area=true, kws...)
@btime trial2()

trial3 = () -> old_acceptance(ntrials, nbins; altitude=altitude, save_events=false, kws...)
@btime trial3()

println("Finished benchmark")

# Showing trials are comparable with # reflected events
println("N reflected:")
println("Anew accepted ≈ ", sum(trial1().rcount))
println("Anew (simplearea) accepted ≈ ", sum(trial2().rcount))
a_old = trial3()
println("Aold accepted ≈ ", Int(round(a_old.ntrials*sum(a_old.rAΩ) / a_old.gAΩ[1])))

# Last benchmark: June 11, 2025
# Starting benchmark with 8 threads...
# 1.196 s (2452334 allocations: 184.25 MiB)
# 219.253 ms (2529992 allocations: 130.51 MiB)
# Anew accepted ≈ 24
# Anew (simplearea) accepted ≈ 625
# Aold accepted ≈ 264