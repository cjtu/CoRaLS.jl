using Profile
using CoRaLS
using Unitful: km, sr, EeV, MHz, m, μV, cm, ustrip

ntrials = 100
nbins = 3
altitude = 100.0km
trigger = LPDA(Nant=8, Ntrig=3, θ0=-90.0, altitude=altitude, skyfrac=0)
kws = Dict(:ice_depth=>10.0m,:trigger=>trigger,:min_energy=>5.0EeV,:max_energy=>50.0EeV,
           :ν_min=>150MHz,:ν_max=>1000MHz,:dν=>30MHz,:min_count=>50,:max_tries=>100)
region = create_region("polar:south,-80")
sc = CircularOrbit(altitude)

# Run from vs code
@profview acceptance(ntrials, nbins, region=region, spacecraft=sc; simple_area=true, kws...)
# @profview_allocs acceptance(ntrials, nbins, region=region, spacecraft=sc; simple_area=true, kws...)


# Run from julia repl
# @profview acceptance(ntrials, nbins, region=region, spacecraft=sc; simple_area=true, kws...)
# Profile.print()

# More granular profiling
# using JET
# Ecr = 10.0EeV
# ψ = 30.0
# Drego = 10m
# Dvacuum = 50.0km
# @report_opt regolith_field(ARW(), Ecr, ψ, Drego, Dvacuum)

# ν, E = regolith_field(ARW(), Ecr, ψ, Drego, Dvacuum)
# iceroughness = GaussianIceRoughness(0.05m)
# θ_ice = 5.0
# NXmax = 1.1
# # @report_call ice_roughness(iceroughness, ν, E, θ_ice, NXmax)