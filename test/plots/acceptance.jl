using Profile
using CoRaLS
using Unitful: km, sr, EeV, MHz, m, μV, cm, ustrip

ntrials = Int(1e5)  # Baseline test: 100_000
nbins = 20  # Baseline test: 20
ice_depth = 20.0m
trigger = magnitude_trigger(1e-5μV/m)
altitude = 50km

@profview A = acceptance(ntrials, nbins, region=SouthPolePSR, spacecraft=FixedPlatform(-90.0, 0.0, altitude), min_energy=1.0EeV, max_energy=600.0EeV)
