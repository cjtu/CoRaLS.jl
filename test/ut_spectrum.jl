using Test
using Unitful: eV, km, sr, yr, ustrip

include("../src/spectrum.jl")

"""
Use python code from "The energy spectrum" Jupyter tutorial on the Auger open 
data website https://www.auger.org/opendata/ (Accessed 07/2024). Uses data
from Auger Open Data summary.zip (https://doi.org/10.5281/zenodo.10488964)
"""
# Raw counts from Fig. 7 of Auger 2020 paper
raw_count = [83143, 47500, 28657, 17843, 12435, 8715, 6050, 4111, 2620, 1691, 991, 624, 372, 156, 83, 24, 9, 6, 0, 0]

# Procedure from the Auger open data notebook
log_E_min = 0.4
E_bins = 20
E_bin_size = 0.1
log_E_max = log_E_min + E_bins * E_bin_size
log_bins = range(log_E_min, log_E_max, length = E_bins + 1)
log_bin_centers = log_bins[1:end-1] .+ 0.05
bins = 10 .^ log_bins .* 1e18
bin_energy = 10 .^ log_bin_centers .* 1e18
bin_width = diff(bins)


@testset verbose = true "spectrum.jl" begin

    @testset "Test auger_spectrum_2020 flux" begin
        # Compute flux from the Auger spectrum
        # Note: the Auger data is Jraw, but the spectral fit computes J (which we want)
        # Technically we're comparing two different things but works out to same
        # shape and < 5% differnce so parameterization looks good
        J_bins = auger_spectrum_2020.(bins .* eV)  # [1 / (km^2 sr yr eV)]

        # Apply effective sampling area-time and trapezoidal integration over energy to get counts
        count_rate = bin_width * eV .* (J_bins[1:end-1] .+ J_bins[2:end]) ./ 2 # [1 / (km^2 sr yr)]
        count = count_rate * 60400 * km^2 * sr * yr  # Auger 2020 paper

        # 5% tolerance, some offset probably due to Jraw -> J conversion
        @test count ≈ raw_count rtol = 0.05  

    end

end