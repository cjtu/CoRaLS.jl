using CoRaLS
using PyPlot
using Unitful: eV, km, sr, yr, ustrip

function plot_auger_2020_comparison()
    # Raw counts from Fig. 7 of Auger 2020 paper
    raw_count = [83143, 47500, 28657, 17843, 12435, 8715, 6050, 4111, 2620, 1691, 991, 624, 372, 156, 83, 24, 9, 6, 0, 0]

    # Init arrays
    log_E_min = 0.4
    E_bins = 20
    E_bin_size = 0.1
    log_E_max = log_E_min + E_bins * E_bin_size
    log_bins = range(log_E_min, log_E_max, length = E_bins + 1)
    log_bin_centers = log_bins[1:end-1] .+ 0.05
    bins = 10 .^ log_bins .* 1e18
    bin_energy = 10 .^ log_bin_centers .* 1e18
    bin_width = diff(bins)

    # Auger spectrum parameterization
    J_auger = auger_spectrum_2020.(bin_energy .* eV)  # [1 / (km^2 sr yr eV)]

    # Trapezoidal integration over energy to get counts
    # Note: the Auger data is Jraw, but the spectral fit computes J
    J_bins = auger_spectrum_2020.(bins .* eV)  # [1 / (km^2 sr yr eV)]
    count_rate = bin_width * eV .* (J_bins[1:end-1] .+ J_bins[2:end]) ./ 2 # [1 / (km^2 sr yr)]
    count = floor.(Int, count_rate * 60400 * km^2 * sr * yr)  # effective experiment sampling area and time

    # Plot
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["axes.axisbelow"] = true
    rcParams["axes.grid"] = true
    rcParams["grid.linestyle"] = "dashed"
    rcParams["xtick.direction"] = "in"
    rcParams["ytick.direction"] = "in"

    fig, axs = plt.subplots(2, 1, figsize=(6, 8))
    axs[1].loglog(bin_energy, ustrip.(J_auger), label="Auger 2020 parameterization")
    axs[1].set_ylabel("J (E) [1 / (km^2 yr sr eV)]")
    axs[1].set_ylim(1e-25, 1e-17)
    axs[1].legend()
    axs[2].scatter(bin_energy, raw_count, color="r",  label="Raw counts, Fig. 7, Auger 2020")
    axs[2].loglog(bin_energy, count, color="k", linestyle="--", label="Counts integrated from J")
    [axs[2].annotate(rc, (E, rc), rotation=30, color="r") for (E, rc) in zip(bin_energy, raw_count)]
    [axs[2].annotate(c, (E, c), ha="right", va="top", rotation=30) for (E, c) in zip(bin_energy, count)]
    axs[2].set_xlabel("Energy [eV]")
    axs[2].set_ylabel("Events per bin")
    axs[2].set_xlim(1.5e18, 2e20)
    axs[2].set_ylim(5e-1, 1e6)
    axs[2].legend()
    fig.savefig("$(@__DIR__)/../figs/auger_comparison.png")
end

plot_auger_2020_comparison()