using CoRaLS
using PyPlot
using Unitful: eV, km, sr, yr, ustrip
using Test

"""
    plot_failure_modes(AΩ::Acceptance)

Plot the failure modes of the CoRaLS acceptance calculation.

Arguments:
- `AΩ::Acceptance`: An object containing the acceptance calculation results.

This function generates a stacked bar chart illustrating the different failure modes of the acceptance calculation.
"""
function plot_failure_modes(AΩ::Acceptance)

    # create the figure
    fig, axs = plt.subplots(2, figsize=(10, 7))

    # the labels and colors
    labels = ["TIR", "XmaxAfterIce", "NoXmax", "Upgoing", "NotVisible", "OffTarget"]
    cm = get_cmap(:tab20)
    colors = [cm(i) for i in 1:7]


    # Count trials that failed
    ntrials = AΩ.ntrials * (length(AΩ.energies) - 1)
    sizes_direct = [0 for _ in 1:7]
    sizes_reflected = [0 for _ in 1:7]
    for signal in AΩ.failed_direct
        sizes_direct[Int(signal)] += 1
    end
    for signal in AΩ.failed_reflected
        sizes_reflected[Int(signal)] += 1
    end


    # Horizontal bar charts of failure modes in percent of total runs
    for i in 1:7
        dpct = round(sizes_direct[i] / ntrials * 100, digits=2)
        axs[1].barh(labels[i], sizes_direct[i] / ntrials * 100, color=colors[i], label="$(dpct)% $(labels[i])")
        rpct = round(sizes_reflected[i] / ntrials * 100, digits=2)
        axs[2].barh(labels[i], sizes_reflected[i] / ntrials * 100, color=colors[i], label="$(rpct)% $(labels[i])%")
    end
    axs[1].set_xlabel("Direct failure Rate [%]")
    axs[1].set_xlim([0, 10])
    axs[1].legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0)
    axs[2].set_xlabel("Reflected failure Rate [%]")
    axs[2].set_xlim([0, 10])
    axs[2].legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0)
    fig.tight_layout()
    return fig, axs

end