using CoRaLS
using Test
using NumericalIntegration
using PyPlot
using LaTeXStrings
using Unitful: MHz, m, km, EeV, μV, NoUnits

function plot_field_models()

    # create a 1x2 set of axes
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))

    # we sample the full range of angles
    Ψ = 0.0:0.1:180.0

    # the maximum frequencies that we calculate for
    ν_max = [50., 80., 150., 300., 500., 1000]MHz

    # get the colors that we sequentially use
    colors = plt.cm.get_cmap("magma")(range(0, stop=1, length=length(ν_max)+1))

    # the left plot is 10EeV and the right plot is 100EeV
    energies = [10.0, 100.0]EeV

    # the distance in regolith and vacuum
    Drego = 6.0m
    Dvacuum = 25.0km

    # produce the two plots in a loop
    for iax = 1:2

        ax = axes[iax] # the current axis that we plot

        # we now produce a curve for each maximum frequency
        for inu = 1:length(ν_max)

            # plot the FORTE model
            ax.plot(Ψ, [integrate(regolith_field(FORTE{GaussianProfile}, energies[iax], deg2rad(ψ), Drego, Dvacuum;
                                                 ν_min=0.0MHz, ν_max=ν_max[inu])...) for ψ in Ψ] ./
                                                     (μV/m) .|> NoUnits,
                    c=colors[inu, :], ls="solid", label="$(Int(ν_max[inu] / MHz)) MHz")

            # plot the JAM model
            ax.plot(Ψ, [integrate(regolith_field(JAM(), energies[iax], deg2rad(ψ), Drego, Dvacuum;
                                                 ν_min=0.0MHz, ν_max=ν_max[inu])...) for ψ in Ψ] ./
                                                     (μV/m) .|> NoUnits,
                    c=colors[inu, :], ls="dashed")

        end # end loop over frequencies


    end

    # set some axes limits
    axes[1].set(ylim=[10., 14_000], xlim=[0, 180.], yscale="log")
    axes[2].set(ylim=[100., 140_000], xlim=[0., 180.], yscale="log")

    # set some titles
    axes[1].set_title("$(energies[1] / EeV) EeV", fontsize=10)
    axes[2].set_title("$(energies[2] / EeV) EeV", fontsize=10)

    # set the labels
    [ax.set(xlabel="Off-Axis Angle (deg)", ylabel=L"Electric Field Strength ($\mu$V/m)") for ax in axes]

    # turn on the grid
    [ax.grid(which="both", ls="dashed") for ax in axes]

    # turn on the legend
    [ax.legend() for ax in axes]

    # and a title
    fig.suptitle("FORTE (solid) vs. JAM (dashed) [$(Int(Drego / m))m + $(Int(Dvacuum / km))km]")

    # and save
    fig.savefig("figs/efield_comparison.pdf")
    fig.savefig("figs/efield_comparison.png")

    return true

end

@testset "efield.jl" begin
    plot_field_models()
end
