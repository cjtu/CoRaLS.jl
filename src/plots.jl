using PyPlot

const colordirect = "#F88D07"
const colorrefl = "#0772F8"

"""
`plots.jl` - Visualization Tools for CoRaLS Simulation Data

Key Functionalities:
- Plotting incident angles at the ice and the surface, along with the corresponding Fresnel coefficients.
- Visualizing polarization angles for direct and reflected events, using both standard and polar coordinates.
- Generating histograms for off-axis angles of both direct and reflected events.
- Displaying direct and reflected acceptance of the CoRaLS project through line plots.
- Creating differential event spectra plots, which are crucial for analyzing the frequency and distribution of cosmic ray events over different energy levels.
"""


"""
    plot_incident_angles(AΩ::Acceptance)

Plot the incident angles at the ice and the surface for both direct and reflected events, along with the corresponding Fresnel coefficients.

Arguments:
- `AΩ::Acceptance`: An object containing the acceptance calculation results.

Generates a two-part plot. The left plot shows histograms of incident angles at the surface for direct and reflected events. The right plot shows incident angles at the ice layer for reflected events. Behind these histograms, the function plots the transmission and reflection coefficients as a function of the incident angle.
"""
function plot_incident_angles(AΩ)

    # create the figure
    fig, axes = plt.subplots(1, 2, figsize=(7, 3))

    # the bins for the tx histogram
    bins = 0:1:90

    # create the right-hand axes for the Fresnel coeffs
    raxes = [axes[1].twinx(), axes[2].twinx()]

    # histogram the observed incident angles
    axes[1].hist(AΩ.direct.θ_i, histtype="step",
        label=L"Direct $\theta_i$",
        color=colordirect, density=true, bins=bins, zorder=3)
    axes[1].hist(AΩ.reflected.θ_i, histtype="step",
        label=L"Refl. $\theta_i$", ls="dashed",
        color=colorrefl, density=true, bins=bins, zorder=3)

    # calculate the transmission coefficients at the surface
    Nsurf = regolith_index(StrangwayIndex(), 0.0m)
    Θ = 0:0.1:90.0 # incident angles at the surface
    T = [surface_transmission(GaussianRoughness(2.0),
        MixedFieldDivergence(),
        θ, Nsurf, 6.0m, 100.0km) for θ in deg2rad.(Θ)]

    # plot the transmission coefficients behind
    raxes[1].plot(Θ, [t[1] for t in T], alpha=0.6, zorder=0,
        color="grey", ls="dashed", label=L"Direct $t_{\parallel}$")
    raxes[1].plot(Θ, [t[2] for t in T], alpha=0.6, zorder=0,
        color="darkgrey", ls="dotted", label=L"Direct $t_{\perp}$")

    # histogram the observed incident angles at the ice
    axes[2].hist(AΩ.reflected.θ_ice, histtype="step",
        label=L"$\theta_{\mathrm{ice}}$",
        color=colorrefl, density=true, bins=bins, zorder=3)

    # calculate the reflection coefficients behind
    Nreg = regolith_index(StrangwayIndex(), 6.0m)
    Nice = 1.305
    Θ = 0.0:0.01:rad2deg(asin(Nice / Nreg)) # incident angles at the ice
    rpar = [fresnel_coeffs[1] for θ in deg2rad.(Θ)]
    rperp = [fresnel_coeffs[2] for θ in deg2rad.(Θ)]

    # plot the reflection coefficients behind
    raxes[2].plot(Θ, abs.(rpar), alpha=0.6, zorder=0,
        color="grey", ls="dashed", label=L"$|r_{\parallel}|$")
    raxes[2].plot(Θ, abs.(rperp), alpha=0.6, zorder=0,
        color="darkgrey", ls="dotted", label=L"$|r_{\perp}|$")

    # set the axis ticks and ranges
    axes[1].set(xlim=[0, 60.0], ylim=[1e-5, 0.9])
    axes[1].set_xticks([0, 15.0, 30.0, 45.0, 60.0])
    axes[2].set(xlim=[0, 60.0], ylim=[1e-5, 0.9])
    axes[2].set_xticks([0, 15.0, 30.0, 45.0, 60.0])
    [rax.set_ylim([0.0, 1.1]) for rax in raxes]
    [ax.set_yscale("log") for ax in axes]

    # labels
    axes[1].set_xlabel("Incident Angle at Surface [deg]")
    axes[1].set_ylabel("Event Density")
    raxes[1].set_ylabel("Fresnel Coeff.")
    axes[2].set_xlabel("Incident Angle at Ice [deg]")
    axes[2].set_ylabel("Event Density")
    raxes[2].set_ylabel("Fresnel Coeff.")

    # leends
    axes[1].legend(loc="center left")
    raxes[2].legend(loc="upper left")


    # grids
    [ax.set_axisbelow(true) for ax in axes]
    [ax.grid(which="both", ls="dashed") for ax in axes]
    fig.tight_layout()

end

"""
    plot_polarization_angle(AΩ::Acceptance)

Plot histograms of polarization angles for direct and reflected events.

Arguments:
- `AΩ::Acceptance`: An object containing the acceptance calculation results.
"""
function plot_polarization_angle(AΩ)

    # create the figure - plot this in polar coordinates
    fig = plt.figure(figsize=(6, 3))
    axes = [fig.add_subplot(1, 2, 1, polar=true),
        fig.add_subplot(1, 2, 2, polar=true)]

    # an array of colors
    colors = [colordirect, colorrefl]

    # axes[1].hist(AΩ.direct.θpol .|> deg2rad, density=true, color=colordirect, histtype="step")
    # axes[2].hist(AΩ.reflected.θpol .|> deg2rad, density=true, color=colorrefl, histtype="step")

    axes[1].hist(AΩ.reflected.θpol .|> deg2rad, density=true, color=colorrefl, histtype="step")
    axes[2].hist(AΩ.reflected.θpolsub .|> deg2rad, density=true, color="navy", histtype="step")

    # for (i, events) = enumerate([AΩ.direct, AΩ.reflected])

    #     # compute the polarization angle - 0 == H-Pol, 90 = V-Pol
    #     χ = deg2rad.(polarization_angle(events))

    #     # histogram the polarization angle
    #     # axes[i].hist(χ[trig.(events) .== true], bins=45, density=true, color=colors[i], histtype="step")
    #     axes[i].hist(χ, density=true, color=colors[i], histtype="step")

    # end
    # end loop over events

    # and make it pretty
    [ax.set(xlabel="Polarization Angle [deg]") for ax in axes]
    [ax.set_xticklabels(["+H", "", "+V", "", "-H", "", "-V", ""]) for ax in axes]
    axes[1].set_title("Top of Ice")
    axes[2].set_title("Bottom of Ice")
    # fig.suptitle("UH MC")
    fig.tight_layout()

    return fig, axes

end

"""
    plot_offaxis_angle(AΩ::Acceptance)

Generate a plot showing the off-axis angle distribution of direct and reflected events in the CoRaLS simulations.

Arguments:
- `AΩ::Acceptance`: An object containing the acceptance calculation results.
"""
function plot_offaxis_angle(AΩ)

    # create the figure
    fig, ax = plt.subplots(figsize=(4, 4))

    # plot the histograms
    ax.hist(AΩ.direct.ψ, histtype="step", color=colordirect, alpha=0.3)
    ax.hist(AΩ.reflected.ψ, histtype="step", color=colorrefl, alpha=0.3)
    ax.hist(AΩ.direct.ψ[AΩ.direct.triggered.==true],
        histtype="step", color=colordirect, label="Direct")
    ax.hist(AΩ.reflected.ψ[AΩ.reflected.triggered.==true],
        histtype="step", color=colorrefl, label="Reflected")

    # and make it pretty
    ax.set(xlabel="Off-Axis Angle [deg]", ylabel="Counts",
        xticks=[0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0],
        yscale="log", xlim=[0.0, 180.0])

    ax.legend(loc="upper right")
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")

    return fig, ax
end

"""
    plot_acceptance(AΩ::Acceptance)

Plot the direct and reflected acceptance of the CoRaLS project.

Arguments:
- `AΩ::Acceptance`: An object containing the acceptance calculation results.

This function generates a line plot illustrating the acceptance rates for both direct and reflected events across different energy levels.
"""
function plot_acceptance(AΩ::Acceptance)

    # create the figure
    fig, ax = plt.subplots(figsize=(4, 4))


    # plot the lines
    ax.plot(18.0 .+ log10.(0.5 * (AΩ.energies[1:end-1] + AΩ.energies[2:end]) ./ 1.0EeV),
        AΩ.dAΩ ./ km^2 ./ sr,
        color=colordirect, label="Direct")
    ax.plot(18.0 .+ log10.(0.5 * (AΩ.energies[1:end-1] + AΩ.energies[2:end]) ./ 1.0EeV),
        AΩ.rAΩ ./ km^2 ./ sr,
        color=colorrefl, label="Reflected")

    # and finally pretty it up
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.set(xlabel=L"Energy [$\log_{10}$(eV)]",
        ylabel=L"Acceptance [km$^2$ sr]",
        yscale="log")
    ax.legend(loc="upper left")

    return fig, ax

end

"""
    plot_differential_spectrum(AΩ::Acceptance, T::Unitful.Time)

Create a plot of the differential event spectrum of CoRaLS events given an acceptance object and a mission duration.

Arguments:
- `AΩ::Acceptance`: An object containing the acceptance calculation results.
- `T::Unitful.Time`: The total observation time of the experiment.

This function plots the estimated number of events per energy bin over the duration of the mission.
"""
function plot_differential_spectrum(AΩ::Acceptance, T)

    # create the figure
    fig, ax = plt.subplots(figsize=(4, 4))

    # the fraction of events that have PSRs
    psr_frac = psr_fraction(AΩ.altitude)

    # calculate the differential event spectrums
    # devents = differential_spectrum(AΩ.energies, psr_frac*AΩ.dAΩ, 2yr)
    revents = differential_spectrum(AΩ.energies, (30_000 / 20_000) * psr_frac * AΩ.rAΩ, 2yr)
    PIRevents = differential_spectrum(AΩ.energies, 6.3 * psr_frac * AΩ.rAΩ, 2yr)

    # plot the direct events
    # ax.hist(18.0 .+ log10.(AΩ.energies[1:end-1] ./ 1.0EeV),
    #          18.0 .+ log10.(AΩ.energies ./ 1.0EeV),
    #          weights=devents,
    #          histtype="step", color="#F88D07",
    #         label="Direct, N=$(Int(round(sum(devents))))")
    ax.hist(18.0 .+ log10.(AΩ.energies[1:end-1] ./ 1.0EeV),
        18.0 .+ log10.(AΩ.energies ./ 1.0EeV),
        weights=PIRevents,
        histtype="step", color="#F88D07",
        label="PIRs, N=$(Int(round(sum(PIRevents))))")

    # plot the reflected events
    ax.hist(18.0 .+ log10.(AΩ.energies[1:end-1] ./ 1.0EeV),
        18.0 .+ log10.(AΩ.energies ./ 1.0EeV),
        weights=revents,
        histtype="step", color="#0772F8",
        label="PSRs, N=$(Int(round(sum(revents))))")

    # and finally pretty it up
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.set(xlabel=L"Energy [$\log_{10}$(eV)]",
        xlim=[18.0 + log10(AΩ.energies[1] / EeV), 18.0 + log10(AΩ.energies[end] / EeV)],
        ylabel="Events per bin / $(T)",
        xticks=[17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0],
        yscale="log")
    ax.legend(loc="upper right")

    return fig, ax

end
