using CoRaLS
using PyPlot


function test_fresnel_reflection()

    # we test the full range of incident angles
    θ_i = 0:0.2:90

    # use representative values for regolith and ice
    nrego = regolith_index(StrangwayIndex(), 6.0m)
    nice = 1.305

    # make the plot
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.plot(θ_i, [CoRaLS.fresnel_rpar(deg2rad(θ), nrego, nice) for θ in θ_i], label=L"$r_{\parallel}$")
    ax.plot(θ_i, [CoRaLS.fresnel_rperp(deg2rad(θ), nrego, nice) for θ in θ_i], label=L"$r_{\perp}$")
    ax.axhline(0., c="k", ls="dotted")
    ax.set(xlabel="Incident Angle (deg)", ylabel="Fresnel Coefficient", title=L"Regolith $\rightarrow$ Ice Reflection")
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.legend()

end

function test_fresnel_transmission()

    # use representative values for surface regolith and vacuum
    ni = regolith_index(StrangwayIndex(), 0.0m)
    nt = 1.0

    # calculate the TIR angle
    θ_c = asin(nt / ni)

    # we test the full range of incident angles
    θ_i = deg2rad(0.01):deg2rad(0.05):θ_c

    # make the plot
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.plot(rad2deg.(θ_i), [CoRaLS.fresnel_tpar(θ, ni, nt) for θ in θ_i], label=L"$t_{\parallel}$")
    ax.plot(rad2deg.(θ_i), [CoRaLS.fresnel_tperp(θ, ni, nt) for θ in θ_i], label=L"$t_{\perp}$")
    ax.axhline(0., c="k", ls="dotted")
    ax.set(xlabel="Incident Angle (deg)", ylabel="Fresnel Coefficient", title=L"Regolith $\rightarrow$ Vacuum Transmission")
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.legend()

end

function test_divergence_coefficient()

    # use representative values for surface regolith and vacuum
    ni = regolith_index(StrangwayIndex(), 0.0m)
    nt = 1.0

    # calculate the TIR angle
    θ_c = asin(nt / ni)

    # we test the full range of incident angles
    θ_i = deg2rad(0.01):deg2rad(0.05):θ_c

    # make the plot
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.plot(rad2deg.(θ_i), [CoRaLS.divergence_tpar(FarFieldDivergence(), θ_i[i], ni, 6.0m, 25.0km) for i in 1:length(θ_i)],
            label=L"$T_{\parallel}$")
    ax.plot(rad2deg.(θ_i), [CoRaLS.divergence_tperp(FarFieldDivergence(), θ_i[i], ni, 6.0m, 25.0km) for i in 1:length(θ_i)],
            label=L"$T_{\perp}$")
    ax.axvline(rad2deg(θ_c), c="k", ls="dotted")
    ax.set(xlabel="Incident Angle (deg)", ylabel="Fresnel Coefficient w/ Far-Field Divergence", ylim=[0., 1.0],
           title=L"Regolith $\rightarrow$ Vacuum Transmission")
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.legend(loc="upper right")

end
