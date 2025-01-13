"""
This file performs various comparisons between different configurations
of the CoRaLS MC.
"""

using CoRaLS
using ProgressMeter
using PyPlot
using LaTeXStrings
using DelimitedFiles
using Unitful: km, sr, EeV, MHz, m, μV, cm

# the number of trials that we use for each run
ntrials = 100_000
#ntrials = 1_000

# and the number of bins in energy for each run
nbins = 20

"""
Baseline acceptance curve.

This isn't the *correct* curve. It is the original curve before starting
the more sophisticated MC studies in August 2021.
"""
function baseline_test()

    # calculate with a constant refractive index
    AΩ = acceptance(ntrials, nbins;
                    altitude=20.0km, ice_depth=6.0m)

    # and now create the plot
    fig, ax = plt.subplots(figsize=(4, 4))

    # plot the direct simulations
    ax.plot(18.0 .+ log10.(0.5*(AΩ.energies[1:end-1] + AΩ.energies[2:end]) / 1.0EeV), AΩ.dAΩ / 1km^2 / 1sr,
            c="lightskyblue", ls="solid", label="Direct")

    # plot the reflected simulations
    ax.plot(18.0 .+ log10.(0.5*(AΩ.energies[1:end-1] + AΩ.energies[2:end]) / 1.0EeV), AΩ.rAΩ / 1km^2 / 1sr,
            c="crimson", ls="solid", label="Reflected")

    # and some some info
    ax.set(xlabel=L"Energy [$\log_{10}$(eV)]", ylabel=L"Acceptance [km$^2$ sr]",
           xlim=[17.5, 20.5],
           ylim=[1e0, 2e2],
           xticks=[17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5],
           yscale="log", title="UH Baseline Acceptance: 30km")

    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.legend()
    fig.savefig("test/figs/UH_baseline_acceptance.png")

end

"""
Current best-guess acceptance curve.


"""
function acceptance_test()

    # the altitudes that we plot
    altitudes = [10.0km, 20.0km, 30.0km]

    # colormaps for the direct and reflected events
    direct_cmap = plt.cm.get_cmap("Blues")(range(1, 0, length=length(altitudes)+2))
    refl_cmap = plt.cm.get_cmap("Reds")(range(1, 0, length=length(altitudes)+2))

    # and now create the plot
    fig, ax = plt.subplots(1, 2, sharey=true, figsize=(6, 3))

    # create another plot for the differential number number of events
        
    for i = 1:length(altitudes)
        # calculate with a constant refractive index
        A = acceptance(ntrials, nbins;
                                              altitude=altitudes[i], ice_depth=6.0m,
                                              min_energy=1.0EeV,
                                              ν_max=1000MHz,
                                              dν = 30MHz,
                                              trigger=gaussian_trigger(sqrt(1000 / 300.)*0.5*67μV/m),
                                              save_events=false)

        # plot the direct simulations
        ax[1].plot(18.0 .+ log10.(0.5*(A.energies[1:end-1] + A.energies[2:end]) / 1.0EeV), A.dAΩ / 1km^2 / 1sr,
                   c=direct_cmap[i+1, :], ls="solid", label="$(altitudes[i])")

        # plot the reflected simulations
        ax[2].plot(18.0 .+ log10.(0.5*(A.energies[1:end-1] + A.energies[2:end]) / 1.0EeV), A.rAΩ / 1km^2 / 1sr,
                   c=refl_cmap[i+1, :], ls="solid", label="$(altitudes[i])")


        # calculate the fraction of area that is covered by PSRs
        Apsr = psr_fraction(altitudes[i])

        # get the differential number of direct and reflected events
        devents = sum( differential_spectrum(A.energies, A.dAΩ, 2yr) )
        revents = sum( differential_spectrum(A.energies, Apsr*A.rAΩ, 2yr) )


    end

    # and some some info
    ax[1].set(xlabel=L"Energy [$\log_{10}$(eV)]", ylabel=L"Acceptance [km$^2$ sr]",
              xlim=[18.5, 20.5], yscale="log", title="Direct")
    ax[2].set(xlabel=L"Energy [$\log_{10}$(eV)]", #ylabel=L"Acceptance [km$^2$ sr]",
              xlim=[18.5, 20.5], yscale="log", title="Reflected")

    for a in ax
        a.set_axisbelow(true)
        a.grid(which="both", linestyle="dashed")
        a.legend()
    end

    fig.savefig("test/figs/UH_acceptance.png")

end

"""
Compare the CoRaLS.jl against the JPL MC.
"""
function compare_against_jpl_test()

    # load the JPL acceptances at each energy
    JPLd30 = readdlm("docs/data/JPL_direct_30km.csv", ',', Float64)
    JPLr30 = readdlm("docs/data/JPL_reflected_30km.csv", ',', Float64)

    # calculate the best fit to the JPL simulations
    A = acceptance(ntrials, nbins;
                                          altitude=30.0km, ice_depth=10.0m,
                                          geometrymodel=ScalarGeometry(),
                                          indexmodel=SurfaceDeepIndex(),
                                          fieldmodel=ARW(),
                                          densitymodel=ConstantDensity(),
                                          slope=GaussianSlope(7.6),
                                          trigger=gaussian_trigger())
    energies = A.energies
    dAΩ = A.dAΩ
    rAΩ = A.rAΩ

    # and now create the plot
    fig, ax = plt.subplots(figsize=(4, 4))

    # plot the direct simulations
    ax.plot(JPLd30[:, 1], JPLd30[:, 2], c="blue", ls="dashed", label="JPL: Direct")
    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), dAΩ / 1km^2 / 1sr,
            c="orange", ls="dashed", label="UH: Direct")

    # plot the reflected simulations
    ax.plot(JPLr30[:, 1], JPLr30[:, 2], c="blue", ls="solid", label="JPL: Refl.")
    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), rAΩ / 1km^2 / 1sr,
            c="orange", ls="solid", label="UH: Refl.")

    # and some some info
    ax.set(xlabel=L"Energy [$\log_{10}$(eV)]", ylabel=L"Acceptance [km$^2$ sr]",
           xlim=[18.5, 20.5], yscale="log", ylim=[1e2, 1e6], title="JPL/UH Comparison")

    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.legend()
    fig.savefig("test/figs/JPL_UH_Acceptance.png")

    return energies, dAΩ, rAΩ

end

"""
Compare different refractive models.

We compare a constant refractivity and a "surface+deep" refractivity.
"""
function compare_refractive_index_test()

    # we produce runs for constant refractive index
    # and for a two-part refractive index

    # calculate with a constant refractive index
    Aconst = acceptance(ntrials, nbins;
                                                    indexmodel=ConstantIndex())
    dAΩconst = Aconst.dAΩ
    rAΩconst = Aconst.rAΩ

    # calculate with a deep+surface refractive index
    Adeep = acceptance(ntrials, nbins;
                                                  indexmodel=SurfaceDeepIndex())
    dAΩdeep = Adeep.dAΩ
    rAΩdeep = Adeep.rAΩ

    # calculate with a deep+surface refractive index
    Aexp = acceptance(ntrials, nbins;
                                                  indexmodel=StrangwayIndex())
    dAΩexp = Aexp.dAΩ
    rAΩexp = Aexp.rAΩ
    energies = Aexp.energies

    # create our plot
    fig, ax = plt.subplots(figsize=(8, 4))

    # plot the direct
    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), dAΩconst / 1km^2 / 1sr,
            c="lightskyblue", ls="dashed", label="Direct: Constant")
    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), dAΩdeep / 1km^2 / 1sr,
            c="lightskyblue", ls="dotted", label="Direct: Surface+Deep")
    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), dAΩexp / 1km^2 / 1sr,
            c="lightskyblue", ls="solid", label="Direct: Exponential")

    # plot the reflected
    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), rAΩconst / 1km^2 / 1sr,
            c="crimson", ls="dashed", label="Reflected: Constant")
    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), rAΩdeep / 1km^2 / 1sr,
            c="crimson", ls="dotted", label="Reflected: Surface+Deep")
    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), rAΩexp / 1km^2 / 1sr,
            c="crimson", ls="solid", label="Reflected: Exponential")

    # add some labels
    ax.set(xlabel=L"Energy [$\log_{10}$(eV)]", ylabel=L"Acceptance [km$^2$ sr]",
           title="Regolith Refractive Index Models", yscale="log")

    # add our grid and legend
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.legend()
    fig.savefig("test/figs/refractive_index_acceptance.png")


end

"""
Compare different surface divergence models.

"""
function compare_surface_divergence_test()

    # we produce runs for constant refractive index
    # and for a two-part refractive index

    # create our plot
    fig, ax = plt.subplots(figsize=(8, 4))

    # add some labels
    ax.set(xlabel=L"Energy [$\log_{10}$(eV)]", ylabel=L"Acceptance [km$^2$ sr]",
           title="Surface Divergence Models", yscale="log")

    # calculate with far-field
    Afar = acceptance(ntrials, nbins;
                                                divergencemodel=FarFieldDivergence())
    dAΩfar = Afar.dAΩ
    rAΩfar = Afar.rAΩ
    energies = Afar.energies

    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), dAΩfar / 1km^2 / 1sr,
            c="lightskyblue", ls="dashed", label="Direct: FORTE")
    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), rAΩfar / 1km^2 / 1sr,
            c="crimson", ls="dashed", label="Reflected: FORTE")

    # calculate with mixed-field
    Amixed = acceptance(ntrials, nbins;
                                                    divergencemodel=MixedFieldDivergence())
    dAΩmixed = Amixed.dAΩ
    rAΩmixed = Amixed.rAΩ
    energies = Amixed.energies

    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), dAΩmixed / 1km^2 / 1sr,
            c="lightskyblue", ls="solid", label=L"Direct: $\sqrt{n}$ FORTE")
    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), rAΩmixed / 1km^2 / 1sr,
            c="crimson", ls="solid", label=L"Reflected: $\sqrt{n}$ FORTE")

    # calculate with near-field
    Anear = acceptance(ntrials, nbins;
                                                  divergencemodel=NearFieldDivergence())
    dAΩnear = Anear.dAΩ
    rAΩnear = Anear.rAΩ
    energies = Anear.energies

    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), dAΩnear / 1km^2 / 1sr,
            c="lightskyblue", ls="dotted", label=L"Direct: $n$ FORTE")
    ax.plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), rAΩnear / 1km^2 / 1sr,
            c="crimson", ls="dotted", label=L"Reflected: $n$ FORTE")


    # add our grid and legend
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.legend()
    fig.savefig("test/figs/surface_divergence_tests.png")

end

"""
Compare different bandwidths.

"""
# function compare_bandwidth_test()

#     # we produce runs for constant refractive index
#     # and for a two-part refractive index

#     # create our plot
#     fig, axes = plt.subplots(1, 2, figsize=(8, 4))

#     # add some labels
#     [ax.set(xlabel=L"Energy [$\log_{10}$(eV)]",
#             ylabel=L"Acceptance [km$^2$ sr]",
#             xlim=[17.5, 20.5],
#             ylim=[1e1, 1e5],
#             yscale="log") for ax in axes]

#     # set some titles
#     axes[1].set_title("Direct")
#     axes[2].set_title("Reflected")

#     # use a colormap for the colors
#     colors = plt.cm.get_cmap("Blues_r")(range(0, 1, length=5))

#     # calculate with 300 MHz bandwidth
#     energies, dAΩ300, _, rAΩ300, _ = acceptance(ntrials, nbins;
#                                                 min_energy=0.3EeV,
#                                                 ν_max=300MHz,
#                                                 roughnessmodel=GaussianRoughness(2.0cm),
#                                                 trigger=gaussian_trigger(sqrt(300 / 300.) * 0.5*67μV/m))

#     axes[1].plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), dAΩ300 / 1km^2 / 1sr,
#             c=colors[2, :], ls="solid", label="300 MHz")
#     axes[2].plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), rAΩ300 / 1km^2 / 1sr,
#             c=colors[2, :], ls="solid", label="300 MHz")

#     # calculate with 500 MHz bandwidth
#     energies, dAΩ500, _, rAΩ500, _ = acceptance(ntrials, nbins;
#                                                 min_energy=0.3EeV,
#                                                 ν_max=500MHz,
#                                                 dν=15MHz,
#                                                 roughnessmodel=GaussianRoughness(2.0cm),
#                                                 trigger=gaussian_trigger(sqrt(500 / 300.) * 0.5*67μV/m))

#     axes[1].plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), dAΩ500 / 1km^2 / 1sr,
#             c=colors[3, :], ls="solid", label="500 MHz")
#     axes[2].plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), rAΩ500 / 1km^2 / 1sr,
#                c=colors[3, :], ls="solid", label="500 MHz")

#     # calculate with 1000 MHz bandwidth
#     energies, dAΩ1000, _, rAΩ1000, _ = acceptance(ntrials, nbins;
#                                                   min_energy=0.3EeV,
#                                                   ν_max=1000MHz,
#                                                   dν=30MHz,
#                                                   roughnessmodel=GaussianRoughness(2.0cm),
#                                                   trigger=gaussian_trigger(sqrt(1000 / 300.) * 0.5*67μV/m))

#     axes[1].plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), dAΩ1000 / 1km^2 / 1sr,
#             c=colors[4, :], ls="solid", label="1000 MHz")
#     axes[2].plot(18.0 .+ log10.(0.5*(energies[1:end-1] + energies[2:end]) / 1.0EeV), rAΩ1000 / 1km^2 / 1sr,
#             c=colors[4, :], ls="solid", label="1000 MHz")



#     # add our grid and legend
#     for ax in axes
#         ax.set_axisbelow(true)
#         ax.grid(which="both", linestyle="dashed")
#         ax.legend()
#     end
#     fig.savefig("test/figs/bandwidth_comparison.png")

#     return energies, dAΩ1000, rAΩ1000

# end

"""
Compare the event rate as a function of surface roughness.
"""
function compare_surface_roughness()

    # we use a different set of trials and bins here to make it faster
    #ntrials = 100_000
    ntrials = 1_000
    nbins = 15

    # the roughness values that we calculate for
    Σ = 0.0:0.5:7.0

    # the array where we store the number of events
    events = zeros(length(Σ))

    # calculate the fraction of area that is covered by PSRs
    Apsr = psr_fraction(20.0km)

    # loop over each roughness
    for i = 1:length(Σ)

        # calculate the acceptance at these energies
        A = acceptance(ntrials, nbins;
                                              ν_max=1000MHz,
                                              min_energy=1.0EeV,
                                              altitude=20.0km,
                                              trigger=gaussian_trigger(sqrt(1000 / 300.) * 0.5*67μV/m),
                                              roughnessmodel=GaussianRoughness(Σ[i]))
        dAΩ = A.dAΩ
        rAΩ = A.rAΩ
        energies = A.energies

        # and calculate the total number of events
        events[i] = sum( differential_spectrum(energies, Apsr*rAΩ, 2yr) )

    end

    # create our plot
    fig, ax = plt.subplots(figsize=(4, 4))

    # use a colormap for the colors
    colors = plt.cm.get_cmap("Blues_r")(range(0, 1, length=5))

    # plot the curve
    ax.plot(Σ , events, color=colors[2, :], ls="solid")
    ax.scatter(Σ , events, color=colors[2, :], s=1.8)

    # add some labels
    ax.set(xlabel=L"Surface RMS, $\sigma$",
            ylabel="Number of Events / 2yr",
            title="Surface Roughness (20km, BW: 1000 MHz)")

    # add our grid and legend
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")

    fig.savefig("test/figs/roughness_comparison.png")

end


"""
Compare the event rate as a function of trigger threshold.
"""
function compare_thresholds()

    # we use a different set of trials and bins here to make it faster
    #ntrials = 60_000
    ntrials = 1_000
    nbins = 20

    # the threshold for each simulation
    thresholds = range(10, stop=230, length=15)*(μV/m)

    # the array where we store the number of events
    events = zeros(length(thresholds))

    # calculate the fraction of area that is covered by PSRs
    Apsr = psr_fraction(20.0km)

    # loop over each roughness
    for i = 1:length(thresholds)

        # calculate the acceptance at these energies
        A = acceptance(ntrials, nbins;
                                              min_energy=0.1EeV,
                                              altitude=20.0km,
                                              ν_max=1000MHz,
                                              dν = 30MHz,
                                              trigger=gaussian_trigger(thresholds[i]),
                                              save_events=false)
        dAΩ = A.dAΩ
        rAΩ = A.rAΩ
        energies = A.energies

        # and calculate the total number of events
        events[i] = sum( differential_spectrum(energies, Apsr*rAΩ, 2yr) )

    end

    # create our plot
    fig, ax = plt.subplots(figsize=(4, 4))

    # use a colormap for the colors
    colors = plt.cm.get_cmap("Blues_r")(range(0, 1, length=5))

    # plot the curve
    ax.plot(thresholds ./ (μV/m), events, color=colors[2, :], ls="solid")
    ax.scatter(thresholds ./ (μV/m), events, color=colors[2, :], s=2.0)

    # add some labels
    ax.set(xlabel=L"Gaussian Beam Threshold [$\mu$V/m]",
           ylabel="Number of Events / 2yr")

    # and set a title
    plt.title("1 GHz BW @ 20km w/ Ice Roughness")

    # add our grid and legend
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")

    fig.savefig("test/figs/threshold_comparison.png")

end


"""
Compare the event rate as a function of ice refractive index
"""
function compare_ice_depth()

    # we use a different set of trials and bins here to make it faster
    #ntrials = 180_000 # 300_000
    ntrials = 1_000 # 300_000
    nbins = 10 # 16

    # scan done on Oct. 8th
    # N=6, baseline=10m
    # N=6, threshold=110m
    # N=8, baseline=18m-20m
    # N=8, threshold=120m
    # N=12, baseline=30m
    # N=12, threshold=140m
    # N=16, baseline=40m
    # N=16, threshold=160m

    # the ice-depth for each simulation
    ice_depths = [5.0, 7.5, 10.0, 12.5, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]m
    # ice_depths = [5.0, 10.0, 15.0, 80.0, 100.0, 110.0]m
    # ice_depths = [20.0, 30.0, 40.0, 50.0, 120.0, 130.0, 140.0, 150.0, 160.0]m # N=12
    # ice_depths = [30.0, 40.0, 50.0, 140.0, 150.0, 160.0, 170.0]m # N=16
    # ice_depths = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]m

    # the number of total antennas we use in each case
    Nants = [8, 12, 16]
    Ntrigs = [3, 5, 6]

    # use a colormap for the colors
    colors = plt.cm.get_cmap("Blues")(range(0, 1, length=length(Nants)+1))

    # create our plot
    fig, ax = plt.subplots(figsize=(4, 4))

    # loop over each antenna configurations
    for i = 1:length(Nants)

        # the array where we store the number of events
        events = zeros(length(ice_depths))

        # loop over each roughness
        for j = 1:length(ice_depths)

            # calculate the acceptance at these energies
            AΩ = acceptance(ntrials, nbins;
                            altitude=20.0km,
                            ice_depth=ice_depths[j],
                            min_energy=0.5EeV, max_energy=500EeV,
                            trigger=LPDA(Nant=Nants[i], Ntrig=Ntrigs[i],
                                         ν_max=800MHz,
                                         θ0=-48.0,
                                         altitude=20.0km),
                            ν_min=150MHz, ν_max=800MHz,
                            save_events=false)

            # and calculate the total number of events
            events[j] = sum( differential_spectrum(AΩ.energies, AΩ.rAΩ, 2yr) )

        end # end loop over ice depths

        scale = 1.0
        if Nants[i] == 8
            scale = 1.15
        end

        # plot the curve
        ax.plot(ice_depths ./ m, scale .* events, color=colors[i+1, :],
                ls="solid", label=label="N=$(Nants[i])")
        ax.scatter(ice_depths ./ m, scale .* events, color=colors[i+1, :], s=2.0)

    end # end loop over number of antennas

    # add some labels
    ax.set(xlabel="Ice Depth [m]",
           ylabel="Number of Events in PSRs / 2yr")

    # force the bottom to be zero
    ax.set_ylim(bottom=0.)

    # set manual ticks
    ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70])

    # and set a title
    # plt.title("1 GHz BW @ 20km w/ Ice Roughness")

    # horizontal dashed grey line at our baseline
    ax.axhline(275, ls="dashed", color="grey", label="Requirement")

    # 5m => 1 Gyr
    # 7m => 2 Gyr
    # 10m => 3 Gyr
    # add vertical lines for the various ages
    ax.axvline(5, ls="dotted", color="orange")
    ax.axvline(7, ls="dotted", color="orange")
    ax.axvline(10, ls="dotted", color="orange")



    # add our grid and legend
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")

    # and enable the legend
    ax.legend(loc="upper right")

    fig.savefig("test/figs/ice_depth_comparison.png")

end

"""
Compare the event rate as a function of ice refractive index
"""
function compare_ice_purity()

    # we use a different set of trials and bins here to make it faster
    #ntrials = 300_000
    ntrials = 1_000
    nbins = 16

    # the ice-depth for each simulation
    # ice_depths = [5.0, 7.5, 10.0, 12.5, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0]m
    purity = range(0, 1, length=10)

    # the array where we store the number of events
    events = zeros(length(purity))

    # the constants that we use to calculate the mixed ice index
    Nrego = regolith_index(StrangwayIndex(), 5.0m)
    Nice = 1.305

    # loop over each roughness
    for i = 1:length(purity)

        # calculate the refractive index of the ice at this purity
        Ntrue = exp(purity[i]*log(Nrego) + (1.0 - purity[i])*log(Nice))

        # calculate the acceptance at these energies
        AΩ = acceptance(ntrials, nbins;
                        altitude=10.0km,
                        ice_depth=5.0m,
                        Nice=Ntrue,
                        min_energy=0.5EeV, max_energy=500EeV,
                        trigger=LPDA(Nant=8, Ntrig=3, ν_max=800MHz,
                                     θ0=-48.0,
                                     altitude=10.0km),
                        ν_min=150MHz, ν_max=800MHz,
                        save_events=false)

        # and calculate the total number of events
        events[i] = sum( differential_spectrum(AΩ.energies, AΩ.rAΩ, 2yr) )

    end

    # create our plot
    fig, ax = plt.subplots(figsize=(4, 4))

    # use a colormap for the colors
    colors = plt.cm.get_cmap("Blues_r")(range(0, 1, length=5))

    # plot the curve
    ax.plot(100*purity, events, color=colors[2, :], ls="solid")
    ax.scatter(100*purity, events, color=colors[2, :], s=2.0)

    # add some labels
    ax.set(xlabel="Volume Fraction of Regolith in Ice Layer [%]",
           ylabel="Number of Events / 2yr")

    # force the bottom to be zero
    ax.set_ylim(bottom=0.)

    # and set a title
    plt.title("Ice Purity")

    # add our grid and legend
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")

    fig.savefig("test/figs/ice_purity_comparison.png")

end


"""
Compare the event rate as a function of beam pointing.
"""
function compare_beam_center()

    # we use a different set of trials and bins here to make it faster
    #ntrials = 250_000
    ntrials = 1_000
    nbins = 20

    # the number of antennas on the SC
    Ns = [6, 8, 12]
    Ntrig = [2, 3, 4]

    # the threshold for each simulation
    θs = [-90., -75, -60, -45, -36., -30., -15., -0.]
    θs = range(-90, stop=0, length=15)

    # run a large simulation at a range of energies to then later apply
    # our different trigger models to
    AΩ = acceptance(ntrials, nbins;
                    min_energy=0.5EeV, max_energy=500EeV,
                    trigger=LPDA(Nant=12, Ntrig=3, ν_max=800MHz),
                    ν_min=150MHz, ν_max=800MHz,
                    save_events=true)

    # create our plot
    fig, ax = plt.subplots(figsize=(4, 4))

    # use a colormap for the colors
    colors = plt.cm.get_cmap("Blues")(range(0, 1, length=length(Ns)+2))

    # loop over each number of antennas - this is one line
    for j = 1:length(Ns)

        # the array where we store the number of events
        events = zeros(length(θs))

        # and loop over beam-tilt
        @showprogress 1 "Retriggering N=$(Ns[j])..." for i = 1:length(θs)

            # construct a new trigger at this location
            trig = LPDA(Nant=Ns[j], Ntrig=Ntrig[j], θ0=θs[i], SNR=4.0)

            # apply the trigger and recalculate the effective area
            AΩeff = retrigger(AΩ, trig).rAΩ
            # and calculate the total number of events
            events[i] = sum( differential_spectrum(AΩ.energies, AΩeff, 2yr) )

        end # end loop over beam tilt

        # plot the curves
        ax.plot(θs, events, color=colors[j+1, :], ls="solid", label="N=$(Ns[j])")
        # ax.scatter(θs, events, color=colors[2, :], s=2.5)

    end

    # put a vertical line at the horizon
    ax.axvline(horizon_angle(20.0km) |> rad2deg, c="k", ls="dotted", label="Horizon", zorder=-3)

    # add some labels
    ax.set(xlabel="Beam Center [deg]",
           ylabel="Number of Events / 2yr")
    ax.set_xticks([-90., -75, -60, -45, -30., -15., -0.])
    ax.set_xlim([-90, 0])
    ax.legend()

    # and set a title
    plt.title("Total LPDAs (150 MHz - 800 MHz)")

    # add our grid and legend
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")

    fig.savefig("test/figs/beam_center_comparison_4sigma.png")

end

"""
Compare the event rate as a function of beam pointing.
"""
function compare_bandwidth()

    # we use a different set of trials and bins here to make it faster
    #ntrials = 60_000
    ntrials = 1_000
    nbins = 15

    # the minimum and maximum frequencies
    ν_mins = [150, 200]MHz
    ν_maxs = [300, 400, 500, 600, 700, 800, 900]MHz

    # calculate the fraction of area that is covered by PSRs
    Apsr = psr_fraction(20.0km)

    # run a large simulation at a rnnnange of energies to then later apply
    # our different trigger models to
    AΩ = acceptance(ntrials, nbins;
                    min_energy=0.5EeV, max_energy=500EeV,
                    trigger=LPDA(),
                    ν_min=100MHz, ν_max=900MHz,
                    save_events=true)

    # create our plot
    fig, ax = plt.subplots(figsize=(4, 4))

    # use a colormap for the colors
    colors = plt.cm.get_cmap("Blues")(range(0, 1, length=length(ν_mins)+2))

    # loop over each maximum bandwidth
    for j = 1:length(ν_mins)

        # where we store the reflected event counts
        events = zeros(length(ν_maxs))

        @showprogress 1 "Retriggering..."  for i = 1:length(ν_maxs)
            # construct a new trigger at this location
            trig = LPDA(Nant=4, θ0=-36.0, SNR=4.0, ν_min=ν_mins[j], ν_max=ν_maxs[i])

            # apply the trigger and recalculate the effective area
            AΩeff = retrigger(AΩ, trig).rAΩ

            # and calculate the total number of events
            events[i] = 1.5 * sum( differential_spectrum(AΩ.energies, Apsr*AΩeff, 2yr) )

            println("$(ν_mins[j]):$(ν_maxs[i]) == $(events[i])")

        end

        # and plot it
        ax.plot(ν_maxs ./ MHz, events, c=colors[j+1, :], label=L"$\nu_{\mathrm{min}}=$"*" $(ν_mins[j])")

    end

    ax.set(xlabel="Maximum Frequency (MHz)", ylabel="PSR Events", title="Bandwidth Scan",
           xticks=[300, 450, 600, 750, 900])
    ax.set_axisbelow(true)
    ax.grid(which="both", ls="dashed")
    ax.legend()
    plt.show()

    fig.savefig("test/figs/bandwidth_comparison.png")

end

"""
Compare the event rate as a function of beam threshold and width
"""
function scan_width_threshold()

    # we use a different set of trials and bins here to make it faster
    #ntrials = 90_000
    ntrials = 1_000
    nbins = 12

    # the beamwidth for each simulation
    BWs = [20., 25., 30., 40., 50.]
    # and the threshold
    thresholds = [20, 30, 40, 50] .* (μV/m)

    # calculate the fraction of area that is covered by PSRs
    Apsr = psr_fraction(20.0km)

    # create our plot
    fig, ax = plt.subplots(figsize=(4, 4))

    # use a colormap for the colors
    colors = plt.cm.get_cmap("Blues_r")(range(0, 1, length=length(thresholds)+2))

    # loop over each roughness
    for ith = 1:length(thresholds)

        # the array where we store the number of events
        events = zeros(length(BWs))

        # loop over each beamwidth
        for ibw = 1:length(BWs)

            # calculate the acceptance at these energies
            AΩ = acceptance(ntrials, nbins;
                            trigger=gaussian_trigger(thresholds[ith], BWs[ibw]),
                            ν_min=100MHz, ν_max=1000MHz,
                            save_events=false)

            # and calculate the total number of events
            events[ibw] = sum( differential_spectrum(AΩ.energies, Apsr*AΩ.rAΩ, 2yr) )

        end

        # for this threshold, we now have a plot of beamwidth
        ax.plot(BWs, events,
                label="Threshold: $(Int(thresholds[ith] / (μV/m)))" * L" $\mu$V/m",
                color=colors[ith+1, :])
    end

    # add some labels
    ax.set(xlabel="Gaussian Beam Width [deg]",
           ylabel="Number of Events / 2yr")
    ax.legend()

    # and set a title
    # plt.title("1 GHz BW @ 20km w/ Ice Roughness")

    # add our grid and legend
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")

    fig.savefig("test/figs/beam_width_threshold_scan.png")

end

"""
Scan through payload altitude and tilt angle.
"""
function scan_altitude_tilt()

    # we use a different set of trials and bins here to make it faster
    #ntrials = 800_000
    ntrials = 1_000
    nbins = 18

    # the range of altitudes
    altitudes = [10., 20., 25., 30., 40., 50.]km

    # the tilt angles
    θ0s = [-90., -75., -60., -55., -50., -45., -40., -30.]
    θ0s = range(-90, -30, length=30)

    # create our plot
    fig, ax = plt.subplots(figsize=(4.5, 4))

    # use a colormap for the colors
    colors = plt.cm.get_cmap("Blues_r")(range(0, 1, length=length(altitudes)+1))

    # loop over the altitudes
    for i = 1:length(altitudes)

        # the array to store the altitudes
        events = zeros(length(θ0s))

        # run the acceptance simulation
        AΩ = acceptance(ntrials, nbins;
                        altitude=altitudes[i],
                        min_energy=0.5EeV, max_energy=500EeV,
                        trigger=LPDA(Nant=8, Ntrig=3, ν_max=800MHz,
                                     θ0=θ0s[1],
                                     altitude=altitudes[i]),
                        ν_min=150MHz, ν_max=800MHz,
                        save_events=true)

        # loop over the beam tilting
        for j = 1:length(θ0s)

            # construct the trigger for this simulation
            trig = LPDA(Nant=8, Ntrig=3, ν_max=800MHz,
                        θ0=θ0s[j],
                        altitude=altitudes[i])

            # retrigger the simulation to get the new effective area
            AΩeff = retrigger(AΩ, trig).rAΩ

            # calculate the number of events
            events[j] = sum(differential_spectrum(AΩ.energies, AΩeff, 2yr))

            # and print the values
            println("h=$(altitudes[i]), θ0=$(θ0s[j]) | $(events[j])")

        end

        # and plot this curve
        ax.plot(θ0s, events, label="h=$(altitudes[i])", color=colors[i, :])

    end

    # add some labels
    ax.set(xlabel="Beam Tilt Below Horizon [deg]",
           ylabel="PSR Events / 2yr")
    ax.set_xticks([-90., -75, -60, -45, -30.])
    ax.set_xlim([-90, -30])
    ax.legend()

    # and set a title
    plt.title("N=8 LPDAs")

    # add our grid and legend
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")

    fig.savefig("test/figs/altitude_tilt_scan.png")
end

"""
Compare the event rate as a function of orbit altitude.
"""
function compare_altitude()

    # we use a different set of trials and bins here to make it faster
    #ntrials = 350_000
    ntrials = 10_000
    nbins = 14

    # the altitudes that we calculate for
    altitudes = [20.0, 30.0, 40.0, 60.0, 80.0]km

    # the number of total antennas we use in each case
    Nants = [3, 8, 12, 16]
    Ntrigs = [1, 3, 5, 6]

    # use a colormap for the colors
    colors = plt.cm.get_cmap("Blues_r")(range(0, 1, length=length(Nants)+1))

    # create our plot
    fig, ax = plt.subplots(figsize=(4, 4))

    # loop over each antenna configurations
    for i = 1:length(Nants)

        # the array where we store the number of events
        events = zeros(length(altitudes))

        for j = 1:length(altitudes)

            # calculate the acceptance at these energies
            AΩ = acceptance(ntrials, nbins;
                            altitude=altitudes[j],
                            ice_depth=5.0m,
                            min_energy=0.5EeV, max_energy=500EeV,
                            trigger=LPDA(Nant=Nants[i], Ntrig=Ntrigs[i], ν_max=800MHz,
                                         θ0=-40.0, σθ=39.1, σϕ=28.4,
                                         altitude=altitudes[j]),
                            ν_min=150MHz, ν_max=800MHz,
                            save_events=false)

            # and calculate the total number of events
            events[j] = sum( differential_spectrum(AΩ.energies, AΩ.rAΩ, 2yr) )

        end

        # and plot the curve for this configuration
        ax.plot(altitudes ./ km, events, color=colors[i, :], ls="solid", label="N=$(Nants[i])")
        ax.scatter(altitudes ./ km, events, color=colors[i, :], s=2.5)
    end

    # add some labels
    ax.set(xlabel="Orbit Altitude [km]",
           ylabel="Number of Events / 2yr",
           xticks=altitudes ./ km)

    # force the bottom to be zero
    # ax.set_ylim(bottom=0.)

    # and set a title
    # plt.title("1 GHz BW @ 20km w/ Ice Roughness")

    # add our grid and legend
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")

    ax.legend(loc="upper right")

    fig.savefig("test/figs/altitude_comparison.png")

end

# function compare_altitude()

#     # we use a different set of trials and bins here to make it faster
#     ntrials = 200_000
#     nbins = 15

#     # the array where we store the number of events
#     events = zeros(length(altitudes))

#     # loop over each roughness
#     for i = 1:length(altitudes)

#         # calculate the fraction of area that is covered by PSRs
#         Apsr = psr_fraction(altitudes[i])

#         # calculate the acceptance at these energies
#         energies, dAΩ, _, rAΩ, _ = acceptance(ntrials, nbins;
#                                               min_energy=1.0EeV,
#                                               altitude=altitudes[i],
#                                               ν_max=1000MHz,
#                                               dν = 30MHz,
#                                               trigger=gaussian_trigger(sqrt(1000 / 300.) * 0.5*67μV/m),
#                                               save_events=false)

#         # and calculate the total number of events
#         events[i] = sum( differential_spectrum(energies, Apsr*rAΩ, 2yr) )

#     end

#     # create our plot
#     fig, ax = plt.subplots(figsize=(4, 4))

#     # use a colormap for the colors
#     colors = plt.cm.get_cmap("Blues_r")(range(0, 1, length=5))

#     # plot the curve
#     ax.plot(altitudes ./ km, events, color=colors[2, :], ls="solid")
#     ax.scatter(altitudes ./ km, events, color=colors[2, :], s=2.0)

#     # add some labels
#     ax.set(xlabel="Altitude [km]",
#            xlim=[0., 80.],
#            ylabel="Number of Events / 2yr")

#     # add our grid and legend
#     ax.set_axisbelow(true)
#     ax.grid(which="both", linestyle="dashed")

#     fig.savefig("test/figs/altitude_comparison.png")

# end

function compare_ice_age()

    # we use a different set of trials and bins here to make it faster
    #ntrials = 400_000
    ntrials = 1_000
    nbins = 26

    # the ice depths that we plot
    ice_depths = [5, 8, 10]m

    # the label for each depth
    ice_age = ["1 Gyr", "2 Gyr", "3 Gyr"]

    # create the figure
    fig, ax = plt.subplots(figsize=(4, 4))

    colors = plt.cm.get_cmap("Blues_r")(range(0, 1, length=5))

    # loop over each ice depth
    for i = 1:length(ice_depths)

        # calculate the acceptance
        AΩ = acceptance(ntrials, nbins;
                        altitude=20.0km,
                        ice_depth=ice_depths[i],
                        min_energy=0.5EeV, max_energy=500EeV,
                        trigger=LPDA(Nant=8, Ntrig=3, ν_max=800MHz,
                                     θ0=-30.0, σθ=28.4, σϕ=39.1,
                                     altitude=20.0km),
                        ν_min=150MHz, ν_max=800MHz,
                        save_events=false)

        # calculate the differential event spectrums
        events = differential_spectrum(AΩ.energies, AΩ.rAΩ, 2yr)

        # and plot it
        ax.hist(18.0 .+ log10.(AΩ.energies[1:end-1] ./ 1.0EeV),
                18.0 .+ log10.(AΩ.energies ./ 1.0EeV),
                weights=events,
                histtype="step", color=colors[i, :],
                label="$(ice_age[i]): N=$(Int64(round(sum(events))))")

    end

    # and finally pretty it up
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.set(xlabel=L"Energy [$\log_{10}$(eV)]",
           xlim=[17.5, 20.5],
           ylabel="Events per bin / 2yr",
           ylim=[1e0, 2e2],
           xticks=[17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5],
           yscale="log")
    ax.legend(loc="upper right")

    fig.savefig("test/figs/ice_age_comparison.png") #

end
