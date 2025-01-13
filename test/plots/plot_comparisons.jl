using CoRaLS
using Unitful: km, sr, EeV, MHz, m, μV, cm, ustrip
using PyPlot

# Compare acceptance for different detector geometry
function plot_detector_geoms()
    # Setup params
    altitude = 50km
    region = AllPSR()
    sc = CircularOrbit(altitude)
    kws = Dict(:ice_depth=>6.0m, :min_energy=>1.0EeV, :max_energy=>200.0EeV,
              :ν_max=>1000MHz, :dν=>30MHz, :region=>region, :spacecraft=>sc,
              :min_count=>5, :max_tries=>1000000)
    trigger_mag = magnitude_trigger(67μV/m)
    trigger_octo = LPDA(Nant=8, Ntrig=3, ν_max=800MHz, θ0=-48.0, altitude=altitude)
    trigger_nadir = LPDA(Nant=4, Ntrig=4, ν_max=800MHz, θ0=-90.0, altitude=altitude)
    trigger_solo = LPDA(Nant=1, Ntrig=1, ν_max=800MHz, θ0=-90.0, altitude=altitude)

    # Run acceptance and collect events
    ntrials = 10000
    nbins = 8
    acc_nadir = acceptance(ntrials, nbins; trigger=trigger_nadir, save_events=true, kws...)
    @info "Saved events:", length(acc_nadir.reflected), length(acc_nadir.direct)
    @info "nadir x4 counts (r,d):", acc_nadir.rcount, acc_nadir.dcount

    # retrigger saved events on each other detector
    acc_mag = retrigger(acc_nadir, trigger_mag)
    @info "magnitude counts (r, d):", acc_mag.rcount, acc_mag.dcount
    acc_octo = retrigger(acc_nadir, trigger_octo)
    acc_solo = retrigger(acc_nadir, trigger_solo)

    # Make plots
    blues = get_cmap(:Blues_r)
    oranges = get_cmap(:Oranges_r)

    fig, ax = plt.subplots(figsize=(8, 8))
    f, ax = plot_acceptance(acc_mag; min_count=2, ax=ax, name="thresh=67μV/m", colorrefl=blues(0.1), colordirect=oranges(0.1))
    f, ax = plot_acceptance(acc_octo; min_count=2, ax=ax, name="N=8,trig=3,dip=48°", linestyle="--", colorrefl=blues(0.2), colordirect=oranges(0.2))
    f, ax = plot_acceptance(acc_nadir; min_count=2, ax=ax, name="N=4,trig=4,dip=90°", linestyle="-.", colorrefl=blues(0.4), colordirect=oranges(0.4))
    f, ax = plot_acceptance(acc_solo; min_count=2, ax=ax, name="N=1,trig=1,dip=90°", linestyle=":", colorrefl=blues(0.6), colordirect=oranges(0.6))
    ax.set_title("50km polar orbit, AllPSR, ice_depth=6m")
    fig.savefig("$(@__DIR__)/../figs/detector_geom_comparison.png")
end

plot_detector_geoms()