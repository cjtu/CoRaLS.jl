using Unitful: EeV, sr, ustrip
using Carlo
using HDF5

"""
Input parameters for the Monte Carlo acceptance calculation.
"""
mutable struct MC <: AbstractMC
    energy::typeof(1.0EeV)
    dcount::Int
    rcount::Int
    region::Region
    spacecraft::Spacecraft
    trigger
    kwargs::AbstractDict
end

"""
Constructor for MC. Runs on initial and restarted runs.
"""
function MC(params::AbstractDict)
    energy = params[:energy]  * EeV
    region = create_region(params[:region])
    sc = create_spacecraft(params[:spacecraft])
    trigger =  trigger_received()  #create_trigger(get(params, :trigger))
    kwargs = get(params, :kwargs, Dict())  # Simulation params
    return MC(energy, 0, 0, region, sc, trigger, kwargs)
end

"""
Initialize a Carlo run. Runs only at the start of a run (not on restarts).
"""
function Carlo.init!(mc::MC, ctx::MCContext, params::AbstractDict)
    # mc.dcount = 0
    # mc.rcount = 0
    return nothing
end

"""
Main Monte Carlo loop. Pass all params to simulation and keep track of triggers.
"""
function Carlo.sweep!(mc::MC, ctx::MCContext)
    direct, reflected = throw_cosmicray(mc.energy, mc.region, mc.spacecraft; mc.kwargs...)
    mc.dcount += mc.trigger(direct)
    mc.rcount += mc.trigger(reflected)
    return nothing
end

"""
Compute measured quantities with Carlo. Keeps track of error propagation.
"""
function Carlo.measure!(mc::MC , ctx::MCContext)
    ntrials = ctx.sweeps
    gAΩ = ustrip(pi * sr * region_area(WholeMoonRegion()))  # [km^2 sr]
    dAΩ = gAΩ * (mc.dcount / ntrials)  
    rAΩ = gAΩ * (mc.rcount / ntrials)

    measure!(ctx, :dAΩ, dAΩ)
    measure!(ctx, :rAΩ, rAΩ)

    return nothing
end

# TODO: Figure out how to make final spectrum of events a carlo evaluable (w/ error propagation)
function Carlo.register_evaluables(::Type{MC}, eval::AbstractEvaluator, params::AbstractDict)
    return nothing
end

function Carlo.write_checkpoint(mc::MC, out::HDF5.Group)
    out["dcount"] = mc.dcount
    out["rcount"] = mc.rcount
    return nothing
end

function Carlo.read_checkpoint!(mc::MC, in::HDF5.Group)
    mc.dcount = read(in, "dcount")
    mc.rcount = read(in, "rcount")
    return nothing
end
