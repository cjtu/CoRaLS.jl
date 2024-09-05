using CoRaLS
using PyPlot
using Unitful: eV, km, sr, yr, ustrip
using Test
using LinearAlgebra

function spherical_to_cartesian(theta, phi, r)
    return (r * [sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)])
end

function plot_orbit_sampling()

    # Plot
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["axes.axisbelow"] = true
    rcParams["axes.grid"] = true
    rcParams["grid.linestyle"] = "dashed"
    rcParams["xtick.direction"] = "in"
    rcParams["ytick.direction"] = "in"
    
    fig = plt.figure(figsize=(8, 8))

    ntrials = 1000
    ax = fig.add_subplot(111, projection="3d")            
    SCs = sample_orbit(parse_orbit(), ntrials)
    # SCs[:, 3] .+= 1737.4e3  # Add in Rmoon to altitude column
    # println(SCs[])
    
    for i = 1:ntrials
        sc = spherical_to_cartesian(deg2rad(SCs[i,1]+90), deg2rad(SCs[i,2]+180), 1737.4e3 + SCs[i,3])
        sc /= norm(sc)
        ax.scatter(sc[1], sc[2], sc[3], color="k", s=20, alpha=0.5)
    end
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim([-1,1])
    ax.set_aspect("equal")
    plt.tight_layout()
    fig.savefig("$(@__DIR__)/../figs/orbit.png")

end


function plot_spherical_sampling()

    # Plot
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["axes.axisbelow"] = true
    rcParams["axes.grid"] = true
    rcParams["grid.linestyle"] = "dashed"
    rcParams["xtick.direction"] = "in"
    rcParams["ytick.direction"] = "in"
    
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    for i = 1:1000
        surface = randn(3)
        surface = surface / norm(surface)
    
        ax.scatter(surface[1], surface[2], surface[3], color="k", s=20, alpha=0.5)
    end
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim([-1,1])
    ax.set_aspect("equal")
    plt.tight_layout()
    fig.savefig("$(@__DIR__)/../figs/geom.png")

end
plot_orbit_sampling()
plot_spherical_sampling()