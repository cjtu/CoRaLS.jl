module CoRaLS

using Reexport

using Unitful: g, cm, m, km, sr, μV, V, eV, GeV, EeV, Hz, MHz, W, K, NoUnits

# exports from geometry.jl
export Rmoon, random_point_on_cap, horizon_angle, random_direction, spherical_cap_area

# exports from spectrum.jl
export auger_spectrum, yr, cr_spectrum

# exports from utils.jl
export rician, polarization_angle, retrigger

# exports from fresnel.jl
export FarFieldDivergence, NearFieldDivergence, MixedFieldDivergence

# exports from ice.jl
export NoIceRoughness, GaussianIceRoughness

# exports from sky.jl
export sky_temperature

# exports from regolith.jl
export regolith_density, regolith_index
export ConstantIndex, SurfaceDeepIndex, StrangwayIndex, ConstantDensity, StrangwayDensity

# exports from surface.jl
export NoSlope, GaussianSlope, NoRoughness, GaussianRoughness

# exports from cosmic_ray.jl
export estimate_Xmax, charge_excess_profile

# exports from efield.jl
export regolith_field, JAM, FORTE, GaisserHillasProfile, GaussianProfile, ARW

# exports from detector.jl
export LPDA, ANITA

# exports from simulate.jl
export throw_cosmicray, ScalarGeometry, RaytracedGeometry, VectorGeometry

# exports from trigger.jl
export magnitude_trigger, gaussian_trigger, rician_trigger

# exports from moon.jl
export psr_area, psr_fraction

# exports from acceptance.jl
export acceptance, differential_spectrum, trials_passed

# exports from plots.jl
export plot_differential_spectrum, plot_offaxis_angle, plot_acceptance

# import all of our code
include("constants.jl")
include("utils.jl")
include("geometry.jl")
include("spectrum.jl")
include("cosmic_ray.jl")
include("fresnel.jl")
include("ice.jl")
include("sky.jl")
include("regolith.jl")
include("efield.jl")
include("trigger.jl")
include("moon.jl")
include("raytrace.jl")
include("simulate.jl")
include("detector.jl")
include("acceptance.jl")
include("surface.jl")
include("plots.jl")

end
