module CoRaLS

using Reexport
using Unitful: g, cm, m, km, sr, μV, V, eV, GeV, EeV, Hz, MHz, W, K, NoUnits
"""
    Module CoRaLS

The `CoRaLS` (Cosmic Ray Lunar Sounder) module is designed for analyzing and simulating cosmic ray interactions, particularly in lunar contexts. It provides a comprehensive toolkit for modeling cosmic ray paths, energy spectra, and their interactions with lunar regolith and other materials.

## Usage
Typical usage of the module involves importing it and then utilizing its various functions to model cosmic ray behavior, calculate geometric properties, or simulate detection scenarios.

## Included Files and Their Purpose
- `constants.jl`: Defines a comprehensive set of constants used across the module. This may include physical constants, lunar parameters, and other universally applicable values.
- `utils.jl`: Provides utility functions that support various operations throughout the module. These might include mathematical helpers, data formatting functions, or common computational routines.
- `geometry.jl`: Contains functions related to geometric calculations and transformations, crucial for modeling and simulating spatial relationships and paths in cosmic ray studies.
- `spectrum.jl`: Offers models and functions for spectral analysis, particularly for analyzing the energy spectrum of cosmic rays and other astrophysical phenomena.
- `cosmic_ray.jl`: Dedicated to modeling and simulating cosmic ray behavior, including their generation, propagation, and interactions with different materials.
- `fresnel.jl`: Implements Fresnel equations and related calculations, typically used for understanding wave propagation and reflection.
- `ice.jl`: Focuses on the properties and behaviors of ice, especially in extraterrestrial environments, which is essential in certain cosmic ray detection scenarios.
- `sky.jl`: Contains models or functions related to the study of the sky, possibly including sky background radiation calculations or celestial body temperature models.
- `regolith.jl`: Addresses the properties of lunar regolith, such as its density and refractive index, and how these properties affect cosmic ray measurements.
- `efield.jl`: Deals with electric field calculations and models, particularly in the context of cosmic ray detection and lunar surface interactions.
- `trigger.jl`: Provides various triggering mechanisms for detecting cosmic rays or other relevant events, crucial in data acquisition systems.
- `moon.jl`: Offers functions specific to lunar studies, including calculations related to the moon's surface, environment, and its interactions with cosmic rays.
- `raytrace.jl`: Implements ray tracing algorithms, an important tool in visualizing and analyzing the paths of cosmic rays and other radiation.
- `simulate.jl`: Contains simulation tools and models for cosmic ray interactions, detector responses, and other related phenomena.
- `detector.jl`: Focuses on models and functions related to cosmic ray detectors, including their design, response, and data analysis.
- `acceptance.jl`: Provides tools for calculating the acceptance of cosmic ray detectors, an essential factor in understanding and comparing detector efficiencies.
- `surface.jl`: Includes models and functions related to surface properties, which might be critical in understanding detector-environment interactions.
- `plots.jl`: Offers a suite of plotting and visualization tools tailored for representing data and results from cosmic ray studies and related simulations.

"""

# exports from geometry.jl
export Rmoon, random_point_on_cap, horizon_angle, random_direction, spherical_cap_area, random_point_on_sphere
export WholeMoonRegion, CircularRegion, PolarRegion, CustomRegion, SouthPolePSR, NorthPolePSR, AllPSR, create_region, is_in_region, region_area
export FixedPlatform, CircularOrbit, SampledPositions, create_spacecraft, get_position

# exports from spectrum.jl
export auger_spectrum_2021, auger_spectrum_2020, auger_spectrum, sample_auger, yr, cr_spectrum

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
export acceptance, differential_spectrum, trials_passed, old_acceptance

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
