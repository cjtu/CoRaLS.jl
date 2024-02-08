# Modules used in acceptance.jl

Several modules are utilized by `acceptance.jl` to perform its tasks. This section provides documentation for these modules.

## Constants

- `c_0`: The speed of light in vacuum, defined as 2.99792458 × 10^8 meters per second.
- `μ_0`: The magnetic constant (also known as the permeability of free space), valued at approximately 1.25663706212 × 10^-6 henrys per meter.
- `k_b`: Boltzmann's constant, used in statistical mechanics and thermodynamics, with a value of about 8.617333 × 10^-5 electronvolts per kelvin.

## cosmic_ray.jl

This file deals with modeling and analysis of cosmic rays. It includes functions to estimate key parameters of cosmic ray interactions, such as the depth of shower maximum (Xmax), and to propagate cosmic rays to their point of maximum interaction. The file also provides tools to model the charge excess profile of cosmic ray-induced air showers.

```@docs
CoRaLS.estimate_Xmax
CoRaLS.propagate_to_Xmax
CoRaLS.charge_excess_profile
```

## fresnel.jl

Fresnel reflections and transmissions related to cosmic ray detection. It includes abstract types for divergence models and specific implementations for different divergence scenarios. It also contains functions for calculating Fresnel reflection and transmission coefficients and their modifications for spherical wave divergence at the surface.

```@docs
CoRaLS.DivergenceModel
CoRaLS.FarFieldDivergence
CoRaLS.fresnel_rpar
CoRaLS.fresnel_tpar
CoRaLS.divergence_tperp
```

## geometry.jl

Contains functions for modeling the geometry of cosmic ray interactions with the lunar surface. This includes functions to calculate the horizon angle, random vectors, and spherical cap areas, as well as functions to convert between spherical and Cartesian coordinates and to intersect rays with spheres.

```@docs
CoRaLS.Rmoon
CoRaLS.random_north_pole_point
CoRaLS.random_south_pole_point
CoRaLS.random_angles_on_cap
CoRaLS.random_point_on_cap
CoRaLS.spherical_to_cartesian
CoRaLS.cartesian_to_spherical
CoRaLS.horizon_angle
CoRaLS.random_vector
CoRaLS.random_direction
CoRaLS.spherical_cap_area
CoRaLS.intersect_with_sphere
CoRaLS.propagate_and_refract
```

## ice.jl

Models the surface roughness of ice, which is a critical factor in simulating the reflection and transmission of radio waves on ice surfaces. It includes abstract and concrete types representing different models of ice roughness and functions to apply these models to simulated electric fields.

```@docs
CoRaLS.IceRoughnessModel
CoRaLS.NoIceRoughness
CoRaLS.GaussianIceRoughness
CoRaLS.ice_roughness
CoRaLS.ice_roughness
```

## moon.jl

Contains functionality related to lunar geography, specifically focusing on permanently shadowed regions (PSRs) on the moon's surface. It includes functions to calculate the visible area of PSRs from a given altitude, convert geographical coordinates to stereographic projections, and determine if a surface point impacts a PSR.

```@docs
CoRaLS.psr_area
CoRaLS.psr_fraction
CoRaLS.latlon_to_xy
CoRaLS.point_impacts_psr
```

## plots.jl

Visualization Tools for CoRaLS Simulation Data

```@docs
CoRaLS.plot_incident_angles
CoRaLS.plot_polarization_angle
CoRaLS.plot_offaxis_angle
CoRaLS.plot_acceptance
CoRaLS.plot_differential_spectrum

```

## regolith.jl

Contains functions for modeling the regolith of the lunar surface. This includes abstract and concrete types representing different models of regolith density and functions to apply these models to simulated electric fields.

```@docs
CoRaLS.regolith_density
CoRaLS.create_density_lut
CoRaLS.regolith_index
CoRaLS.RegolithIndex
CoRaLS.ConstantIndex
CoRaLS.SurfaceDeepIndex
CoRaLS.StrangwayIndex
CoRaLS.RegolithDensity
CoRaLS.ConstantDensity
CoRaLS.StrangwayDensity
CoRaLS.OldIncorrectDensity

```

## spectrum.jl

```@docs
CoRaLS.sample_power_law
CoRaLS.auger_spectrum
CoRaLS.sample_auger
```

## surface.jl

Contains functions for modeling surface roughness and slope. This includes the Gaussian roughness and slope models, as well as functions to generate random surface normals and calculate surface transmission.

```@docs
CoRaLS.RoughnessModel
CoRaLS.NoRoughness
CoRaLS.GaussianRoughness
CoRaLS.SlopeModel
CoRaLS.NoSlope
CoRaLS.GaussianSlope
CoRaLS.random_surface_normal
CoRaLS.surface_transmission
```

## trigger.jl

Contains functions describing the trigger mechanisms for cosmic ray detection. This includes the magnitude trigger, Gaussian trigger, and Rician trigger.

```@docs
CoRaLS.magnitude_trigger
CoRaLS.gaussian_trigger
CoRaLS.rician_trigger
```
