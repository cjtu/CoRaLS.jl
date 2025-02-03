# Modules used in acceptance.jl

The main sub-modules called by `acceptance.jl`.

## cosmic_ray.jl

Models cosmic rays and interactions. Computes depth of shower maximum (Xmax), and charge excess profile of cosmic ray-induced showers.

```@autodocs
Modules = [CoRaLS]
Pages = ["cosmic_ray.jl"]
```

## efield.jl

Model and scale the electric field from a supplied cosmic ray event lookup table file.

```@autodocs
Modules = [CoRaLS]
Pages = ["efield.jl"]
```


## fresnel.jl

Computes fresnel coefficients and wave divergence at the surface.

```@autodocs
Modules = [CoRaLS]
Pages = ["fresnel.jl"]
```

## geometry.jl

Model the geometry of the Moon and local

```@autodocs
Modules = [CoRaLS]
Pages = ["geometry.jl"]
```

## ice.jl

Models interactions of the RF signal with a variably rough ice layer.

```@autodocs
Modules = [CoRaLS]
Pages = ["ice.jl"]
```

## regolith.jl

Models the lunar regolith and propagation of RF signal through it.

```@autodocs
Modules = [CoRaLS]
Pages = ["regolith.jl"]
```

## spectrum.jl

Cosmic ray spectrum (AUGER) calculations.

```@autodocs
Modules = [CoRaLS]
Pages = ["spectrum.jl"]
```

## surface.jl

Models surface roughness and slope, including generating random surface normals and calculating surface transmission.

```@autodocs
Modules = [CoRaLS]
Pages = ["surface.jl"]
```

## trigger.jl

Simple triggers for cosmic ray detection.

```@autodocs
Modules = [CoRaLS]
Pages = ["trigger.jl"]
```

## detector.jl

Detector models for complex triggering with antennas and SNR levels.

```@autodocs
Modules = [CoRaLS]
Pages = ["detector.jl"]
```

## Sky

Computes the sky temperature for noise calculations in `detector.jl`.

```@autodocs
Modules = [CoRaLS]
Pages = ["sky.jl"]
```

## utils.jl

Miscellaneous utilities.

```@autodocs
Modules = [CoRaLS]
Pages = ["utils.jl"]
```

## plots.jl

Visualization Tools for CoRaLS Simulation Data

```@autodocs
Modules = [CoRaLS]
Pages = ["plots.jl"]
```