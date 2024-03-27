import Unitful: m, km, s, H
"""
    constants.jl

This file defines fundamental physical constants used across the `CoRaLS` module. These constants are essential for various calculations and simulations involving cosmic rays, electromagnetic fields, and other related phenomena.

## Constants
- `Rmoon`: The polar radius of the Moon, approximated as 1737.4 kilometers.
- `c_0`: The speed of light in vacuum, defined as 2.99792458 × 10^8 meters per second.
- `μ_0`: The magnetic constant (also known as the permeability of free space), valued at approximately 1.25663706212 × 10^-6 henrys per meter.
- `k_b`: Boltzmann's constant, used in statistical mechanics and thermodynamics, with a value of about 8.617333 × 10^-5 electronvolts per kelvin.
"""
# Spherical approximation to the polar radius of the Moon. From: https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
const Rmoon = 1736km

# Speed of light in vacuum
const c_0 = 2.99792458e8m / s

# Magnetic constant (permeability of free space)
const μ_0 = 1.25663706212e-6H / m

# Boltzmann's constant
const k_b = 8.617333e-5eV / K
