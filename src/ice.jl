"""
Abstract type for ice surface roughness model.
"""
abstract type IceRoughnessModel end

"""
No ice roughness (smooth)
"""
struct NoIceRoughness <: IceRoughnessModel end

"""
Gaussian ice surface roughness model.

This is characterized by Gaussian "sigma" of the
angular distribution of the surface in degrees.
"""
struct GaussianIceRoughness <: IceRoughnessModel
    σ::typeof(0.0cm)
end

"""
    ice_roughness(::NoIceRoughness, ν, E, θ_i)

Apply roughness to a simulated electric field given the frequencies,
the electric field, and the angle of incidence at the surface. This
returns the frequencies and electric field back to the caller.

Without roughness, we just return ν and E as given.
"""
function ice_roughness(::NoIceRoughness, ν, E, θ_i)
    return ν, E
end

"""
    ice_roughness(::GaussianIceRoughness, ν, E, θ_i)

Apply roughness to a simulated electric field given the frequencies,
the electric field, and the angle of incidence at the surface. This
returns the frequencies and electric field back to the caller.

This is based on the "Diffuse reflection by rough surfaces: an introduction"
by Sylvain, Pg. 671
"""
function ice_roughness(roughness::GaussianIceRoughness, ν, E, θ_i)

    # calculate the wavelength for each frequency
    λ = c_0 ./ ν

    # calculate `g` - which is the parameter to our exponential model
    g = ( 4π .* (roughness.σ ./ λ) .* cos(θ_i) ).^2.0

    # apply the roughness and return
    return ν, E .* exp.( - g / 2.0 )
end
