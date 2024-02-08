"""
    ice.jl

Models the surface roughness of ice, which is a critical factor in simulating the reflection and transmission of radio waves on ice surfaces. It includes abstract and concrete types representing different models of ice roughness and functions to apply these models to simulated electric fields.

## Main Components

- `IceRoughnessModel`: An abstract type for general ice roughness models.

- `NoIceRoughness`: A model representing no surface roughness (smooth ice).

- `GaussianIceRoughness`: A model representing Gaussian ice surface roughness, characterized by a standard deviation of surface angles.

- `ice_roughness`: Functions to apply roughness models to electric fields, modifying them based on the specified roughness model.

"""
abstract type IceRoughnessModel end

"""
    NoIceRoughness <: IceRoughnessModel

Represents a model of ice surface with no roughness (smooth surface). This model assumes a perfectly smooth ice surface, affecting the reflection and transmission of radio waves accordingly.
"""
struct NoIceRoughness <: IceRoughnessModel end

"""
    GaussianIceRoughness(σ::typeof(0.0cm)) <: IceRoughnessModel

A model representing Gaussian ice surface roughness, characterized by a standard deviation (σ) of the angular distribution of the surface.

# Arguments
- `σ`: Standard deviation of the surface angles, representing the degree of roughness.
"""
struct GaussianIceRoughness <: IceRoughnessModel
    σ::typeof(0.0cm)
end

"""
    ice_roughness(model::IceRoughnessModel, ν, E, θ_i)

Apply the specified ice roughness model to a simulated electric field. The function modifies the electric field based on the roughness characteristics of the model.

# Arguments
- `model`: An instance of `IceRoughnessModel` to apply.
- `ν`: Frequencies of the electric field.
- `E`: Electric field amplitudes.
- `θ_i`: Angle of incidence at the ice surface.

# Returns
- Modified frequencies and electric field based on the roughness model. If the model is `NoIceRoughness`, the function returns the given frequencies and electric field. If the model is `GaussianIceRoughness`, the function modifies the electric field based on the "Diffuse reflection by rough surfaces: an introduction" by Sylvain, Pg. 671.
"""
function ice_roughness(::NoIceRoughness, ν, E, θ_i)
    return ν, E
end

function ice_roughness(roughness::GaussianIceRoughness, ν, E, θ_i)

    # calculate the wavelength for each frequency
    λ = c_0 ./ ν

    # calculate `g` - which is the parameter to our exponential model
    g = (4π .* (roughness.σ ./ λ) .* cos(θ_i)) .^ 2.0

    # apply the roughness and return
    return ν, E .* exp.(-g / 2.0)
end
