using NPZ
using Interpolations
using Unitful: g, cm, m, C, MeV, GeV, V, MHz, NoUnits

"""
A function to load and construct the ARW E-field interpolator.
"""
function load_ARW()

    # load the data file
    data = npzread("$(@__DIR__)/../data/Efield_LUT_1e20_onesided_20210724.npz")

    # construct the interpolator
    itp = LinearInterpolation((data["view_angle_array"], data["frq_list"] ./ 1e6),
                              transpose(data["RE_clean"]))

    # construct the function that does the interpolation
    return (θ, ν) -> itp(θ, (ν ./ MHz) .|> NoUnits)
end

"""
The data file for the practical and accurate calculation.
"""
const ARW_RE = load_ARW()

"""
This type is the base for all instances of electric field models
for the cosmic ray emission.
"""
abstract type FieldModel end

"""
This type is the base for all longitudinal profile models.
"""
abstract type ProfileModel end

"""
Gaussian charge excess profile.
"""
struct GaussianProfile <: ProfileModel end

"""
Gaisser-Hillas charge excess profile.
"""
struct GaisserHillasProfile <: ProfileModel end

"""
The FORTE (N. Lehtinen, PRD 2004) field model.

This requires an argument to determine what charge profile
to use - currently Gaussian or GaisserHillas
"""
struct FORTE{T<:ProfileModel} <: FieldModel end

"""
The JAM field model.

This is model is *wrong* and should never be used for production
results. The angular width is *way* too wide.

Parameterization from Alvarez-muniz et al arxiv:astro-ph/0512337 2005.
"""
struct JAM <: FieldModel end

"""
The "ARW" model from Practical & Accurate Calculations...
"""
struct ARW <: FieldModel end

"""
    regolith_field(::ARW, Ecr, θ,
                        Drego, Dvacuum;
                        n = sqrt(3.0), ν_min=50.0MHz, ν_max=300.0MHz, dν=10.0MHz,
                        Latten = 35.0m)

Calculate the integrated electric field at the payload using the ARW model.

"""
function regolith_field(::ARW, Ecr, θ,
                        Drego, Dvacuum;
                        n = sqrt(3.0), ν_min=30.0MHz, ν_max=300.0MHz,
                        density = 1.8g/cm^3, kwargs...)

    # TODO: Currently dν is fixed to 10 MHz by other parts of the code
    dν = 10MHz

    # calculate the array of frequencies that we calculate at
    ν = ν_min:dν:ν_max

    # calculate the actual Cherenkov angle for the desired efield
    θtrue = acos(1.0 / n)

    # calculate the Cherenkov angle for the simulation interpolants
    θsim = acos(1.0 / 2.02)

    # and now map (θ - θtrue) into the θsim space
    θactual = clamp(θsim + (θ - θtrue), 0., π)

    # evaluate R|E| at the θ and across a frequency band
    RE = ARW_RE(θactual, ν) * (1V/1MHz)

    # scale by the energy - the interpolator was done at 1e20 eV
    RE *= ((Ecr / 100.0EeV) |> NoUnits)

    # correct this simulation using the new shower profiles
    RE *= 1.1222

    # scale by the distance to the shower and explicitly force to
    # normal electric-field spectral units
    E = (RE ./ (Drego + Dvacuum)) .|> V/m/MHz

    # get the attenuation length at this density
    Latten = attenuation_length(ν, n, density)

    # and apply the attenuation length due to proapgation in regolith
    E .*= exp.( - Drego ./ Latten )

    return ν, E

end

"""
    regolith_field(::FORTE, Ecr, θ,
                        Drego, Dvacuum;
                        n = sqrt(3.0), ν_min=50.0MHz, ν_max=300.0MHz, dν=10.0MHz,
                        Latten = 35.0m)

Calculate the integrated electric field at the payload using the FORTE model.

This integrates a Gaisser-Hillas fit to a longitudinal charge excess profile
to determine the electric field.
"""
function regolith_field(::FORTE{GaisserHillasProfile}, Ecr, θ,
                        Drego, Dvacuum;
                        n = sqrt(3.0), ν_min=30.0MHz, ν_max=300.0MHz, dν=10.0MHz,
                        Latten = 35.0m, kwargs...)

    error("Gaisser-Hillas not implemented!")
end

"""
    regolith_field(::FORTE, Ecr, θ,
                        Drego, Dvacuum;
                        n = sqrt(3.0), ν_min=50.0MHz, ν_max=300.0MHz, dν=10.0MHz,
                        Latten = 35.0m)

Calculate the integrated electric field at the payload using the FORTE model.

This currently uses the Gaussian approxiation to the shower as presented in
the Appendix of the FORTE paper.
"""
function regolith_field(::FORTE{GaussianProfile}, Ecr, θ,
                        Drego, Dvacuum;
                        n = sqrt(3.0), ν_min=30.0MHz, ν_max=300.0MHz, dν=10.0MHz,
                        density = 1.8g/cm^3, kwargs...)

    # calculate the total charge in the shower for a particle with this energy
    Q = 0.2 * 0.62 * ((Ecr / 1GeV) |> NoUnits) * (1.602176634e-19C)

    # we assume the regolith has a relative permittivity of 1
    μ_r = 1.0

    # Assume a fixed shower length for now (1.32)
    # TODO: Incorporate the acutal shower length for this energy
    L = 1.32m

    # calculate the array of frequencies that we calculate at
    ν = ν_min:dν:ν_max

    # and the associated wave vectors
    k = 2π .* n .* ν / c_0

    # construct the functional form of the E-field
    RE = sqrt(2π) * (μ_0 .* μ_r .* Q .* L .* ν) .*
        sin(θ) .* exp.( -0.5( k .* L .* (cos(θ) .- (1.0/n))) .^ 2.0 )

    # scale by the distance to the shower and explicitly force to
    # normal electric-field spectral units
    E = (RE ./ (Drego + Dvacuum)) .|> V/m/MHz

    # get the attenuation length at this density
    Latten = attenuation_length(ν, n, density)

    # and apply the attenuation length due to proapgation in regolith
    E .*= exp.( - Drego ./ Latten )

    return ν, E

end

"""
    regolith_field(::JAM, Ecr, θ,
                        Drego, Dvacuum;
                        n = sqrt(2.8), ν_min=30.0MHz, ν_max=300.0MHz, dν=10.0MHz,
                        Latten = 35.0m)

Calculate the integrated electric field at the payload using the JAM model.

*Do not use this* - this vastly overestimates the width of the "beam",
and allows significantly more events to pass.
"""
function regolith_field(::JAM, Ecr, θ,
                        Drego, Dvacuum;
                        n = sqrt(3.0), ν_min=30.0MHz, ν_max=300.0MHz, dν=10.0MHz,
                        density = 1.8g/cm^3, kwargs...)

    # constants given in Alvarez-muniz
    X0 = 22.59g/cm^2; RM = 11.70g/cm^2
    ρ = 1.80g/cm^3
    Ec = 40.0MeV
    kL = 23.74
    kR = 1.41
    alpha = 1.23
    beta = 2.70
    kDelta = 17.78
    kE = 3.30e-16V/cm/MHz

    # the cherenkov angle in this medium
    θc = acos(1.0 / n)

    # calculate the array of frequencies that we calculate at
    ν = ν_min:dν:ν_max

    # the width of the angular distribution of the electric field (radians)
    Δθ = ( (c_0 ./ ν) * (ρ / (kDelta * X0)) .|> Unitful.NoUnits ) * (1.0 / sqrt(n*n - 1.))

    # the turnover frequency due to the longitudinal development
    # this is not currently used since we the longitudinal shower
    # development is so compact in the regolith (I think?)
    νL = (ρ / (kL * X0)) * (c_0 / abs(1.0 - n*cos(θ)))

    # this is the turnover frequency due to the radial development
    # which we consider to be dominant for regolith
    νR = (ρ / (kR * RM)) * (c_0 / sqrt(n*n - 1.0))

    # calculate the R*(electric field) at the Cherenkov angle
    REc = kE * (Ecr / Ec |> NoUnits) .* (X0 / ρ) .* (ν ./ 1MHz) .* sin(θc) .* (1.0 ./ (1.0 .+ (ν ./ νR).^alpha))

    # modulate this by the angular distribution
    RE = REc .* (sin(θ) ./ sin(θc)) .* exp.(- ((θ - θc) ./ Δθ).^2.)

    # and divide out the distance to get the efield at the payload
    E = RE ./ (Drego + Dvacuum)

    # get the attenuation length at this density
    Latten = attenuation_length(ν, n, density)

    # and apply the attenuation length due to proapgation in regolith
    E .*= exp.( - Drego ./ Latten )

    # Peter's code includes the ZHSfac - let's apply it here
    # TODO: Check if this actually needed in this formalism.
    E *= sqrt(2.) / 2.

    # and integrate this over frequency and we are done
    return ν, E

end

"""
    attenuation_length(ν, n, density)

Calculate the attenuation length at frequencies `ν` at a given density.
"""
function attenuation_length(ν, n, density; tanδnorm=6.5e-4)

    # 1e-3 is Peter's estimate for the density-normalized loss tangent
    # for polar region temperatures. Closer to 7e-4 for PSRs.
    # here we scale it up by the density
    # assuming a density of 1.27g/cm^3 at the surface.
    tanδ = tanδnorm * ( density / (1.0g/cm^3) ) |> NoUnits

    # we need the wavelength for the attenuation length calculation
    λ = c_0 ./ ν

    # since tanδ is small, we use an approximation to calculate
    # the attenuation length. I got this from PG on Slack.
    Latten = λ ./ (π * n * tanδ)

    return Latten

end
