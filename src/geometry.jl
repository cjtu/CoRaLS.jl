using Rotations
using StaticArrays
using Distributions
using ToggleableAsserts
using Unitful
import Unitful: km

"""
Spherical approximation to the polar radius of the Moon.
"""
const Rmoon = 1737.4km

# Define Moon region types and functions
abstract type Region end

struct WholeMoonRegion <: Region end

struct CircularRegion <: Region
    center_lat::Float64
    center_lon::Float64
    radius::Float64
end

struct PolarRegion <: Region
    pole::Symbol  # :north or :south
    angle::Float64  # degrees from pole
end

struct CustomRegion <: Region
    criteria::Function
    area::Float64  # must be in [km^2]
end

# PSR Definitions from the following paper giving PSR area in south and north.
# Ref: Mazarico, E., G. A. Neumann, D. E. Smith, M. T. Zuber, and M. H. Torrence (2011) “Illumination conditions of the lunar polar regions using LOLA topography.” Icarus, 211, no. 2: 1066-1081. https://doi.org/10.1016/j.icarus.2010.10.030 
# Chance of hitting PSR = area of PSR / area 10 of spherical cap 10 degrees from pole (equiprobable for C.R. to hit anywhere on pole)
#  TODO: to improve accuracy, could update this to lookup table of lat/lons in PSRs

SouthPolePSR = CustomRegion(
    (lat, lon) -> (lat <= -80) && rand() < 0.0557,  # 5.57% chance in PSR
    1.6055e4  # [km^2] Area of South pole PSRs (Mazarico et al. 2011)
)

NorthPolePSR = CustomRegion(
    (lat, lon) -> (lat >= 80) && rand() < 0.0446,  # 4.46% chance in PSR
    1.2866e4  # [km^2] Area of North pole PSRs (Mazarico et al. 2011)
)

AllPSR = CustomRegion(
    (lat, lon) -> (abs(lat) >= 80) && rand() < 0.0502,  # 5.02% chance in PSR
    2.8921e4  # [km^2] Total PSR area (Mazarico et al. 2011)
)

function create_region(config::String)
    if config == "whole_moon"
        return WholeMoonRegion()
    elseif startswith(config, "circular:")
        lat, lon, radius = parse.(Float64, split(config[10:end], ","))
        return CircularRegion(lat, lon, radius)
    elseif startswith(config, "polar:")
        pole, angle = split(config[7:end], ",")
        return PolarRegion(Symbol(pole), parse(Float64, angle))
    else
        throw(ArgumentError("Invalid region configuration"))
    end
end

# region location filtering
function is_in_region(surface, region::WholeMoonRegion)
    return true  # whole moon includes all surface points
end

function is_in_region(surface, region::CircularRegion)
    lat, lon = deg2rad.(cartesian_to_latlon(surface))
    center_lat, center_lon = deg2rad.(region.center_lat, region.center_lon)
    
    dlat = (lat - center_lat)
    dlon = (lon - center_lon)
    a = sin(dlat/2)^2 + cos(lat) * cos(center_lat) * sin(dlon/2)^2
    gc = 2 * atan(sqrt(a), sqrt(1-a)) * Rmoon
    
    return gc <= region.radius
end

function is_in_region(surface, region::PolarRegion)
    lat, _ = cartesian_to_latlon(surface)
    if region.pole == :north
        return lat >= 90 - region.angle
    else
        return lat <= -90 + region.angle
    end
end

function is_in_region(surface, region::CustomRegion)
    lat, lon = cartesian_to_latlon(surface)
    return region.criteria(lat, lon)
end

# Return area of region
function region_area(region::WholeMoonRegion)
    return 4 * π * Rmoon^2
end

function region_area(region::CircularRegion)
    return 2 * π * Rmoon^2 * (1 - cos(region.radius / Rmoon))
end

function region_area(region::PolarRegion)
    return 2 * π * Rmoon^2 * (1 - cos(region.angle * π / 180))
end

function region_area(region::CustomRegion)
    return region.area * km^2
end

# Define spacecraft types and location sampling
abstract type Spacecraft end

struct FixedPlatform <: Spacecraft
    lat::Float64
    lon::Float64
    altitude::typeof(1.0km)
end

struct CircularOrbit <: Spacecraft 
    altitude::typeof(1.0km)
end

struct SampledOrbit <: Spacecraft
    positions
end

function create_spacecraft(config::String)
    if startswith(config, "fixed:")
        lat, lon, alt = parse.(Float64, split(config[7:end], ","))
        return FixedPlatform(lat, lon, alt * km)
    elseif startswith(config, "orbit:")
        alt = parse(Float64, split(config, ":")[2])
        return CircularOrbit(alt * km)  # TODO: strip km from end if given
    elseif startswith(config, "file:")
        filename = config[6:end]
        return parse_orbit(filename)
    else
        throw(ArgumentError("Invalid spacecraft configuration"))
    end
end

"""
    parse_orbit(fname)

Parse orbital info from a CSV file with columns: time, longitude, latitude, altitude.

Returns SampledOrbit matrix of positions to randomly sample.
"""
function parse_orbit(fname="lro_orbit_1yr_2010.csv")
    data = readdlm("$(@__DIR__)/../data/$(fname)", ',', skipstart=1)
    return SampledOrbit(data[:, 2:end])
end

# Not yet implemented - no use for orbital datetimes yet
# function parse_orbit(fname="lro_orbit_1yr_2010.csv", fmt=dateformat"yyyy-mm-dd HH:MM:SS.ssssss \UTC")
#     data = readdlm("$(@__DIR__)/../data/$(fname)", ',', skipstart=1)
#     datetime = [Dates.DateTime(dt, fmt) for dt in data[:, 1]]
#     return SampledOrbit(data[:, 2:end])
# end

function get_position(spacecraft::FixedPlatform)
    return latlon_to_cartesian(spacecraft.lat, spacecraft.lon, Rmoon + spacecraft.altitude)
end

function get_position(spacecraft::CircularOrbit)
    return random_point_on_cicular_orbit(Rmoon + spacecraft.altitude)
end

function get_position(spacecraft::SampledOrbit)
    idx = rand(1:size(spacecraft.positions)[1])
    lon, lat, alt = spacecraft.positions[idx, :]  # [deg, deg, km]
    return latlon_to_cartesian(lat, lon, Rmoon + alt*km)
end

function is_visible(surface, spacecraft)
    surface_normal = normalize(surface)  
    to_spacecraft = normalize(spacecraft - surface)
    return dot(surface_normal, to_spacecraft)  > 0  # central angle < 90 deg
end

"""
    random_point_on_sphere()

Return cartesian point (x, y, z) randomly distributed on sphere of radius r.

This uses the triple gaussian method which is slightly faster than the 
spherical coords method. https://mathworld.wolfram.com/SpherePointPicking.html
"""
function random_point_on_sphere(r=1)
    xyz = SVector{3, Float64}(randn(), randn(), randn())
    return r * xyz / norm(xyz)
end

"""
Return cartesian vector to spacecraft at altitude.

Note: this samples uniformly in latitude which does not give uniformly
distributed points on the sphere (overdensity of points on the poles which 
is expected for circular polar orbits).
"""
function random_point_on_cicular_orbit(r=Rmoon)
    φ = π * rand()  # Co-latitute uniform on [0, π]
    λ = 2π * rand() # Longitude uniform on [0, 2π]
    return spherical_to_cartesian(φ, λ, r)
end
"""
    random_north_pole_point()

Draw a random point on the North lunar pole
"""
random_north_pole_point() = random_point_on_cap(deg2rad(10.0); r=Rmoon)

"""
    random_south_pole_point()

Draw a random point on the South lunar pole
"""
function random_south_pole_point()
    # sample the (theta, phi) angles on the North pole
    theta, phi = random_angles_on_cap(deg2rad(10.0))

    # rotate this to be at the south pole
    theta = π - theta

    # and return them in Cartesian coordinates
    return spherical_to_cartesian(theta, phi, Rmoon)
end

"""
Generate a (theta, phi) pair on a spherical cap uniformly in solid angle.

All angles are in radians.
"""
function random_angles_on_cap(theta_max; theta_min=0.0, phi_min=0.0, phi_max=2 * pi)

    # draw a random polar angle uniformly in solid angle
    theta = acos(rand(Uniform(cos(theta_max), cos(theta_min))))

    # draw an azimuthal angle uniformaly in the desired range
    phi = rand(Uniform(phi_min, phi_max))

    return theta, phi

end

"""
Sample a random point on a spherical cap uniformly in solid angle.

We assume that the spherical cap is centered on the +z-pole of the sphere.
"""
function random_point_on_cap(theta_max; r=Rmoon, kwargs...)

    # sample the (theta, phi) angles on the cap
    theta, phi = random_angles_on_cap(theta_max; kwargs...)

    # and return them in Cartesian coordinates
    return spherical_to_cartesian(theta, phi, r)

end


"""
Convert a spherical point (theta, phi, r) into Cartesian coordinates.
"""
function spherical_to_cartesian(theta, phi, r)
    return SVector{3}(r * [sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)])
end

"""
Convert a cartesian point (x, y, z) into a spherical point (theta, phi, r).
"""
function cartesian_to_spherical(point::SVector{3})
    x, y, z = point
    r = sqrt(x * x + y * y + z * z)
    return SVector{3}([acos(z / r), atan(y, x), r])
end

"""
Convert a cartesian point (x, y, z) into (lat, lon)).
"""
function cartesian_to_latlon(point::SVector{3})
    x, y, z = point
    lat = rad2deg(asin(z / norm(point)))
    lon = rad2deg(atan(y, x))
    return lat, lon
end

"""
Convert a (lat, lon) to cartesian (x, y, z) vector at radius r.
"""
function latlon_to_cartesian(lat::Float64, lon::Float64)
    lat_rad = deg2rad(lat)
    lon_rad = deg2rad(lon)
    x = cos(lat_rad) * cos(lon_rad)
    y = cos(lat_rad) * sin(lon_rad)
    z = sin(lat_rad)
    return SVector{3}(x, y, z)
end

function latlon_to_cartesian(lat::Float64, lon::Float64, r::Union{Float64, typeof(1.0km)})
    return r * latlon_to_cartesian(lat, lon)
end

"""
    horizon_angle(altitude; R=Rmoon)

Calculate the angle below the horizontal of the horizon from a given altitude.

This returns "negative" angles, in radians, for "below" the horizontal.

# Arguments
- `height`: The altitude above the surface (with units).`

"""
function horizon_angle(altitude; R=Rmoon)
    return -acos(R / (R + altitude))
end

"""
    random_vector()

Generate a random Cartesian direction distributed uniformly on the unit-sphere
with the pole of the sphere aligned with `normal`
"""
function random_vector(normal=SA[0.0, 0.0, 1.0])

    # check that normal is accurately normalized
    @toggled_assert norm(normal) ≈ 1.0

    # draw a random polar angle uniformly in solid angle
    theta = acos(rand(Uniform(0.0, 1.0)))

    # draw an azimuthal angle uniformaly in the desired range
    phi = rand(Uniform(0.0, 2π))

    # this is the vector in the z-hat space
    direction = spherical_to_cartesian(theta + π / 2.0, phi, 1.0)

    # we must now rotate this so that z-hat is aligned with normal
    # we do this with an axis angle representation

    # this is our z-hat vector in the original coordinate system
    zhat = SA[0.0, 0.0, 1.0]

    # construct the angle between z-hat and normal - already normalized
    θ = acos(zhat ⋅ normal)

    # and construct the *right-handed* axis
    axis = normal × zhat

    # have to check that axis can be properly normalized
    if norm(axis) < 1e-6
        # if this is true, we are in the same coordinate system
        return direction
    end

    # otherwide, normalize, do the rotation, and return
    return AngleAxis(-θ, (axis / norm(axis))...) * direction


end

"""
    random_direction()

Generate a random Cartesian direction distributed uniformly on the unit-sphere
with the pole of the sphere aligned with `normal` weighted by cos^2(\theta)
"""
function random_direction(normal=SA[0.0, 0.0, 1.0])

    # check that normal is accurately normalized
    @toggled_assert norm(normal) ≈ 1.0

    # draw a random polar angle uniformly in solid angle
    theta = acos(sqrt(rand(Uniform(0.0, 1.0))))

    # draw an azimuthal angle uniformaly in the desired range
    phi = rand(Uniform(0.0, 2π))

    # this is the vector in the z-hat space
    direction = spherical_to_cartesian(theta + π / 2.0, phi, 1.0)

    # we must now rotate this so that z-hat is aligned with normal
    # we do this with an axis angle representation

    # this is our z-hat vector in the original coordinate system
    zhat = SA[0.0, 0.0, 1.0]

    # construct the angle between z-hat and normal - already normalized
    θ = acos(zhat ⋅ normal)

    # and construct the *right-handed* axis
    axis = normal × zhat

    # have to check that axis can be properly normalized
    if norm(axis) < 1e-6
        # if this is true, we are in the same coordinate system
        return direction
    end

    # otherwide, normalize, do the rotation, and return
    return AngleAxis(-θ, (axis / norm(axis))...) * direction


end

"""
    spherical_cap_area(theta, r=Rmoon)

Calculate the area of a spherical cap with central angle `theta`.

Returned units will be the square of the units of `r`.
"""
function spherical_cap_area(theta, r=Rmoon)
    return 2π * r * r * (1.0 - cos(theta))
end


"""
    intersect_with_sphere(start, direction, radius)

Propagate a vector from `start` to `radius` along `direction`.

Assumes `direction` is normalized.

See the below link for the formalism behind this implementation:
    https://www.scratchapixel.com/lessons/3d-basic-rendering/
        minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection

"""
function intersect_with_sphere(start, direction, radius)

    # the distance from the origin to the perpendicular radial line
    LD = start ⋅ direction

    # the squared distance from the origin that we intersect that line
    d2 = (start ⋅ start) - LD * LD

    # if d^2 is greater than radius^2, then we never intersect
    if d2 > radius * radius
        return nothing
    end

    # compute the offset from the central intersection
    thc = sqrt(radius * radius - d2)

    # there are two solutions. -LD +/- thc

    #  -LD + thc is the first solution
    #  -LD - thc is the second solution

    # when we are INSIDE the sphere, we want the far/second solution
    # as this is the solution along the ray path.
    # when we are OUTSIDE the sphere, we want the closest/first solution
    # as that is the location of first intersection

    # we use the radius of the starting location to determine which solution
    # we pick so this method works from inside and outside
    t = -LD - sign(norm(start) - radius) * thc

    # `t` is now the length alongth `direction` that we must step to reach the intersection
    intersection = start + t * direction

    # a sanity check that we are actually at the correct radius
    @toggled_assert (norm(intersection) - radius) / radius < 1e-6 "Radius of intersection is invalid"

    return intersection


end

"""
    propagate_and_refract(start, antenna, direction, scale, indexmodel)

Propagate from `start` along `direction` to `antenna`, shifting `direction`
by scale along the normal using the `indexmodel` for the refractive index.

This function is for use during minimization of the refractive ray path
and is not intended for public use (the API may change).
"""
function propagate_and_refract(start, direction, antenna, scale, indexmodel)

    # scale the z-coordinate to change the incident angle at the surface
    emit = direction + (scale[1] .* (start / norm(start)))
    emit /= norm(emit)

    # propagate it to the surface
    surface = intersect_with_sphere(start, emit, Rmoon)

    # we need the surface normal, not the actual surface location
    normal = surface / norm(surface)

    # calculate the incident angle at the surface - they are already normalized
    θ_i = acos(normal ⋅ emit)

    # check for TIR - send wave along the surface
    if θ_i > asin(1.0 / regolith_index(indexmodel, 0.0m))
        # return emit, surface, nothing, nothing, 0.
        θ_r = pi / 2.0
    else
        # use Snell's law to calculate the actual refracted angle
        θ_r = asin(regolith_index(indexmodel, 0.0m) * sin(θ_i))
    end

    # and calculate the actual refracted angle to the payload
    θ_true = acos(((antenna - surface) / norm(antenna - surface)) ⋅ normal)

    # and return the residuals between the two
    return emit, surface, θ_i, θ_r, abs(θ_r - θ_true)

end
