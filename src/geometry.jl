using Rotations
using StaticArrays
using Distributions
using ToggleableAsserts
using Unitful
import Unitful: km

"""
Spherical approximation to the polar radius of the Moon.

Note: PG's MC had 1736.0km - TODO check.
"""
const Rmoon = 1737.4km

"""
    random_point_on_moon([n=1])

Return random latitude (-90, 90) and longitude (-180, 180) on the Moon. If n
if given, randomly draw n samples, return n x 3 array.

Samples gaussian in cartesian x,y,z, then divides by hypotenuse to put on unit 
sphere, then multiplies by Rmoon to place on the lunar surface.
"""
function random_point_on_moon()
    xyz = randn(3) 
    surface = Rmoon * xyz / norm(xyz)
    return surface
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
function cartesian_to_spherical(x, y, z)
    r = sqrt(x * x + y * y + z * z)
    return SVector{3}([acos(z / r),
        atan(y, x),
        r])
end

"""
Convert a cartesian point (x, y, z) into a lat (-90, 90), lon (-180, 180) point.
"""
function cartesian_to_latlon(x, y, z)
    theta, phi, _ = cartesian_to_spherical(x, y, z)
    φ, λ = rad2deg.([theta-pi/2, phi])  # lat (-90, 90), lon (-180, 180)
    return φ, λ
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
