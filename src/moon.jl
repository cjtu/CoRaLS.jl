using DelimitedFiles
using Unitful: km

# load the list of PSRs
const PSRlist = readdlm("$(@__DIR__)/../data/PSR_list.csv", ',', comments=true)

# load the North and South Pole PSR maps
# const NorthPoleMap = load("$(@__DIR__)/../data/NAC_POLE_PSR_NORTH.TIF")
# const SouthPoleMap = load("$(@__DIR__)/../data/NAC_POLE_PSR_SOUTH.TIF")

"""
    psr_area(altitude)

Calculate the area of CoRaLS observing area covered by lunar
permanently-shadowed regions.
"""
function psr_area(altitude)

    # calculate the horizon angle at this altitude
    theta = rad2deg(-horizon_angle(altitude))

    # calculate the total area of PSR viewable at this altitude
    area = sum( PSRlist[abs.(PSRlist[:, 1]) .> (90.0 - theta), 3] )km^2

    return area

end

"""
    psr_fraction(altitude)

Calculate the fraction of CoRaLS observing area covered by lunar
permanently-shadowed regions.
"""
function psr_fraction(altitude)

    # the area of the orbital band observed by CoRaLS
    total = 2.0*(CoRaLS.spherical_cap_area(π/2.0) -
        CoRaLS.spherical_cap_area(π/2. + horizon_angle(altitude)))

    return psr_area(altitude) / total

end

"""
Identify a lunar pole, North or South.
"""
abstract type LunarPole end
struct NorthPole <: LunarPole end
struct SouthPole <: LunarPole end

"""
    latlon_to_xy(::NorthPole, lat, lon)

Convert a (latitude, longitude) in degrees into (x, y) in meters
to North Pole stereographic projection, as used by the LROC NAC.

Formulas from:
http://pds.lroc.asu.edu/data/LRO-L-LROC-5-RDR-V1.0/LROLRC_2001/DOCUMENT/RDRSIS.PDF
"""
function latlon_to_xy(::NorthPole, lat, lon)
    lonp = 0.0  # central longitude in degrees
    x = 2 * Rmoon * tand(45.0 - lat/2.0) * sind(lon - lonp)
    y = -2 * Rmoon * tand(45.0 - lat/2.0) * cosd(lon - lonp)

    return x, y
end

"""
    latlon_to_xy(::NorthPole, lat, lon)

Convert a (latitude, longitude) in degrees into (x, y) in meters
to South Pole stereographic projection, as used by the LROC NAC.

Formulas from:
http://pds.lroc.asu.edu/data/LRO-L-LROC-5-RDR-V1.0/LROLRC_2001/DOCUMENT/RDRSIS.PDF
"""
function latlon_to_xy(::SouthPole, lat, lon)
    lonp = 0.0  # central longitude in degrees
    x = 2 * Rmoon * tand(45.0 + lat/2.0) * sind(lon - lonp)
    y = 2 * Rmoon * tand(45.0 + lat/2.0) * cosd(lon - lonp)

    return x, y
end

"""
    point_impacts_psr(surface)

Check whether a point on the surface impacts a PSR
"""
function point_impacts_psr(surface)

    # convert this point to (theta, phi)
    theta, phi, _ = cartesian_to_spherical(surface...)

    # calculate the latitude and longitude in this coordinate system
    lat = (π/2.0 - theta) |> rad2deg
    lon = phi |> rad2deg

    println("$(lat), $(lon)")

    # the limits of our PSR maps
    dlow = -304.00560802km
    dhigh = +304.00560802km

    # the total number of pixels
    totpix = size(CoRaLS.NorthPoleMap)[1]

    # if we are in the North Pole
    if lat > 0
        # get the coordinates in PS
        x, y = latlon_to_xy(NorthPole(), lat, lon)
    else # we are on the South Pole
        x, y = latlon_to_xy(SouthPole(), lat, lon)
    end

    println("$(x), $(y)")

    # convert it into indices into each map
    ix = Int(round(0 + (x - dlow)*((totpix - 0) / (dhigh - dlow))))
    iy = Int(round(0 + (y - dlow)*((totpix - 0) / (dhigh - dlow))))

    println("$(ix), $(iy)")

    # and finally check if the map has any PSRs there
    if lat > 0
        return NorthPoleMap[ix, iy] > 0.0
    else
        return SouthPoleMap[ix, iy] > 0.0
    end

end
