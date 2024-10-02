using Dates
using DelimitedFiles
"""
    orbit.jl

Samples payload location from a given orbit.
"""
abstract type AbstractOrbit end
struct CircularOrbit <: AbstractOrbit end
"""
    struct Orbit
        
Store orbital date and location (lon, lat, alt) information.
"""
struct Orbit <: AbstractOrbit
    length::Int64
    datetime::Array{Dates.DateTime, 1}
    location::Array{Float32, 2} # [longitude, latitude, altitude]
end


"""
    parse_orbit(fname)

Parse orbital info from a CSV file with columns: time, longitude, latitude, altitude.

Time column is parsed with datetime format given in fmt.
"""
function parse_orbit(fname="lro_orbit_1yr_2010.csv", fmt=dateformat"yyyy-mm-dd HH:MM:SS.ssssss \UTC")
    Data = readdlm("$(@__DIR__)/../data/$(fname)", ',', skipstart=1)

    datetime = [Dates.DateTime(dt, fmt) for dt in Data[:, 1]]
    return Orbit(size(Data, 1), datetime, Data[:, 2:end])
end


"""
    sample_orbit(orbit::CircularOrbit, n::Integer; altitude::Union{Real, Unitful.Length} = 20km)

Sample `n` random positions on a circular orbit around a celestial body.

# Arguments
- `orbit::CircularOrbit`: The circular orbit to sample from.
- `n::Integer`: The number of positions to sample.
- `altitude::Union{Real, Unitful.Length} = 20km`: The altitude of the orbit. Defaults to 20km.

# Returns
- `Matrix{Float64}`: A 3×n matrix where each column represents a position [φ, λ, altitude].
  φ is the latitude in radians, λ is the longitude in radians, and altitude is in meters.
```
"""
function sample_orbit(orbit::CircularOrbit, n::Integer; altitude::Union{Real, Unitful.Length} = 20km)
    # Convert altitude to km if it's a Unitful.Length
    alt_m = altitude isa Unitful.Length ? ustrip(u"km", altitude) : Float64(altitude)

    # Sample uniform distributions for latitude and longitude
    φ = rad2deg.(rand(Uniform(0, π), n))  # latitude
    λ = rad2deg.(rand(Uniform(0, 2π), n))      # longitude

    # Create and return the matrix of positions
    return [φ λ fill(alt_m, 1, n)']
end

"""
 sample_orbit(orbit::Orbit, n::Integer)

Sample `n` random points from a pre-defined orbit.

# Arguments
- `orbit::Orbit`: The orbit to sample from. Must have `length` and `location` fields.
- `n::Integer`: The number of positions to sample.

# Returns
- `Matrix`: A subset of `orbit.location` containing `n` randomly sampled rows.
```
"""
function sample_orbit(orbit::Orbit, n::Integer)
    idx = rand(1:orbit.length, n)
    return orbit.location[idx, :]
end

