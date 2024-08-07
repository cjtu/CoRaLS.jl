using Dates
using DelimitedFiles
"""
    orbit.jl

Samples payload location from a given orbit.
"""

"""
    struct Orbit
        
Store orbital date and location (lon, lat, alt) information.
"""
struct Orbit
    length::Int64
    datetime::Array{Dates.DateTime, 1}
    location::Array{Float32, 2} # [longitude, latitude, altitude]
end


"""
    parse_orbit(fname)

Parse orbital info from a CSV file with columns: time, longitude, latitude, altitude.
"""
function parse_orbit(fname="lro_orbit_1yr_2010.csv", fmt=dateformat"yyyy-mm-dd HH:MM:SS.ssssss \UTC")
    Data = readdlm("$(@__DIR__)/../data/$(fname)", ',', skipstart=1)

    datetime = [Dates.DateTime(dt, fmt) for dt in Data[:, 1]]
    return Orbit(size(Data, 1), datetime, Data[:, 2:end])
end


"""
    sample_orbit(orbit, n)

Sample `n` random points from the orbit.
"""
function sample_orbit(orbit::Orbit, n::Int64)
    idx = rand(1:orbit.length, n)
    return orbit.location[idx, :]
end


# """
#     pt_in_fov(locations, latitude, longitude)

# Return only the payload locations Is point on surface vi.
# """
# function pt_in_fov()
#    return 0
# end
