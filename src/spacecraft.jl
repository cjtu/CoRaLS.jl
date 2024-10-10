"""
    Module spacecraft

Defines the orbital platform observation geometry.
"""
abstract type Spacecraft end

struct FixedPlatform <: Spacecraft
    lat::Float64
    lon::Float64
end

struct CircularOrbit <: Spacecraft end

struct RandomPositions <: Spacecraft
    positions::Vector{SVector{3,Float64}}
end

function create_spacecraft(config::String)
    if startswith(config, "fixed:")
        lat, lon = parse.(Float64, split(config[7:end], ","))
        return FixedPlatform(lat, lon)
    elseif config == "orbit"
        return CircularOrbit()
    elseif startswith(config, "file:")
        filename = config[6:end]
        positions = [SVector{3,Float64}(parse.(Float64, split(line))) for line in eachline(filename)]
        return RandomPositions(positions)
    else
        throw(ArgumentError("Invalid spacecraft configuration"))
    end
end

function get_position(spacecraft::FixedPlatform, altitude::Float64)
    direction = latlon_to_direction(spacecraft.lat, spacecraft.lon)
    return direction * (Rmoon + altitude)
end

function get_position(spacecraft::CircularOrbit, altitude::Float64)
    return random_orbit_position(altitude)
end

function get_position(spacecraft::RandomPositions, altitude::Float64)
    position = rand(spacecraft.positions)
    return normalize(position) * (Rmoon + altitude)
end
