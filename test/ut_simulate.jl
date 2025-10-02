using CoRaLS
using Test
using LinearAlgebra

@testset "Emission Vector at Xmax" begin
    # Non-parallel geometry: antenna offset in x
    Rmoon = 1737.4 # km
    depth = 1.0 # km
    origin = [0.0, 0.0, Rmoon]
    Xmax = [0.0, 0.0, Rmoon - depth]
    antenna = [10.0, 0.0, Rmoon + 100.0] # 10 km offset in x
    axis = [0.0, 0.0, 1.0]
    view = (antenna - origin) / norm(antenna - origin)
    normal_Xmax = Xmax / norm(Xmax)
    proj_Xmax = view - (view ⋅ normal_Xmax) * normal_Xmax
    proj_Xmax /= norm(proj_Xmax)
    θ_emit = π/4 # 45 degrees
    emit = cos(θ_emit) * normal_Xmax + sin(θ_emit) * proj_Xmax
    emit /= norm(emit)
    # Should be normalized
    @test abs(norm(emit) - 1.0) < 1e-10
    # Should have a positive z component, but not parallel
    @test emit[3] < 1.0 && emit[3] > 0.5

    # Edge case: view parallel to normal_Xmax
    antenna2 = [0.0, 0.0, Rmoon + 100.0]
    view2 = (antenna2 - origin) / norm(antenna2 - origin)
    proj_Xmax2 = view2 - (view2 ⋅ normal_Xmax) * normal_Xmax
    # Should be a zero vector
    @test norm(proj_Xmax2) < 1e-10
end
