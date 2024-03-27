using Test
import Unitful: km, m, s, H, eV, K

using CoRaLS

@testset verbose = true "constants.jl" begin
    @test Rmoon ≈ 1736.0km
    @test c_0 ≈ 2.99792458e8m / s
    @test μ_0 ≈ 1.25663706212e-6H / m
    @test k_b ≈ 8.617333e-5eV / K
end