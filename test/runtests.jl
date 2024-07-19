using Aqua
using NonConvexPenalizedRegression
using Random
using Statistics
using Test

include("test_gaussian.jl")
include("test_binomial.jl")
include("test_std.jl")
include("test_scad.jl")

@testset "Aqua.jl" begin
  Aqua.test_all(NonConvexPenalizedRegression)
end
