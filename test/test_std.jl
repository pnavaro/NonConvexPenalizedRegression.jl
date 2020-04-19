@testset "std standardizes correctly" begin

import NCVREG: standardize

rng = MersenneTwister(1234)
n = 20
p = 5
l = 5

X = randn(rng,(n,p))
XX, center, scale = standardize(X)
@test all(mean(XX, dims=1)  .< fill(1e-13, p))
@test all([ c' * c for c in eachcol(XX)] .- fill(n, p) .< 1e-13 )

@test all(center .- vec(mean(X,dims=1)) .< 1e-13)
@test all(scale  .- vec(std(X,dims=1))  .< 1e-13)

end
