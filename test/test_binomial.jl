@testset "Binomial (logistic regression)" begin

    using NonConvexPenalizedRegression
    using RData
    
    data = load(joinpath(@__DIR__, "data/Heart.rda"))
    X  = data["Heart"]["X"]
    y  = data["Heart"]["y"]

end
