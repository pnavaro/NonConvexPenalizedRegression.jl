"""
    gaussian_loss( r )

Gaussian loss 
"""
function gaussian_loss(r::Vector{Float64})::Float64
    l = 0.0
    for i in eachindex(r)
        l += r[i]^2
    end
    l
end

