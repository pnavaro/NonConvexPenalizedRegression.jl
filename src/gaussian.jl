"""
    cdfit_gaussian( X, y, penalty, λ, eps, max_iter, γ, multiplier, α, dfmax) 

Coordinate descent for gaussian models

reference : https://github.com/pbreheny/ncvreg/blob/master/src/gaussian.c
"""
function cdfit_gaussian(X, y, penalty, λ, eps, max_iter, γ, m, α, dfmax)

    violations::Int64 = 0

    n, p = size(X)
    @assert n == length(y)
    L = length(λ)

    β = zeros(Float64, (p, L))
    loss = zeros(Float64, L)
    iter = zeros(Int64, L)

    a = zeros(Float64, p) # Beta from previous iteration

    tot_iter = 0

    r = copy(y)
    z = zeros(Float64, p)
    for j in eachindex(z)
        z[j] = crossprod(X, r, j) / n
    end

    e1 = [false for i = 1:p]
    e2 = [false for i = 1:p]

    # If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
    rss = gaussian_loss(r)

    lstart = 1

    sdy = sqrt(rss / n)

    cutoff::Float64 = 0.0

    # Path
    for l = lstart:L
        if l > 1

            # Assign a
            for j in eachindex(a)
                a[j] = β[j, l]
            end

            # Check dfmax
            nv = 0
            for j in eachindex(a)
                if (a[j] != 0)
                    nv += 1
                end
            end
            if ((nv > dfmax) || (tot_iter == max_iter))
                for ll = 1:L
                    inter[l] = -1
                end
                break
            end

            # Determine eligible set
            penalty == :Lasso && (cutoff = 2 * λ[l] - lam[l-1])
            penalty == :MCP && (cutoff = λ[l] + γ / (γ - 1) * (λ[l] - λ[l-1]))
            penalty == :SCAD && (cutoff = λ[l] + γ / (γ - 2) * (λ[l] - λ[l-1]))

            for j in eachindex(z)
                (abs(z[j]) > (cutoff * α * m[j])) && (e2[j] = 1)
            end

        else

            # Determine eligible set
            lmax = 0
            for j in eachindex(z)
                if (abs(z[j]) > lmax)
                    lmax = abs(z[j])
                end
            end

            penalty == :Lasso && (cutoff = 2 * λ[l] - lmax)
            penalty == :MCP && (cutoff = λ[l] + γ / (γ - 1) * (λ[l] - lmax))
            penalty == :SCAD && (cutoff = λ[l] + γ / (γ - 2) * (λ[l] - lmax))

            for j in eachindex(z)
                (abs(z[j]) > (cutoff * α * m[j])) && (e2[j] = 1)
            end

        end

        while (tot_iter < max_iter)

            while (tot_iter < max_iter)
                while (tot_iter < max_iter)

                    # Solve over the active set

                    iter[l] += 1
                    tot_iter += 1
                    maxChange::Float64 = 0

                    for j in eachindex(z)
                        if e1[j]

                            z[j] = crossprod(X, r, j) / n + a[j]

                            # Update β_j
                            l1 = λ[l] * m[j] * α
                            l2 = λ[l] * m[j] * (1 - α)

                            penalty == :MCP && (β[j, l] = mcp(z[j], l1, l2, γ))
                            penalty == :SCAD && (β[j, l] = scad(z[j], l1, l2, γ))
                            penalty == :Lasso && (β[j, l] = lasso(z[j], l1, l2))

                            #Update r
                            shift::Float64 = β[j, l] - a[j]

                            if shift != 0

                                for i in eachindex(r)
                                    r[i] = r[i] - shift * X[i, j]
                                end

                                (abs(shift) > maxChange) && (maxChange = abs(shift))

                            end

                        end

                    end

                    # Check for convergence
                    for j in eachindex(a)
                        a[j] = β[j, l]
                    end

                    (maxChange < eps * sdy) && break

                end

                # Scan for violations in strong set
                violations = 0

                for j in eachindex(e1)
                    if !e1[j] && e2[j]

                        z[j] = crossprod(X, r, j) / n

                        # Update β_j
                        l1 = λ[l] * m[j] * α
                        l2 = λ[l] * m[j] * (1 - α)

                        penalty == :MCP && (β[j, l] = mcp(z[j], l1, l2, γ))
                        penalty == :SCAD && (β[j, l] = scad(z[j], l1, l2, γ))
                        penalty == :Lasso && (β[j, l] = lasso(z[j], l1, l2))

                        # If something enters the eligible set, 
                        # update eligible set & residuals
                        if (β[j, l] != 0)

                            e1[j] = e2[j] = true
                            for i in eachindex(r)
                                r[i] = r[i] - β[j, l] * X[i, j]
                            end
                            a[j] = β[j, l]
                            violations += 1

                        end
                    end
                end

                violations == 0 && break

            end

            # Scan for violations in rest
            violations = 0

            for j in eachindex(e2)
                if (e2[j] == 0)

                    z[j] = crossprod(X, r, j) / n

                    # Update β_j
                    l1 = λ[l] * m[j] * α
                    l2 = λ[l] * m[j] * (1 - α)

                    penalty == :MCP && (β[j, l] = mcp(z[j], l1, l2, γ))
                    penalty == :SCAD && (β[j, l] = scad(z[j], l1, l2, γ))
                    penalty == :Lasso && (β[j, l] = lasso(z[j], l1, l2))

                    # If something enters the eligible set, update eligible set & residuals

                    if (β[j, l] != 0)

                        e1[j] = e2[j] = true
                        for i in eachindex(r)
                            r[i] -= β[j, l] * X[i, j]
                        end
                        a[j] = β[j, l]
                        violations += 1

                    end

                end

            end

            if violations == 0
                break
            end

        end

        loss[l] = gaussian_loss(r)

    end

    β, loss, iter

end

