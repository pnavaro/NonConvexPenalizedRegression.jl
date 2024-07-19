function reg_binomial_scad(X, y)

    gamma = 3.7
    XX, center, scale = standardize(X)


end

function p_binomial(eta)
    if eta > 10
        return 1.0
    elseif eta < -10
        return 0.0
    else
        return exp(eta) / (1 + exp(eta))
    end
end


fmax2(x, y) = (x < y) ? y : x

"""
$(SIGNATURES)

Coordinate descent for binomial models
"""
function cdfit_glm(
    X,
    y,
    family,
    penalty,
    lambda,
    eps,
    max_iter,
    gamma,
    multiplier,
    alpha,
    dfmax,
    user,
    warn,
)

    # Lengths/dimensions
    n = length(y)
    p = length(X) / n
    L = length(lambda)

    lam = lambda
    m = multiplier
    tot_iter = 0

    # Outcome
    beta0 = zeros(L)
    b0 = beta0
    beta = zeros(L, p)
    b = beta
    Dev = zeros(L)
    Eta = zeros(L * n)
    iter = zeros(Int, L)

    # Intermediate quantities
    a0 = 0          # Beta0 from previous iteration
    a = zeros(p)    # Beta from previous iteration
    r = zeros(n)
    w = zeros(n)
    s = zeros(n)
    z = zeros(p)
    eta = zeros(n)
    e1 = zeros(Int, p)
    e2 = zeros(Int, p)

    # Initialization
    ybar = sum(y) / n
    a0 = b0[1] = log(ybar / (1 - ybar))
    nullDev = 0.0
    if family == "binomial"
        a0 = b0[1] = log(ybar / (1 - ybar))
        for i in eachindex(y)
            nullDev -= 2 * y[i] * log(ybar) + 2 * (1 - y[i]) * log(1 - ybar)
        end
    elseif family == "poisson"
        a0 = b0[1] = log(ybar)
        for i in eachindex(y)
            if y[i] != 0
                nullDev += 2 * (y[i] * log(y[i] / ybar) + ybar - y[i])
            else
                nullDev += 2 * ybar
            end
        end
    end
    s .= y .- ybar
    eta .= a0
    z .= crossprod(X, s, n, j) ./ n

    # If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
    if user
        lstart = 0
    else
        lstart = 1
        Dev[1] = nullDev
        Eta .= eta
    end

    # Path
    for l = lstart+1:L

        if l != 0
            # Assign a, a0
            a0 = b0[l]
            for j in eachindex(a)
                a[j] = b[j, l]
            end

            # Check dfmax
            nv = 0
            for j in eachindex(a)
                if a[j] != 0
                    nv += 1
                end
            end
            if (nv > dfmax) || (tot_iter == max_iter)
                for ll = l+1:L
                    iter[ll] = missing
                end
                break
            end

            # Determine eligible set
            if penalty == "lasso"
                cutoff = 2 * lam[l] - lam[l-1]
            end
            if penalty == "MCP"
                cutoff = lam[l] + gamma / (gamma - 1) * (lam[l] - lam[l-1])
            end
            if penalty == "SCAD"
                cutoff = lam[l] + gamma / (gamma - 2) * (lam[l] - lam[l-1])
            end
            for j in eachindex(z)
                if fabs(z[j]) > (cutoff * alpha * m[j])
                    e2[j] = 1
                end
            end
        else

            # Determine eligible set
            lmax = 0.0
            for j in eachindex(z)
                if (fabs(z[j]) > lmax)
                    lmax = fabs(z[j])
                end
            end
            if penalty == "lasso"
                cutoff = 2 * lam[l] - lmax
            end
            if penalty == "MCP"
                cutoff = lam[l] + gamma / (gamma - 1) * (lam[l] - lmax)
            end
            if penalty == "SCAD"
                cutoff = lam[l] + gamma / (gamma - 2) * (lam[l] - lmax)
            end
            for j in eachindex(z)
                if fabs(z[j]) > (cutoff * alpha * m[j])
                    e2[j] = 1
                end
            end

        end

        while tot_iter < max_iter
            while tot_iter < max_iter
                while tot_iter < max_iter
                    iter[l] += 1
                    tot_iter += 1

                    # Approximate L
                    REAL(Dev)[l] = 0
                    if family == "binomial"
                        v = 0.25
                        for i in eachindex(eta)
                            mu = p_binomial(eta[i])
                            w[i] = fmax2(mu * (1 - mu), 0.0001)
                            s[i] = y[i] - mu
                            r[i] = s[i] / w[i]
                            if y[i] == 1
                                Dev[l] = Dev[l] - log(mu)
                            end
                            if y[i] == 0
                                Dev[l] = Dev[l] - log(1 - mu)
                            end
                        end
                    elseif family == "poisson"
                        for i in eachindex(w)
                            mu = exp(eta[i])
                            w[i] = mu
                            s[i] = y[i] - mu
                            r[i] = s[i] / w[i]
                            if y[i] != 0
                                Dev[l] += y[i] * log(y[i] / mu)
                            end
                        end
                    end

                    # Check for saturation
                    if (REAL(Dev)[l] / nullDev < 0.01)
                        warn && @warn("Model saturated; exiting...")
                        for ll = l+1:L
                            iter[ll] = missing
                        end
                        tot_iter = max_iter
                        break
                    end

                    # Intercept
                    xwr = crossprod(w, r, 0)
                    xwx = sum(w, n)
                    b0[l] = xwr / xwx + a0
                    for i in eachindex(r)
                        si = b0[l] - a0
                        r[i] -= si
                        eta[i] += si
                    end
                    max_change = fabs(si) * xwx / n

                    # Covariates
                    for j in eachindex(e1)

                        if e1[j]

                            # Calculate u, v
                            xwr = wcrossprod(X, r, w, j)
                            xwx = wsqsum(X, w, j)
                            u = xwr / n + (xwx / n) * a[j]
                            v = xwx / n

                            # Update b_j
                            l1 = lam[l] * m[j] * alpha
                            l2 = lam[l] * m[j] * (1 - alpha)
                            if penalty == "MCP"
                                b[l*p+j] = MCP(u, l1, l2, gamma, v)
                            end
                            if penalty == "SCAD"
                                b[l*p+j] = SCAD(u, l1, l2, gamma, v)
                            end
                            if penalty == "lasso"
                                b[l*p+j] = lasso(u, l1, l2, v)
                            end

                            # Update r
                            shift = b[l*p+j] - a[j]
                            if shift != 0
                                for i in eachindex(r)
                                    si = shift * X[j*n+i]
                                    r[i] -= si
                                    eta[i] += si
                                end
                                if fabs(shift) * sqrt(v) > max_change
                                    max_change = fabs(shift) * sqrt(v)
                                end
                            end
                        end
                    end

                    # Check for convergence
                    a0 = b0[l]
                    for j in eachindex(a)
                        a[j] = b[j, l]
                    end
                    if max_change < eps
                        break
                    end
                end

                # Scan for violations in strong set
                violations = 0
                for j in eachindex(z)
                    if e1[j] == 0 && e2[j] == 1
                        z[j] = crossprod(X, s, j) / n
                        l1 = lam[l] * m[j] * alpha
                        if fabs(z[j]) > l1
                            e1[j] = e2[j] = 1
                            violations += 1
                        end
                    end
                end
                if (violations == 0)
                    break
                end
            end

            # Scan for violations in rest
            violations = 0
            for j in eachindex(z)
                if e2[j] == 0
                    z[j] = crossprod(X, s, n, j) / n
                    l1 = lam[l] * m[j] * alpha
                    if fabs(z[j]) > l1
                        e1[j] = e2[j] = 1
                        violations += 1
                    end
                end
            end
            if violations == 0
                for i in eachindex(eta)
                    Eta[i, l] = eta[i]
                end
                break
            end
        end
    end

    return beta0, beta, Dev, Eta, iter

end
