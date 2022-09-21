"""
Coordinate descent for poisson models
"""
function cdfit_poisson( X, y, penalty, lambda, eps, max_iter, gamma, 
        multiplier, alpha, dfmax, user, warn)

    # Declarations
    n = length(y)
    p = length(X)/n
    L = length(lambda)
    beta0 = zeros(L)
    b0 = beta0
    beta = zeros(L*p)
    b = beta
    Dev = zeros(L)
    Eta = zeros(L*n)
    iter = zeros(Int, L)
    a = zeros(p)    # Beta from previous iteration
    a0 = 0                    # Beta0 from previous iteration
    lam = lambda
    tot_iter = 0
    m = multiplier
    r = zeros(n)
    w = zeros(n)
    s = zeros(n)
    z = zeros(p)
    eta = zeros(n)
    e1 = zeros(Int, p)
    e2 = zeros(Int, p)

    # Initialization
    ybar = sum(y)/n
    a0 = b0[0] = log(ybar)
    nullDev = 0
    for i in eachindex(y)
        if y[i] != 0
            nullDev += y[i]*log(y[i]/ybar) 
        end 
    end
    for i in eachindex(s)
        s[i] = y[i] - ybar 
    end
    for i in eachindex(eta)
        eta[i] = a0 
    end
    for j in eachindex(z)
        z[j] = crossprod(X, s, j)/n 
    end

    # If lam[0]=lam_max, skip lam[0] -- closed form sol'n available
    if user
        lstart = 0
    else
        lstart = 1
        Dev[0] = nullDev
        for i in eachindex(eta)
            Eta[i] = eta[i] 
        end
    end

    # Path
    for l=lstart:L-1
        if l != 0
            # Assign a, a0
            a0 = b0[l-1]
            for j in eachindex(a)
                a[j] = b[(l-1)*p+j] 
            end

            # Check dfmax
            nv = 0
            for j in eachindex(a)
                if a[j] != 0 
                    nv += 1 
                end
            end
            if (nv > dfmax) | (tot_iter == max_iter)
              for ll = l:L-1
                  iter[ll] = missing 
              end
              break
            end

            # Determine eligible set
            if penalty == "lasso" 
                cutoff = 2*lam[l] - lam[l-1] 
            end

            if penalty == "MCP"
                cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lam[l-1]) 
            end

            if penalty == "SCAD" 
                cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lam[l-1]) 
            end

            for j in eachindex(z)
                if abs(z[j]) > (cutoff * alpha * m[j])
                    e2[j] = 1 
                end 
            end

        else 

            # Determine eligible set
            lmax = 0
            for j in eachindex(z)
                if abs(z[j]) > lmax
                    lmax = abs(z[j])
                end 
            end
            if penalty == "lasso"
                cutoff = 2*lam[l] - lmax 
            end
            if penalty == "MCP"
                cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lmax) 
            end
            if penalty == "SCAD"
                cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lmax) 
            end
            for j in eachindex(z)
                if abs(z[j]) > (cutoff * alpha * m[j])
                    e2[j] = 1 
                end 
            end
        end

        while tot_iter < max_iter 
            while tot_iter < max_iter
                while tot_iter < max_iter 

                    iter[l] += 1
                    tot_iter += 1
                    Dev[l] = 0
                    maxChange = 0
                    for i in eachindex(eta)
                        mu = exp(eta[i])
                        w[i] = mu
                        s[i] = y[i] - mu
                        r[i] = s[i]/w[i]
                        if y[i]!=0
                            Dev[l] += y[i]*log(y[i]/mu)
                        end
                    end
                    if Dev[l]/nullDev < .01
                        if warn
                            @warn "Model saturated exiting..."
                        end
                        for ll=l:L-1
                            iter[ll] = missing
                        end
                        tot_iter = max_iter
                        break
                    end

                    # Intercept
                    xwr = crossprod(w, r, 0)
                    xwx = sum(w)
                    b0[l] = xwr/xwx + a0
                    for i in eachindex(r)
                        si = b0[l] - a0
                        r[i] -= si
                        eta[i] += si
                    end

                    # Covariates
                    for j in eachindex(e1)

                        if e1[j]

                            # Calculate u, v
                            xwr = wcrossprod(X, r, w, j)
                            xwx = wsqsum(X, w, j)
                            u = xwr/n + (xwx/n)*a[j]
                            v = xwx/n

                            # Update b_j
                            l1 = lam[l] * m[j] * alpha
                            l2 = lam[l] * m[j] * (1-alpha)
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
                            if shift !=0
                                for i in eachindex(eta)
                                    si = shift*X[j*n+i]
                                    r[i] -= si
                                    eta[i] += si
                                end
                                if abs(shift)*sqrt(v) > maxChange
                                    maxChange = abs(shift)*sqrt(v)
                                end
                            end
                        end
                    end

                    # Check for convergence
                    a0 = b0[l]
                    for j in eachindex(a)
                        a[j] = b[l*p+j]
                    end
                    if maxChange < eps 
                        break
                    end
                end

                # Scan for violations in strong set
                violations = 0
                for j in eachindex(z)
                    if e1[j]==0 && e2[j]==1
                        z[j] = crossprod(X, s, j)/n
                        l1 = lam[l] * m[j] * alpha
                        if abs(z[j]) > l1
                            e1[j] = 1
                            e2[j] = 1
                            violations += 1
                        end
                    end
                end
                if violations==0
                    break 
                end
            end

            # Scan for violations in rest
            violations = 0
            for j in eachindex(z)
                if e2[j] == 0
                    z[j] = crossprod(X, s, j)/n
                    l1 = lam[l] * m[j] * alpha
                    if abs(z[j]) > l1
                        e1[j] = 1
                        e2[j] = 1
                        violations += 1
                    end
                end
            end
            if violations==0
                for i in eachindex(eta)
                    Eta[n*l+i] = eta[i]
                end
                break
            end
        end
    end
    return res
end
