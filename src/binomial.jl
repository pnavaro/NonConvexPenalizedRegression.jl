function p_binomial(eta)
  if eta > 10
    return 1
  elseif eta < -10
    return 0
  else 
    return exp(eta)/(1+exp(eta))
  end
end

"""
Weighted cross product of y with jth column of x
"""
function wcrossprod(X, y, w, j)
  nn = n*j
  val=0
  for i in eachindex(y)
      val += X[nn+i]*y[i]*w[i]
  end

  return val

end

fmax2( x, y) = (x < y) ? y : x

"""
Coordinate descent for binomial models
"""
function cdfit_binomial( X, y, penalty, lambda, eps, max_iter, gamma, multiplier, alpha, dfmax, user, warn)

    n = length(y)
    p = length(X)/n
    L = length(lambda)
    beta0 = zeros(L)
    beta = zeros(L*p)
    Dev = zeros(L)
    Eta = zeros(L*n)
    iter = zeros(Int, L)
    a = zeros(p)    # Beta from previous iteration
    a0 = 0          # Beta0 from previous iteration
    lam = lambda
    tot_iter = 0;
    m = multiplier
    r = zeros(n)
    w = zeros(n)
    s = zeros(n)
    z = zeros(p)
    eta = zeros(n)
    e1 = zeros(Int, p)
    e2 = zeros(Int, p);

    # Initialization
    ybar = sum(y) / n
    a0 = b0[0] = log(ybar/(1-ybar))
    nullDev = 0
    for i in eachindex(y) 
        nullDev = nullDev - y[i]*log(ybar) - (1-y[i])*log(1-ybar)
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
    for l in lstart:L-1

        if l != 0
            # Assign a, a0
            a0 = b0[l-1]
            for j in eachindex(a)
                a[j] = b[(l-1)*p+j]
            end

            # Check dfmax
            nv = 0;
            for j in eachindex(a)
              if a[j] != 0 
                nv +=1
              end
            end

            if (nv > dfmax) || (tot_iter == max_iter)
              for ll in l:L-1
                  iter[ll] = missing
              end
              break
            end

            # Determine eligible set
            penalty == "lasso" && (cutoff = 2*lam[l] - lam[l-1])
            penalty == "MCP" && (cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lam[l-1]))
            penalty == "SCAD" && (cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lam[l-1]))

            for j in eachindex(z)
                if abs(z[j]) > (cutoff * alpha * m[j])
                   e2[j] = 1;
                end
            end

        else

            # Determine eligible set
            lmax = 0.0
            for j in eachindex(z) 
                if abs(z[j]) > lmax
                   lmax = abs(z[j])
                end
            end
            penalty == "lasso" && ( cutoff = 2*lam[l] - lmax)
            penalty == "MCP" && ( cutoff = lam[l] + gamma/(gamma-1)*(lam[l] - lmax))
            penalty == "SCAD" && ( cutoff = lam[l] + gamma/(gamma-2)*(lam[l] - lmax))
            for j in eachindex(z)
                if abs(z[j] > (cutoff * alpha * m[j])) 
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
                    for i in eachindex(eta)
                        pi = p_binomial(eta[i])
                        w[i] = fmax2(pi*(1-pi), 0.0001)
                        s[i] = y[i] - pi
                        r[i] = s[i]/w[i]
                        if (y[i]==1) REAL(Dev)[l] = REAL(Dev)[l] - log(pi)
                        if (y[i]==0) REAL(Dev)[l] = REAL(Dev)[l] - log(1-pi)
                    end
                    if REAL(Dev)[l]/nullDev < .01
                        if warn 
                           @warning("Model saturated; exiting...")
                        end
                        for ll =l:L-1
                            iter[ll] = missing
                        end
                        tot_iter = max_iter
                        break
                    end

                    # Intercept
                    xwr = crossprod(w, r, n, 0);
                    xwx = sum(w, n);
                    b0[l] = xwr/xwx + a0;
                    for i in eachindex(r)
                      si = b0[l] - a0;
                      r[i] -= si;
                      eta[i] += si;
                    end

                    maxChange = abs(si)*xwx/n
                    
                    # Covariates
                    for j in eachindex(a)
                      if e1[j]

                        # Calculate u, v
                        xwr = wcrossprod(X, r, w, n, j)
                        xwx = wsqsum(X, w, n, j)
                        u = xwr/n + (xwx/n)*a[j]
                        v = xwx/n;

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
                          for i in eachindex(r)
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
                  if fabs(z[j]) > l1
                    e1[j] = e2[j] = 1;
                    violations +=1
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
              if e2[j]==0
                z[j] = crossprod(X, s, j)/n
                l1 = lam[l] * m[j] * alpha
                if abs(z[j]) > l1
                  e1[j] = e2[j] = 1
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
