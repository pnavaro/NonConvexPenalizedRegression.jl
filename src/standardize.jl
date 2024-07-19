function standardize(X)

    n, p = size(X)

    XX = similar(X)
    c = zeros(Float64, p)
    s = zeros(Float64, p)

    for j = 1:p

        c[j] = 0
        for i = 1:n
            c[j] += X[i, j]
        end
        c[j] = c[j] / n
        for i = 1:n
            XX[i, j] = X[i, j] - c[j]
        end

        s[j] = 0
        for i = 1:n
            s[j] += XX[i, j]^2
        end
        s[j] = sqrt(s[j] / n)
        for i = 1:n
            XX[i, j] = XX[i, j] / s[j]
        end

    end

    any(s .< 1e-6) && throw("Sigular matrix")

    ns = findall(s .> 1e-6)

    copy(XX[:, ns]), c, s

end
