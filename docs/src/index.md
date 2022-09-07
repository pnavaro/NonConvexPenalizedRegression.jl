# NonConvexPenalizedRegression.jl

Regularization Paths for SCAD and MCP Penalized Regression Models

This is partial translationn Julia  of the R package [ncvreg](http://pbreheny.github.io/ncvreg/).  
Only `gaussian` family is translated.

Algorithm is described in **Breheny P and Huang J (2011)** "Coordinate descent algorithms for nonconvex penalized regression, with applications to biological feature selection". *Annals of Applied Statistics*, 5: 232–253

I needed to do regression with SCAD penalty and I can't find it in any Julia package. 
Perharps it is now implemented in [MLJLinearModels.jl](https://github.com/alan-turing-institute/MLJLinearModels.jl).


```julia-repl
julia> using LinearAlgebra 
julia> using Random
julia> using NonConvexPenalizedRegression
julia> using RCall
julia> rng = MersenneTwister(1234);
julia> n, p = 50, 5
(50, 5)

julia> X = randn(rng, n, p)              # feature matrix
50×5 Array{Float64,2}:
  0.867347   -1.22672     0.183976    0.0377058   1.48494
 -0.901744   -0.541716   -1.27635    -0.490009    1.23969
  ⋮
 -1.00978    -1.66323    -0.744522    0.427383   -1.37986
 -0.543805   -0.521229   -0.191176   -0.492253   -0.984217

julia> a0 = collect(1:p)                # ground truths
5-element Array{Int64,1}:
 1
 2
 3
 4
 5

julia> y = X * a0 + 0.1 * randn(n) # generate response
50-element Array{Float64,1}:
   6.411769869798991
  -1.59817739694925
   ⋮
 -11.769008550241434
  -9.107294931708777

julia> XX = hcat(X, randn(rng, n, p))
50×10 Array{Float64,2}:
  0.867347   -1.22672     0.183976    0.0377058  …   1.69129     0.969694    0.222167    -0.953909
  ⋮                                              ⋱
 -0.543805   -0.521229   -0.191176   -0.492253      -2.5788     -0.329958    0.00775707  -0.370354

julia> @rput XX
50×10 Array{Float64,2}:
  0.867347   -1.22672     0.183976    0.0377058  …   1.69129     0.969694    0.222167    -0.953909
  ⋮                                              ⋱
 -0.543805   -0.521229   -0.191176   -0.492253      -2.5788     -0.329958    0.00775707  -0.370354

julia> @rput y
50-element Array{Float64,1}:
   6.411769869798991
  -1.59817739694925
   ⋮
 -11.769008550241434
  -9.107294931708777

julia> R"library(ncvreg)"
RObject{StrSxp}
[1] "ncvreg"    "stats"     "graphics"  "grDevices" "utils"     "datasets"
[7] "methods"   "base"

julia> R"scad <- coef(ncvreg(XX, y, lambda=0.2, penalty='SCAD', eps=.0001))"
RObject{RealSxp}
 (Intercept)           V1           V2           V3           V4           V5
-0.003322903  1.025666431  2.000108987  2.983498545  3.997543804  4.982264144
          V6           V7           V8           V9          V10
 0.000000000  0.000000000  0.000000000  0.000000000  0.000000000

julia> @rget scad
11-element Array{Float64,1}:
 -0.0033229032078964105
  1.0256664305062173
  2.000108987156345
  2.9834985454557157
  3.997543803872521
  4.982264143916537
  0.0
  0.0
  0.0
  0.0
  0.0

julia> println( " R scad = $scad")
 R scad = [-0.0033229032078964105, 1.0256664305062173, 2.000108987156345, 2.9834985454557157, 3.997543803872521, 4.982264143916537, 0.0, 0.0, 0.0, 0.0, 0.0]

julia> λ = [0.2]
1-element Array{Float64,1}:
 0.2

julia> scad = NonConvexPenalizedRegression.coef(SCAD(XX, y, λ))
SCAD([-0.003322960709765954; 1.0256660512338405; … ; 0.0; 0.0])

julia> println( " Julia scad = $scad")
 Julia scad = SCAD([-0.003322960709765954; 1.0256660512338405; 2.00010933635426; 2.983498839847109; 3.99754375703709; 4.982264100245242; 0.0; 0.0; 0.0; 0.0; 0.0])

```

```@autodocs
Modules = [NonConvexPenalizedRegression]
```

```@index
```

