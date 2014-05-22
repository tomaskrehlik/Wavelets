function modwtForward(V::Vector{Float64}, filter::waveletFilter, j::Int)
  (N,) = size(V)
  Wj = nans(N)
  Vj = nans(N)
  
  for t = 0:(N-1) 
    k = t
    Wjt = filter.h[1]*V[k+1]
    Vjt = filter.g[1]*V[k+1]
    for n = 1:(filter.L-1) 
      k = k - 2^(j-1)
      if (k < 0)
        k = k + ceil(-k/N)*N
      end
      Wjt = Wjt + filter.h[n+1]*V[k+1]
      Vjt = Vjt + filter.g[n+1]*V[k+1]
    end
    Wj[t+1] = Wjt
    Vj[t+1] = Vjt    
  end
  
  return (Wj, Vj)
end

function extendSeries(X::Array{Float64}, method::String, length::String, n::Int, j::Int)
  
  (N, ) = size(X)
    
  # determine final length 'n' of series after extension
  if (length == "arbitrary")
    if (n <= N)
     error("Invalid argument: 'n' must be greater than length of series when length='arbitrary'.")
    end
  elseif (length == "powerof2")
    k = N/(2^j)
    if (round(k) == k)
      error("Invalid argument: length of series should not be multiple of 2^j when length='powerof2'.")
    else
      n = ceil(k)*2^j
    end
  elseif (length == "double") 
    n = 2*N
  end

  if isa(X, Vector)
    # extend the series to length 'n'
    (method == "periodic") ? X = repeat(X, outer=[convert(Int,ceil(n/N))])[1:n] : ()
    (method == "reflection") ? X = repeat([X,X[end:-1:1]], outer=[convert(Int,ceil(n/N))])[1:n] : ()
    (method == "zeros") ? X = [X, zeros(N-n)] : ()
    (method == "mean") ? X = [X, mean(X).*ones(N-n)] : ()
    # if(method == "reflection.inverse") X = apply(X, 2, function(x,n,N){rep(c(x,2*x[N]-rev(x)),length=n)}, n=n, N=N)
  else
    # extend the series to length 'n'
    (method == "periodic") ? X = mapslices(x -> repeat(x, outer=[convert(Int,ceil(n/N))])[1:n], X, 1) : ()
    (method == "reflection") ? X = mapslices(x -> repeat(x, outer=[convert(Int,ceil(n/N))])[1:n], [X,X[end:-1:1]], 2) : ()
    (method == "zeros") ? X = mapslices(x -> [x, zeros(N-n)], X, 2) : ()
    (method == "mean") ? X = mapslices(x -> [x, mean(x).*ones(N-n)], X, 2) : ()
    # if(method == "reflection.inverse") X = apply(X, 2, function(x,n,N){rep(c(x,2*x[N]-rev(x)),length=n)}, n=n, N=N)
  end

  return X
end

function modwt(X::Matrix, filter::String, nLevels::Int, boundary::String)
  @assert filter=="haar"
  @assert boundary=="periodic"
  # require("waveletFilters.jl")

  # get wavelet coeficients and length
  # call some function
  filter = haarFilter(nLevels, true)

  # convert X to a matrix
  if isa(X, Vector)
    (N,) = size(X)
    nSeries = 1
  else
    (N, nSeries) = size(X)
  end
  # reflect X for reflection method
  if (boundary == "reflection")
    X = extendSeries(X)
    N = 2*N
  end

  # initialize variables for pyramid algorithm
  nBoundary = zeros(nLevels)
  WCoefs = zeros(N, nLevels, nSeries)
  VCoefs = zeros(N, nLevels, nSeries)
  

  # implement the pyramid algorithm
  for i=1:nSeries
    Vj = X[:,i]
    for j=1:nLevels
      (WCoefs[:,j,i], VCoefs[:,j,i]) = modwtForward(Vj,filter,j)
      Vj = convert(Vector, VCoefs[:,j,i])
      Lj = (2^j-1)*(filter.L-1)+1
      nBoundary[j] = min(Lj,N)
    end
  end
  
  return (WCoefs, VCoefs)
end