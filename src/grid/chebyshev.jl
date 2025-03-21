module BaryChebTools

using ..StaticArrays

export BaryCheb1D, interp1D, interpND, integrate1D, integrateND
######################################################
#---------------- 1D barycheb ------------------------
######################################################

"""
    function barychebinit(n::Integer)

Get Chebyshev nodes of first kind and corresponding barycentric Lagrange interpolation weights. 

Reference: Berrut, J.P. and Trefethen, L.N., 2004. Barycentric lagrange interpolation. SIAM review, 46(3), pp.501-517.
# Arguments
- `n`: order of the Chebyshev interpolation
# Returns
- Chebyshev nodes
- Barycentric Lagrange interpolation weights
"""
function barychebinit(n)
    x = zeros(Float64, n)
    w = similar(x)
    for i in 1:n
        c = (2i - 1)π / (2n)
        x[n - i + 1] = cos(c)
        w[n - i + 1] = (-1)^(i - 1) * sin(c)
    end
    return x, w
end

"""
    function barychebinit(::Type{T}, n) where {T}

Same as [`barychebinit(Tuple{Any})`](@ref), but returns two arrays of type `T`.
"""
function barychebinit(::Type{T}, n) where {T}
    x = zeros(T, n)
    w = similar(x)
    for i in 1:n
        c = (2i - 1)π / (2n)
        x[n - i + 1] = cos(c)
        w[n - i + 1] = (-1)^(i - 1) * sin(c)
    end
    return x, w
end

function vandermonde(x)
    # generate vandermonde matrix for x
    n = length(x)
    vmat = zeros(Float64, (n, n))
    for i in 1:n
        for j in 1:n
            vmat[i,j] = x[i]^(j-1)
        end
    end
    return vmat
end

function vandermonde(::Type{T}, x) where {T}
    # generate vandermonde matrix for x
    n = length(x)
    vmat = zeros(T, (n, n))
    for i in 1:n
        for j in 1:n
            vmat[i,j] = x[i]^(j-1)
        end
    end
    return vmat
end

function invvandermonde(order)
    # generate inverse of vandermonde matrix for cheb x of given order
    n = order
    vmat = zeros(Float64, (n, n))
    for i in 1:n
        for j in 1:n
            c = (2n-2i+1)π/(2n)
            x = cos(c)
            vmat[i,j] = x^(j-1)
        end
    end
    return inv(transpose(vmat))
end

function invvandermonde(::Type{T}, order) where {T}
    # generate inverse of vandermonde matrix for cheb x of given order
    n = order
    vmat = zeros(T, (n, n))
    for i in 1:n
        for j in 1:n
            c = (2n-2i+1)π/(2n)
            x = cos(c)
            vmat[i,j] = x^(j-1)
        end
    end
    return inv(transpose(vmat))
end


@inline function divfactorial(j, i)
    if i == 0
        return 1
    elseif i == 1
        return 1/j
    elseif i == -1
        return j-1
    else
        return factorial(j-1) / factorial(j+i-1)
    end
end


function weightcoef(a, i::Int, n)
    # integrate when i=1; differentiate when i=-1
    b = zeros(Float64, n)
    for j in 1:n
        if j+i-1 > 0
            # b[j] = a^(j+i-1)/(j+i-1)
            b[j] = a^(j+i-1) * divfactorial(j, i)
        elseif j+i-1 == 0
            b[j] = 1
        else
            b[j] = 0
        end
    end
    return b
end


function weightcoef(::Type{T}, a, i::Int, n) where {T}
    # integrate when i=1; differentiate when i=-1
    b = zeros(T, n)
    for j in 1:n
        if j+i-1 > 0
            # b[j] = a^(j+i-1)/(j+i-1)
            b[j] = a^(j+i-1) * divfactorial(j, i)
        elseif j+i-1 == 0
            b[j] = 1
        else
            b[j] = 0
        end
    end
    return b
end



@inline function calcweight(invmat, b)
    return invmat*b
end

"""
    function barycheb(n, x, f, wc, xc)

Barycentric Lagrange interpolation at Chebyshev nodes.

Reference: Berrut, J.P. and Trefethen, L.N., 2004. Barycentric lagrange interpolation. SIAM review, 46(3), pp.501-517.

# Arguments
- `n`: order of the Chebyshev interpolation
- `x`: coordinate to interpolate
- `f`: array of size n, function at the Chebyshev nodes
- `wc`: array of size n, Barycentric Lagrange interpolation weights
- `xc`: array of size n, coordinates of Chebyshev nodes
# Returns
- Interpolation result
"""
function barycheb(n, x, f, wc, xc)
    #print("type of data in cheb $(typeof.([n,x,f,wc,xc]))\n")
    #print("interp: $(f[1])\n")
    for j in 1:n
        if x == xc[j]
            return f[j]
        end    
    end

    num, den = 0.0, 0.0
    for j in 1:n
        q = wc[j] / (x - xc[j])
        num += q * f[j]
        den += q
    end
    return num / den
end

function barychebND(n, xs, f, wc, xc, DIM)
    haseq = false
    eqinds = zeros(Int, DIM)
    for i in 1:DIM
        for j in 1:n
            if xs[i] == xc[j]
                eqinds[i] = j
                haseq = true
            end
        end
    end

    if haseq
        newxs = [xs[i] for i in 1:DIM if eqinds[i]==0]
        newDIM = length(newxs)
        if newDIM == 0
            return f[CartesianIndex(eqinds...)]
        else
            newf = view(f, [(i==0) ? (1:n) : (i) for i in eqinds]...)
            return _barychebND_noneq(n, newxs, newf, wc, xc, newDIM)
        end
    else
        return _barychebND_noneq(n, xs, f, wc, xc, DIM)
    end
end

function _barychebND_noneq(n, xs, f, wc, xc, DIM)
    # deal with the case when there's no xs[i] = xc[j]
    inds = CartesianIndices(NTuple{DIM, Int}(ones(Int, DIM) .* n))
    num, den = 0.0, 0.0
    for (indi, ind) in enumerate(inds)
        q = 1.0
        for i in 1:DIM
            q *= wc[ind[i]] / (xs[i] - xc[ind[i]])
        end
        num += q * f[indi]
        den += q
    end
    return num / den
end

function barycheb2(n, x, f, wc, xc)
    for j in 1:n
        if x == xc[j]
            return f[j, :]
        end    
    end

    den = 0.0
    num = zeros(eltype(f), size(f)[2])
    for j in 1:n
        q = wc[j] / (x - xc[j])
        num += q .* f[j, :]
        den += q
    end
    return num ./ den
end

function chebint(n, a, b, f, invmat)
    wc = weightcoef(b, 1, n) .- weightcoef(a, 1, n)
    intw = calcweight(invmat, wc)
    return sum(intw .* f)
end

function chebdiff(n, x, f, invmat)
    wc = weightcoef(x, -1, n)
    intw = calcweight(invmat, wc)
    return sum(intw .* f)
end

struct BaryCheb1D{N,T}
    # wrapped barycheb 1d grid, x in [-1, 1]
    x::SVector{N, T}
    w::SVector{N, T}
    invmat::SMatrix{N, N, T}
end

function BaryCheb1D(N::Int,T) 
    x, w = barychebinit(T,N)
    invmat = invvandermonde(T,N)

    return BaryCheb1D{N,T}(x, w, invmat)
end
function  BaryCheb1D(N::Int) 
    x, w = barychebinit(Float64,N)
    invmat = invvandermonde(Float64,N)

    return BaryCheb1D{N,Float64}(x, w, invmat)
end 

BaryCheb1D{N}(x, w, invmat) where{N} = BaryCheb1D{N,Float64}(x, w, invmat)


Base.getindex(bc::BaryCheb1D, i) = bc.x[i]

function interp1D(data, xgrid::BaryCheb1D{N,T}, x) where {N,T}
    return barycheb(N, x, data, xgrid.w, xgrid.x)
end

function integrate1D(data, xgrid::BaryCheb1D{N,T}; x1=-1, x2=1) where {N,T}
    return chebint(N, x1, x2, data, xgrid.invmat)
end

function interpND(data, xgrid::BaryCheb1D{N,T}, xs) where {N,T}
    return barychebND(N, xs, data, xgrid.w, xgrid.x, length(xs))
end

function integrateND(data, xgrid::BaryCheb1D{N,T}, x1s, x2s) where {N,T}
    DIM = length(x1s)
    @assert DIM == length(x2s)

    intws = zeros(T, (DIM, N))
    for i in 1:DIM
        wc = weightcoef(x2s[i], 1, N) .- weightcoef(x1s[i], 1, N)
        intws[i, :] = calcweight(xgrid.invmat, wc)
    end

    result = 0.0
    inds = CartesianIndices(NTuple{DIM, Int}(ones(Int, DIM) .* N))
    for (indi, ind) in enumerate(inds)
        w = 1.0
        for i in 1:DIM
            w *= intws[i, ind[i]]
        end
        result += data[indi] * w
    end

    return result
end

function integrateND(data, xgrid::BaryCheb1D{N,T}, DIM) where {N,T}
    @assert N ^ DIM == length(data)

    intws = zeros(T, (DIM, N))
    wc = weightcoef(1.0, 1, N) .- weightcoef(-1.0, 1, N)
    intw = calcweight(xgrid.invmat, wc)

    result = 0.0
    inds = CartesianIndices(NTuple{DIM, Int}(ones(Int, DIM) .* N))
    for (indi, ind) in enumerate(inds)
        w = 1.0
        for i in 1:DIM
            w *= intw[ind[i]]
        end
        result += data[indi] * w
    end

    return result
end

end
