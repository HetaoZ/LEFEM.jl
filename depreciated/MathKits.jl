@doc """
A collection of useful mathematical functions by Julia language. Author: https://gitee.com/iamzhtr

In REPL use `?MathKits.f` get help on the function `f`, for example, `?MathKits.pnpoly!`. 

Available functions: `between`, `betweeneq`, `crossproduct`, `distance_to_segment`, `norm2`, `pnpoly`/`pnpoly!`, `polygon_area`, `product`/`product!`, `root_on_segment`
""" ->
module MathKits
using Markdown, LinearAlgebra, Statistics

# -----------------------------------
# macros
mutable struct Point2d
    x::Float64
    y::Float64
end

# ----------------------------------- 


# functions
function affine_matrix(theta::Float64, d::Vector{Float64}; s::Vector{Float64} = [1.0, 1.0])
    A1 = rotate_matrix(theta) 
    A2 = d
    return A1, A2
end

function between(p::Array, point1::Array, point2::Array)
    ans = true
    for i = 1:length(p)
        if !(point1[i] < p[i] < point2[i])
            ans = false
            break
        end
    end
    return ans
end

function betweeneq(p::Array, point1::Array, point2::Array)
    ans = true
    for i = 1:length(p)
        if !(point1[i] <= p[i] <= point2[i])
            ans = false
            break
        end
    end
    return ans
end

function bi_lin_interp(x::T, y::T, φLD::T, xLD::T, yLD::T, φRD::T, xRD::T, yRD::T, φLU::T, xLU::T, yLU::T, φRU::T, xRU::T, yRU::T) where T <: Real
    @assert !isnan(φLD + φRD + φLU + φRU)
    a = (xRD - xLD) * (yLU - yLD)
    aLD = (x - xLD) * (y - yLD)
    aRD = (xRD - x) * (y - yRD)
    aLU = (x - xLU) * (yLU - y)
    aRU = (xRU - x) * (yRU - y)
    return (φLD * aRU + φRD * aLU + φLU * aRD + φRU * aLD) / a
end

@doc """
`ans = crossproduct(a::Vector{T} where T <: Number, b::Vector{T} where T <: Number)`

Return the cross product of two 3D or 2D vectors.
""" ->
function crossproduct(a::Vector{T} where T <: Number, b::Vector{T} where T <: Number)
    if (length(a), length(b)) == (3, 3)
        return eltype(a)[a[2]*b[3]-b[2]*a[3], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-b[1]*a[2]] 
    elseif (length(a), length(b)) == (2, 2)
        return eltype(a)(a[1]*b[2]-b[1]*a[2])
    else
        error("undef dim")
    end
end


function d2udx2(h::Float64, u1::Float64, u2::Float64, u3::Float64, u4::Float64, u5::Float64)
    return [-1 16 -30 16 -1] * [u1, u2, u3, u4, u5] / (12 * h^2)
end
function d2udx2(h::Float64, u::Vector{Float64})
    return [-1 16 -30 16 -1] * u / (12 * h^2)
end

@doc """
Second order cross difference on a (5x5)-elements stencil

d(du/dx)dy
""" ->
function d2udxdy(hx::Float64, hy::Float64, u::Array{Float64})
    return dudx(hy, [dudx(hx, u[:, j]) for j = 1:size(u, 2)])
end

@doc """
Second order cross difference on a (5x5)-elements stencil

d(du/dy)dx
""" ->
function d2udydx(hx::Float64, hy::Float64, u::Array{Float64})
    return dudx(hx, [dudx(hy, u[i, :]) for i = 1:size(u, 1)])
end


function dudx(h::Float64, u1::Float64, u2::Float64)
    return (u2 - u1)/h
end

@doc """
First order difference on a 4-elements stencil
""" ->
function dudx(h::Float64, u1::Float64, u2::Float64, u3::Float64, u4::Float64)
    return (-2 * u1 - 3 * u2 + 6 * u3 - u4) / (6 * h)
end

function dudx(h::Float64, u::Vector{Float64})
    return dudx(h, Tuple(u)...)
end

@doc """
`ans = distance_to_segment(vertex1::Array{T} where T <: Number, vertex2::Array{T} where T <: Number, point::Array{T} where T <: Number)`

Return the distance of `point` to a finite-length segment from `vertex1` to `vertex2`.
""" ->
function distance_to_segment(point::Array{T} where T <: Number, vertex1::Array{T} where T <: Number, vertex2::Array{T} where T <: Number)
    root, lambda = root_on_segment(point, vertex1, vertex2)
    if betweeneq(root, vertex1, vertex2)        
        return norm(point - root)
    else
        return min(norm(point - vertex1), norm(point - vertex2))
    end
end

function get_normal(v::Vector{Float64})
    return normalize(rotate_matrix(pi/2) * v)
end

function get_volume(xs::Array{Float64}...)
    dim = length(xs)
    if dim == 1
        volume = maximum(xs[1]) - minimum(xs[1])
    elseif dim == 2
        volume = polygon_area(xs[1], xs[2])
    else
        println("You need to calculate the volume by yourself!")
    end
    return volume
end

function lin_interp(x::T, φL::T, xL::T, φR::T, xR::T) where T <: Real
    @assert !isnan(φL + φR)
    a = xR - xL
    aL = x - xL
    aR = xR - x
    return (φL * aR + φR * aL) / a
end


@doc """
point: (x0, y0)

line: (x1, y1) -- (x2, y2)
""" ->
function mirrored_on_line(x0::Float64, y0::Float64, x1::Float64, y1::Float64, x2::Float64, y2::Float64)
    lambda = ( (x0- x1) * (x2 - x1) + (y0 - y1) * (y2 - y1) )  / ( (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) )
    xH = [x1, y1] + [x2 - x1, y2 - y1] * lambda
    xI = 2 * xH - [x0, y0]  
    return xI, xH, lambda
end

@doc """
`ans = norm2(v::Array{T} where T <: Number)`

Return the sum of all squares of element of the array `v`, i.e., `norm2!(v) = norm!(v)^2`.
""" ->
norm2(v::Array{T} where T <: Number) = norm(v)^2

@doc """
```julia
ans = pnpoly(vertices_x::Vector, vertices_y::Vector, point_x::Real, point_y::Real)
```
Find if a point lies within a 2D polygon. Return 0 if it lies inside, else return 1. This algorithm was proposed by W. Randolph Franklin: https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html

`vertices_x`/`vertices_y`: Vector of `x`/`y` of the polygon vertices. Note that `x[end] != x[1]` (very important).
 
`point_x`/`point_y`: `x`/`y` of the point.
""" ->  
function pnpoly(vertices_x::Vector, vertices_y::Vector, point_x::Real, point_y::Real)
    @assert length(vertices_x) == length(vertices_y)
    nvert = length(vertices_x)
    VertX = [vertices_x; vertices_x[1]]
    VertY = [vertices_y; vertices_y[1]]
    if ( (point_x < minimum(VertX)) || (point_x > maximum(VertX)) || (point_y < minimum(VertY)) || (point_y > maximum(VertY)) )
        c = 1
    else
        c = 1
        for i = 1:nvert
            j = i + 1
            if ( (VertY[i] > point_y) ⊻ (VertY[j] > point_y ) ) && ( point_x < (VertX[j] - VertX[i]) * (point_y - VertY[i]) / (VertY[j] - VertY[i]) + VertX[i] )
                # ⊻ : \xor
                c = 1 - c
            end
        end
    end
    # c = 0: inside, 1: outside
    return c
end

@doc """
```julia
ans = pnpoly!(vertices_x::Vector, vertices_y::Vector, point_x::Real, point_y::Real)
```
Find if a point lies within a 2D polygon. Return 0 if it lies inside, else return 1. This algorithm was proposed by W. Randolph Franklin: https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html

`vertices_x`/`vertices_y`: Vector of `x`/`y` of the polygon vertices. Note that `x[end] != x[1]` (very important).
 
`point_x`/`point_y`: `x`/`y` of the point.
""" ->  
function pnpoly!(vertices_x::Vector, vertices_y::Vector, point_x::Real, point_y::Real)
    @assert length(vertices_x) == length(vertices_y)
    nvert = length(vertices_x)
    VertX = [vertices_x; vertices_x[1]]
    VertY = [vertices_y; vertices_y[1]]
    if ( (point_x < minimum(VertX)) || (point_x > maximum(VertX)) || (point_y < minimum(VertY)) || (point_y > maximum(VertY)) )
        c = 1
    else
        c = 1
        for i = 1:nvert
            j = i + 1
            if ( (VertY[i] > point_y) ⊻ (VertY[j] > point_y ) ) && ( point_x < (VertX[j] - VertX[i]) * (point_y - VertY[i]) / (VertY[j] - VertY[i]) + VertX[i] )
                # ⊻ : \xor
                c = 1 - c
            end
        end
    end
    # c = 0: inside, 1: outside
    return c
end
     
#计算任意多边形的面积，顶点按照顺时针或者逆时针方向排列
function polygon_area(x::Array, y::Array)
    n = length(x)
    if n < 3 
        return 0.
    end
    s = x[n]*(y[1] - y[n-1]) + x[1]*(y[2] - y[n])
    for i = 2:n-1
        s += x[i]*(y[i+1] - y[i-1])
    end
    return abs(s) * 0.5
end

@doc """
`ans = product(v::Array{T} where T <: Number)`

Return the product of all elements of the array `v`.
""" ->
function product(v::Array{T} where T <: Number)
    s = 1.
    for a in v
        s *= a
    end
    return eltype(v)(s)
end

@doc """
`ans = product!(v::Array{T} where T <: Number)`

Return the product of all elements of the array `v`.
""" ->
function product!(v::Array{T} where T <: Number)
    s = 1.
    for a in v
        s *= a
    end
    return eltype(v)(s)
end

@doc """
`root, lambda = root_on_segment(point::Array{T} where T <: Number, vertex1::Array{T} where T <: Number, vertex2::Array{T} where T <: Number)`

Return the `root` point and length ratio `lambda` of `point` on the segment from `vertex1` to `vertex2`. The length ratio 'lambda = ||root - vertex1|| / ||vertex2 - vertex1||`. Note that `lambda < 0` or `lambda > 1` when `root` lies out of the segment. 
""" ->
function root_on_segment(point::Array{T} where T <: Number, vertex1::Array{T} where T <: Number, vertex2::Array{T} where T <: Number)
    AC, AB = point - vertex1, vertex2 - vertex1
    lambda = norm(AC .* AB) / norm(AB .* AB)
    root = vertex1 + lambda * AB
    return root, lambda
end

@doc """
arc angle, not degree
""" ->
function rotate_matrix(angle::Real)
    return [cos(angle) -sin(angle); sin(angle) cos(angle)]
end

function shapestar(n::Int; center::Vector = [0, 0], R::Real = 1, r::Real = 1, angle::Real = 0)
    pi = π
    X, Y = Array{Float64}(undef,2*n), Array{Float64}(undef,2*n)
    for i=1:n
        X[2*i-1] = center[1] + R*cos(pi/180*((i-1)*(-360/n)+angle+90))
        Y[2*i-1] = center[2] + R*sin(pi/180*((i-1)*(-360/n)+angle+90))
        X[2*i] = center[1] + r*cos(pi/180*((i-1+0.5)*(-360/n)+angle+90))
        Y[2*i] = center[2] + r*sin(pi/180*((i-1+0.5)*(-360/n)+angle+90))
    end

    return X, Y
end

# ----------------------------------- #
end