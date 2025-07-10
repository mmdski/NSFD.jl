# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct Axis
    length::Float64
    n::Int
    Δ::Float64
    frac::Float64
end

function Axis(length::Float64, n::Int, frac::Float64)
    return Axis(length, n, length / n, frac)
end

Base.getindex(ax::Axis, i::Int) = ax.Δ * (i - ax.frac)

struct GridData
    x_length::Float64
    y_length::Float64
    nx::Int64
    ny::Int64
end

struct Grid
    x::Axis
    y::Axis
end

Base.getindex(grid::Grid, i::Int, j::Int) = grid.x[i], grid.y[j]
