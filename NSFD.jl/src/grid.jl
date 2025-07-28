# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Axis
struct Axis{Offset}
    Δ::Float64
end

function Base.getindex(ax::Axis{O}, i::Int) where {O}
    return ax.Δ * (i - 1 + O)
end

const CAxis = Axis{0.5}
const RAxis = Axis{1.0}

function offset(::Axis{O}) where {O}
    return O
end

# GridAxes
struct GridAxes{XAxis,YAxis}
    x::XAxis
    y::YAxis
end

function GridAxes(::Type{XAxis}, ::Type{YAxis},
                  data::GridData) where {XAxis<:Axis,YAxis<:Axis}
    x_axis = XAxis(data.Δx)
    y_axis = YAxis(data.Δy)
    return GridAxes{XAxis,YAxis}(x_axis, y_axis)
end

function Base.getindex(axes::GridAxes, i::Int, j::Int)
    return axes.x[i], axes.y[j]
end

# AbstractGrid
abstract type AbstractGrid end

function Base.getindex(grid::AbstractGrid, i::Int, j::Int)
    return grid.axes[i, j]
end

function Base.show(io::IO, grid::AbstractGrid)
    nx, ny = grid.data.nx, grid.data.ny
    x_off = offset(grid.axes.x)
    y_off = offset(grid.axes.y)
    Δx = grid.data.Δx
    Δy = grid.data.Δy
    print(io,
          "$(typeof(grid))(size = ($nx, $ny), Δ = ($Δx, $Δy), offset = ($x_off, $y_off))")
    return nothing
end

function indices(grid::AbstractGrid)
    return 0:(grid.data.nx + 1), 0:(grid.data.ny + 1)
end

function interior_indices(grid::AbstractGrid)
    return 1:(grid.data.nx), 1:(grid.data.ny)
end

# u-Grid
# velocity component in x-direction
struct UGrid <: AbstractGrid
    data::GridData
    axes::GridAxes{RAxis,CAxis}
end

function UGrid(data::GridData)
    return UGrid(data, GridAxes(RAxis, CAxis, data))
end

# v-Grid
# velocity component in y-direction
struct VGrid <: AbstractGrid
    data::GridData
    axes::GridAxes{CAxis,RAxis}
end

function VGrid(data::GridData)
    return VGrid(data, GridAxes(CAxis, RAxis, data))
end

# p-Grid
# pressure grid
struct PGrid <: AbstractGrid
    data::GridData
    axes::GridAxes{CAxis,CAxis}
end

function PGrid(data::GridData)
    return PGrid(data, GridAxes(CAxis, CAxis, data))
end

struct StaggeredGrid
    u::UGrid
    v::VGrid
    p::PGrid
end

function StaggeredGrid(grid_data::GridData)
    u = UGrid(grid_data)
    v = VGrid(grid_data)
    p = PGrid(grid_data)
    return StaggeredGrid(u, v, p)
end

function meshgrid(grid::AbstractGrid)
    X = Matrix{Float64}(undef, grid.data.nx, grid.data.ny)
    Y = similar(X)
    for i in 1:(grid.data.nx)
        for j in 1:(grid.data.ny)
            X[i, j], Y[i, j] = grid[i, j]
        end
    end
    return X, Y
end
