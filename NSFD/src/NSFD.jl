# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

module NSFD

export GridData, UGrid, VGrid, PGrid, meshgrid
export Field

struct Axis{Offset}
    Δ::Float64
end

Base.getindex(ax::Axis{O}, i::Int) where {O} = ax.Δ * (i - 1 + O)

CAxis = Axis{0.5}
RAxis = Axis{1.0}

struct GridData
    x_length::Float64
    y_length::Float64
    nx::Int64
    ny::Int64
    Δx::Float64
    Δy::Float64
    function GridData(x_length::Float64, y_length::Float64, nx::Int64, ny::Int64)
        Δx = x_length / nx
        Δy = y_length / ny
        return new(x_length, y_length, nx, ny, Δx, Δy)
    end
end

struct Grid{XAxis,YAxis}
    data::GridData
    x::XAxis
    y::YAxis
end

function Grid(::Type{XAxis}, ::Type{YAxis}, data::GridData) where {XAxis<:Axis,YAxis<:Axis}
    Δx = data.x_length / data.nx
    Δy = data.y_length / data.ny
    x_axis = XAxis(Δx)
    y_axis = YAxis(Δy)
    return Grid{XAxis,YAxis}(data, x_axis, y_axis)
end

Base.getindex(grid::Grid, i::Int, j::Int) = grid.x[i], grid.y[j]

UGrid(data) = Grid(RAxis, CAxis, data)
VGrid(data) = Grid(CAxis, RAxis, data)
PGrid(data) = Grid(CAxis, CAxis, data)

function meshgrid(grid::Grid)
    X = Matrix{Float64}(undef, grid.data.nx, grid.data.ny)
    Y = similar(X)
    for i in 1:(grid.data.nx)
        for j in 1:(grid.data.ny)
            X[i, j], Y[i, j] = grid[i, j]
        end
    end
    return X, Y
end

struct Field
    values::Matrix{Float64}
end

function Field(grid_data::GridData)
    values = zeros(Float64, grid_data.nx + 2, grid_data.ny + 2)
    return Field(values)
end

Base.getindex(f::Field, i::Int, j::Int) = f.data[i + 1, j + 1]
Base.setindex!(f::Field, v::Float64, i::Int, j::Int) = (f.data[i + 1, j + 1] = v)

end
