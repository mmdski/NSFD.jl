# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

abstract type AbstractField end

Base.getindex(f::AbstractField, i::Int, j::Int) = f.data[i + 1, j + 1]
Base.setindex!(f::AbstractField, v::Float64, i::Int, j::Int) = (f.data[i + 1, j + 1] = v)

function interior_shape(f::AbstractField)
    return f.grid.x.n, f.grid.y.n
end

function grid_data(f::AbstractField)
    return GridData(f.grid.x.length, f.grid.y.length, f.grid.x.n, f.grid.y.n)
end

function meshgrid(f::AbstractField)
    nx = f.grid.x.n
    ny = f.grid.y.n
    X = Matrix{Float64}(undef, nx, ny)
    Y = similar(X)
    for i in 1:nx
        for j in 1:ny
            X[i, j], Y[i, j] = f.grid[i, j]
        end
    end
    return X, Y
end

struct UField <: AbstractField
    grid::Grid
    data::Matrix{Float64}
end

function UField(grid_data::GridData)
    data = zeros(Float64, grid_data.nx + 2, grid_data.ny + 2)
    x_axis = Axis(grid_data.x_length, grid_data.nx, 0.0)
    y_axis = Axis(grid_data.y_length, grid_data.ny, 0.5)
    grid = Grid(x_axis, y_axis)
    return UField(grid, data)
end

struct VField <: AbstractField
    grid::Grid
    data::Matrix{Float64}
end

function VField(grid_data::GridData)
    data = zeros(Float64, grid_data.nx + 2, grid_data.ny + 2)
    x_axis = Axis(grid_data.x_length, grid_data.nx, 0.5)
    y_axis = Axis(grid_data.y_length, grid_data.ny, 0.0)
    grid = Grid(x_axis, y_axis)
    return VField(grid, data)
end

struct VelocityField
    u::UField
    v::VField
end

struct PField <: AbstractField
    grid::Grid
    data::Matrix{Float64}
end

function PField(grid_data::GridData)
    data = zeros(Float64, grid_data.nx + 2, grid_data.ny + 2)
    x_axis = Axis(grid_data.x_length, grid_data.nx, 0.5)
    y_axis = Axis(grid_data.y_length, grid_data.ny, 0.5)
    grid = Grid(x_axis, y_axis)
    return PField(grid, data)
end
