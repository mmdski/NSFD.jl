# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# StaggeredField
abstract type AbstractStaggeredScalarField end

function Base.getindex(field::AbstractStaggeredScalarField, i::Int, j::Int)
    return field.values[i + 1, j + 1]
end

function Base.getindex(field::AbstractStaggeredScalarField, I::UnitRange, J::UnitRange)
    return field.values[I .+ 1, J .+ 1]
end

function Base.setindex!(field::AbstractStaggeredScalarField, v::Real, i::Int, j::Int)
    return (field.values[i + 1, j + 1] = v)
end

function Base.size(field::AbstractStaggeredScalarField)
    return (field.grid.data.nx, field.grid.data.ny)
end

function Base.size(field::AbstractStaggeredScalarField, d::Int)
    return size(field)[d]
end

function Base.axes(field::AbstractStaggeredScalarField)
    return Base.OneTo.(size(field))
end

function Base.axes(field::AbstractStaggeredScalarField, d::Int)
    return Base.OneTo(size(field, d))
end

function Base.show(io::IO, field::AbstractStaggeredScalarField)
    nx, ny = field.grid.data.nx, field.grid.data.ny
    return print(io, "$(typeof(field))(size = ($nx, $ny))")
end

function set!(field::AbstractStaggeredScalarField, v::Real)
    return fill!(field.values, v)
end

function set!(field::AbstractStaggeredScalarField, f::Function)
    I, J = indices(field.grid)
    for i in I, j in J
        field[i, j] = f(field.grid[i, j]...)
    end
    return field
end

struct UField <: AbstractStaggeredScalarField
    grid::UGrid
    values::Matrix{Float64}

    function UField(grid_data::GridData)
        grid = UGrid(grid_data)
        values = Matrix{Float64}(undef, grid.data.nx + 2, grid.data.ny + 2)
        return new(grid, values)
    end
end

function UField(grid_data::GridData, init::Float64)
    u = UField(grid_data)
    set!(u, init)
    return u
end

struct VField <: AbstractStaggeredScalarField
    grid::VGrid
    values::Matrix{Float64}

    function VField(grid_data::GridData)
        grid = VGrid(grid_data)
        values = Matrix{Float64}(undef, grid.data.nx + 2, grid.data.ny + 2)
        return new(grid, values)
    end
end

function VField(grid_data::GridData, init::Float64)
    v = VField(grid_data)
    set!(v, init)
    return v
end

struct PField <: AbstractStaggeredScalarField
    grid::PGrid
    values::Matrix{Float64}

    function PField(grid_data::GridData)
        grid = PGrid(grid_data)
        values = Matrix{Float64}(undef, grid.data.nx + 2, grid.data.ny + 2)
        return new(grid, values)
    end
end

function PField(grid_data::GridData, init::Float64)
    p = PField(grid_data)
    set!(p, init)
    return p
end

struct StaggeredVectorField
    x::UField
    y::VField
    function StaggeredVectorField(grid_data::GridData)
        return new(UField(grid_data), VField(grid_data))
    end
end

function StaggeredVectorField(grid_data::GridData, u₀::Tuple{Float64,Float64})
    u = StaggeredVectorField(grid_data)
    set!(u.x, u₀[1])
    set!(u.y, u₀[2])
    return u
end

function Base.getindex(field::StaggeredVectorField, i::Int, j::Int)
    return StaggeredVector(field.x[i, j], field.y[i, j])
end

function Base.setindex!(field::StaggeredVectorField, vec::StaggeredVector, i::Int, j::Int)
    field.x[i, j] = vec.x
    field.y[i, j] = vec.y
    return field
end

function Base.size(field::StaggeredVectorField)
    return (field.x.grid.data.nx, field.x.grid.data.ny)
end

function Base.size(field::StaggeredVectorField, d::Int)
    return size(field)[d]
end

function Base.axes(field::StaggeredVectorField)
    return Base.OneTo.(size(field))
end

function Base.axes(field::StaggeredVectorField, d::Int)
    return Base.OneTo(size(field, d))
end

function Base.show(io::IO, field::StaggeredVectorField)
    nx, ny = field.x.grid.data.nx, field.x.grid.data.ny
    return print(io, "$(typeof(field))(size = ($nx, $ny))")
end

function set!(field::StaggeredVectorField, vec::StaggeredVector)
    set!(field.x, vec.x)
    set!(field.y, vec.y)
    return field
end

function set!(field::StaggeredVectorField, u::Tuple{Real,Real})
    set!(field.x, u[1])
    set!(field.y, u[2])
    return field
end
