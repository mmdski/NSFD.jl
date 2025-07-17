# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# StaggeredField
abstract type AbstractStaggeredField end

function Base.getindex(field::AbstractStaggeredField, i::Int, j::Int)
    return field.values[i + 1, j + 1]
end

function Base.getindex(field::AbstractStaggeredField, I::UnitRange, J::UnitRange)
    return field.values[I .+ 1, J .+ 1]
end

function Base.setindex!(field::AbstractStaggeredField, v::Real, i::Int, j::Int)
    return (field.values[i + 1, j + 1] = v)
end

function Base.size(field::AbstractStaggeredField)
    return (field.grid.data.nx, field.grid.data.ny)
end

function Base.size(field::AbstractStaggeredField, d::Int)
    return size(field)[d]
end

function Base.axes(field::AbstractStaggeredField)
    return Base.OneTo.(size(field))
end

function Base.axes(field::AbstractStaggeredField, d::Int)
    return Base.OneTo(size(field, d))
end

function Base.show(io::IO, field::AbstractStaggeredField)
    nx, ny = field.grid.data.nx, field.grid.data.ny
    return print(io, "$(typeof(field))(size = ($nx, $ny), type = $(eltype(field.values)))")
end

function set!(field::AbstractStaggeredField, v::Real)
    return fill!(field.values, v)
end

function set!(field::AbstractStaggeredField, f::Function)
    I, J = indices(field.grid)
    for i in I, j in J
        field[i, j] = f(field.grid[i, j]...)
    end
    return field
end

struct UField <: AbstractStaggeredField
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

struct VField <: AbstractStaggeredField
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

struct PField <: AbstractStaggeredField
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

struct StaggeredGridFields
    grid_data::GridData
    u::UField
    v::VField
    p::PField

    function StaggeredGridFields(grid_data::GridData)
        u = UField(grid_data)
        v = VField(grid_data)
        p = PField(grid_data)
        return new(grid_data, u, v, p)
    end
end

function StaggeredGridFields(grid_data::GridData, u_init::Float64, v_init::Float64,
                             p_init::Float64)
    fields = StaggeredGridFields(grid_data)
    set!(fields.u, u_init)
    set!(fields.v, v_init)
    set!(fields.p, p_init)
    return fields
end

function StaggeredGridFields(grid_data::GridData, init::Float64)
    return StaggeredGridFields(grid_data, init, init, init)
end
