# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# StaggeredField
abstract type AbstractStaggeredField end

function Base.getindex(field::AbstractStaggeredField, i::Int, j::Int)
    return field.values[i + 1, j + 1]
end

function Base.setindex!(field::AbstractStaggeredField, v::Real, i::Int, j::Int)
    return (field.values[i + 1, j + 1] = v)
end

function Base.show(io::IO, field::AbstractStaggeredField)
    nx, ny = field.grid.data.nx, field.grid.data.ny
    return print(io, "$(typeof(field))(size = ($nx, $ny), type = $(eltype(field.values)))")
end

function set!(field::AbstractStaggeredField, v::Real)
    I, J = interior_indices(field.grid)
    for i in I, j in J
        field[i, j] = v
    end
    return field
end

function set!(field::AbstractStaggeredField, f::Function)
    I, J = interior_indices(field.grid)
    for i in I, j in J
        field[i, j] = f(field.grid[i, j]...)
    end
    return field
end

struct UField <: AbstractStaggeredField
    grid::UGrid
    values::Matrix{Float64}

    function UField(data::GridData, init::Float64)
        grid = UGrid(data)
        values = fill(init, data.nx + 2, data.ny + 2)
        return new(grid, values)
    end
end

function UField(data::GridData)
    return UField(data, 0.0)
end

struct VField <: AbstractStaggeredField
    grid::VGrid
    values::Matrix{Float64}

    function VField(data::GridData, init::Float64)
        grid = VGrid(data)
        values = fill(init, data.nx + 2, data.ny + 2)
        return new(grid, values)
    end
end

function VField(data::GridData)
    return VField(data, 0.0)
end

struct PField <: AbstractStaggeredField
    grid::PGrid
    values::Matrix{Float64}

    function PField(data::GridData, init::Float64)
        grid = PGrid(data)
        values = fill(init, data.nx + 2, data.ny + 2)
        return new(grid, values)
    end
end

function PField(data::GridData)
    return PField(data, 0.0)
end
