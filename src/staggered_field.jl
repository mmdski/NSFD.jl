# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct StaggeredField{CellPos}
    values::Matrix{StaggeredValue{CellPos}}
end

function StaggeredField(::Type{CellPos}, nx::Int, ny::Int) where {CellPos}
    values = Matrix{StaggeredValue{CellPos}}(undef, nx + 2, ny + 2)
    return StaggeredField{CellPos}(values)
end

function Base.getindex(field::StaggeredField, i::Int, j::Int)
    return field.values[i + 1, j + 1]
end

function Base.setindex!(field::StaggeredField{P}, v::StaggeredValue{P}, i::Int,
                        j::Int) where {P}
    return (field.values[i + 1, j + 1] = v)
end

function Base.setindex!(field::StaggeredField{P}, v::Float64, i::Int,
                        j::Int) where {P}
    return (field.values[i + 1, j + 1] = StaggeredValue{P}(v))
end

function Base.show(io::IO, field::StaggeredField)
    return print(io, "$(typeof(field))")
end

function set!(field::StaggeredField{P}, x::StaggeredValue{P}) where {P}
    fill!(field.values, x)
    return field
end

function set!(field::StaggeredField{P}, x::Float64) where {P}
    fill!(field.values, StaggeredValue{P}(x))
    return field
end
