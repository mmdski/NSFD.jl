# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct StaggeredValue{CellPos}
    value::Float64
end

function Base.show(io::IO, v::StaggeredValue{P}) where {P}
    return print(io, P, "(", v.value, ")")
end

function Base.:-(v::StaggeredValue{P}) where {P}
    return StaggeredValue{P}(-v.value)
end

function Base.:+(c::Float64, v::StaggeredValue{P}) where {P}
    return StaggeredValue{P}(c + v.value)
end

function Base.:-(c::Float64, v::StaggeredValue{P}) where {P}
    return c + -v
end

function Base.:+(v::StaggeredValue{P}, c::Float64) where {P}
    return c + v
end

function Base.:-(v::StaggeredValue{P}, c::Float64) where {P}
    return v + -c
end

function Base.:*(c::Float64, v::StaggeredValue{P}) where {P}
    return StaggeredValue{P}(c * v.value)
end

function Base.:/(c::Float64, v::StaggeredValue{P}) where {P}
    return StaggeredValue{P}(c / v.value)
end

function Base.:*(v::StaggeredValue{P}, c::Float64) where {P}
    return c * v
end

function Base.:/(v::StaggeredValue{P}, c::Float64) where {P}
    return StaggeredValue{P}(v.value / c)
end

function Base.:+(a::StaggeredValue{P}, b::StaggeredValue{P}) where {P}
    return StaggeredValue{P}(a.value + b.value)
end

function Base.:-(a::StaggeredValue{P}, b::StaggeredValue{P}) where {P}
    return StaggeredValue{P}(a.value - b.value)
end

function Base.:*(a::StaggeredValue{P}, b::StaggeredValue{P}) where {P}
    return StaggeredValue{P}(a.value * b.value)
end
