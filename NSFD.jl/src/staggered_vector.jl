# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct StaggeredVector
    x::Float64
    y::Float64
end

function Base.:+(u::StaggeredVector, v::StaggeredVector)
    return StaggeredVector(u.x + v.x, u.y + v.y)
end

function Base.:-(u::StaggeredVector, v::StaggeredVector)
    return StaggeredVector(u.x - v.x, u.y - v.y)
end

function Base.:-(u::StaggeredVector, c::Real)
    return StaggeredVector(u.x - c, u.y - c)
end

function Base.:*(c::Real, u::StaggeredVector)
    return StaggeredVector(c * u.x, c * u.y)
end

function Base.:/(u::StaggeredVector, c::Real)
    return StaggeredVector(u.x / c, u.y / c)
end
