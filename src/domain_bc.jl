# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

abstract type BCLocation end
struct Domain <: BCLocation end
struct Internal <: BCLocation end

abstract type BCDirection end
struct East <: BCDirection end
struct West <: BCDirection end
struct North <: BCDirection end
struct South <: BCDirection end

abstract type BoundaryCondition{L<:BCLocation,D<:BCDirection} end

struct NoSlip{L<:BCLocation,D<:BCDirection} <: BoundaryCondition{L,D}
    ut::Float64  # tangential wall speed
end

function NoSlip(::Type{L}, ::Type{D}) where {L<:BCLocation,D<:BCDirection}
    return NoSlip{L,D}(0.0)
end

function NoSlip(::Type{L}, ::Type{D}, ut::Float64) where {L<:BCLocation,D<:BCDirection}
    return NoSlip{L,D}(ut)
end

function apply!(bc::NoSlip{Domain,East}, u::StaggeredField{EC}, v::StaggeredField{CN})
    for j in 1:axes(u, 2)
        u[end, j] = 0.0
        v[end + 1, j] = 2 * bc.ut - v[end, j]
    end
end

function apply!(bc::NoSlip{Domain,West}, u::StaggeredField{EC}, v::StaggeredField{CN})
    for j in 1:axes(u, 2)
        u[0, j] = 0.0
        v[0, j] = 2 * bc.ut - v[1, j]
    end
end

function apply!(bc::NoSlip{Domain,North}, u::StaggeredField{EC}, v::StaggeredField{CN})
    for i in 1:axes(u, 1)
        u[i, end + 1] = 2 * bc.ut - u[i, end]
        v[i, end] = 0.0
    end
end

function apply!(bc::NoSlip{Domain,South}, u::StaggeredField{EC}, v::StaggeredField{CN})
    for i in 1:axes(u, 1)
        u[i, 0] = 2 * bc.ut - u[i, 1]
        v[i, 0] = 0.0
    end
end
