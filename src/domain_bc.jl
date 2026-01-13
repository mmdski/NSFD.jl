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
    for j in axes(u, 2)
        u[end, j] = 0.0
        v[end + 1, j] = 2 * bc.ut - v[end, j]
    end
end

function apply!(bc::NoSlip{Domain,West}, u::StaggeredField{EC}, v::StaggeredField{CN})
    for j in axes(u, 2)
        u[0, j] = 0.0
        v[0, j] = 2 * bc.ut - v[1, j]
    end
end

function apply!(bc::NoSlip{Domain,North}, u::StaggeredField{EC}, v::StaggeredField{CN})
    for i in axes(u, 1)
        u[i, end + 1] = 2 * bc.ut - u[i, end]
        v[i, end] = 0.0
    end
end

function apply!(bc::NoSlip{Domain,South}, u::StaggeredField{EC}, v::StaggeredField{CN})
    for i in axes(u, 1)
        u[i, 0] = 2 * bc.ut - u[i, 1]
        v[i, 0] = 0.0
    end
end

function apply!(nbc::BoundaryCondition{Domain,North},
                sbc::BoundaryCondition{Domain,South},
                ebc::BoundaryCondition{Domain,East},
                wbc::BoundaryCondition{Domain,West},
                u::StaggeredField{EC}, v::StaggeredField{CN})
    apply!(nbc, u, v)
    apply!(sbc, u, v)
    apply!(ebc, u, v)
    apply!(wbc, u, v)
    return
end

function apply!(p::StaggeredField{CC})
    for i in axes(p, 1)
        p[i, 0] = p[i, 1]
        p[i, end + 1] = p[i, end]
    end
    for j in axes(p, 2)
        p[0, j] = p[1, j]
        p[end + 1, j] = p[end, j]
    end
end

function apply!(f::StaggeredField{EC}, u::StaggeredField{EC})
    for j in axes(f, 2)
        f[0, j] = u[0, j]
        f[end, j] = u[end, j]
    end
end

function apply!(g::StaggeredField{CN}, v::StaggeredField{CN})
    for i in axes(g, 1)
        g[i, 0] = v[i, 0]
        g[i, end] = v[i, end]
    end
end

function apply!(f::StaggeredField{EC}, g::StaggeredField{CN},
                u::StaggeredField{EC}, v::StaggeredField{CN})
    apply!(f, u)
    apply!(g, v)
    return
end
