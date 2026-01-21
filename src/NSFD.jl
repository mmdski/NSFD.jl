# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

module NSFD

export
       Geometry, CC, EC, CN, EN,
       fill_fg!, fill_rhs!, fill_uv_next!, StaggeredValue,
       StaggeredField, set!, maxabs,
       FlowFields, FG, interp,
       PressureIterator, DomainBC, Domain, East, West, North, South, NoSlip, apply!,
       compute_δt, compute_f, compute_g, compute_u_next, compute_v_next, compute_rhs,
       compute_p_it, compute_p_res

# cell position
abstract type CellPos end

struct CC <: CellPos end  # center, center. location of pressure
struct EC <: CellPos end  # east, center. location of u component of velocity
struct CN <: CellPos end  # center, north. location of v component of velocity
struct EN <: CellPos end  # east, north. used for interpolated u, v values

include("staggered_value.jl")
include("staggered_field.jl")
include("operators.jl")

# geometry
struct Geometry
    x_length::Float64
    y_length::Float64
    imax::Int
    jmax::Int
    δx::Float64
    δy::Float64
    function Geometry(x_length::Float64, y_length::Float64, imax::Int, jmax::Int)
        return new(x_length, y_length, imax, jmax, x_length / imax, y_length / jmax)
    end
end

struct FlowFields
    u::StaggeredField{EC}
    v::StaggeredField{CN}
    p::StaggeredField{CC}
    function FlowFields(imax::Int, jmax::Int)
        return new(StaggeredField(EC, imax, jmax),
                   StaggeredField(CN, imax, jmax),
                   StaggeredField(CC, imax, jmax))
    end
end

struct FG
    f::StaggeredField{EC}
    g::StaggeredField{CN}
    function FG(imax::Int, jmax::Int)
        return new(StaggeredField(EC, imax, jmax), StaggeredField(CN, imax, jmax))
    end
end

function set!(uvp::FlowFields, u::Float64, v::Float64, p::Float64)
    set!(uvp.u, u)
    set!(uvp.v, v)
    set!(uvp.p, p)
    return uvp
end

function compute_δt(uvp::FlowFields, geom::Geometry, Re::Float64, τ::Float64)
    return τ * min(Re / 2 * (1 / (geom.δx * geom.δx) +
                             1 / (geom.δy * geom.δy))^-1,
                   geom.δx / maxabs(uvp.u),
                   geom.δy / maxabs(uvp.v))
end

function compute_f(u::StaggeredField{EC}, v::StaggeredField{CN},
                   δx::Float64, δy::Float64, δt::Float64,
                   GX::Float64, Re::Float64, γ::Float64,
                   i::Int, j::Int)
    return u[i, j] +
           δt * (1.0 / Re * (lap(u, δx, δy, i, j)) -
                 advect_u(u, v, δx, δy, γ, i, j) + GX)
end

function compute_g(u::StaggeredField{EC}, v::StaggeredField{CN},
                   δx::Float64, δy::Float64, δt::Float64,
                   GY::Float64, Re::Float64, γ::Float64,
                   i::Int, j::Int)
    return v[i, j] +
           δt * (1.0 / Re * (lap(v, δx, δy, i, j)) -
                 advect_v(u, v, δx, δy, γ, i, j) + GY)
end

function compute_rhs(f::StaggeredField{EC}, g::StaggeredField{CN},
                     δx::Float64, δy::Float64, δt::Float64,
                     i::Int, j::Int)
    return 1.0 / δt * div(f, g, δx, δy, i, j)
end

function compute_p_it(p::StaggeredField{CC}, rhs::StaggeredField{CC},
                      δx::Float64, δy::Float64, ω::Float64,
                      i::Int, j::Int)
    return (1.0 - ω) * p[i, j] +
           ω / (2.0 / δx^2 + 2.0 / δy^2) *
           ((p[i + 1, j] + p[i - 1, j]) / δx^2 +
            (p[i, j + 1] + p[i, j - 1]) / δy^2 - rhs[i, j])
end

function compute_p_res(p_it::StaggeredField{CC}, rhs::StaggeredField,
                       δx::Float64, δy::Float64,
                       i::Int, j::Int)
    return lap(p_it, δx, δy, i, j) - rhs[i, j]
end

function compute_u_next(f::StaggeredField{EC}, p::StaggeredField{CC},
                        δx::Float64, δt::Float64,
                        i::Int, j::Int)
    return f[i, j] - δt * ∂x(p, δx, i, j)
end

function compute_v_next(g::StaggeredField{CN}, p::StaggeredField{CC},
                        δy::Float64, δt::Float64,
                        i::Int, j::Int)
    return g[i, j] - δt * ∂y(p, δy, i, j)
end

function fill_fg!(fg::FG, uvp::FlowFields, geom::Geometry,
                  δt::Float64, GX::Float64, GY::Float64,
                  Re::Float64, γ::Float64)
    for j in 1:(geom.jmax)
        for i in 1:(geom.imax - 1)
            fg.f[i, j] = compute_f(uvp.u, uvp.v, geom.δx, geom.δy,
                                   δt, GX, Re, γ, i, j)
        end
    end

    for j in 1:(geom.jmax - 1)
        for i in 1:(geom.imax)
            fg.g[i, j] = compute_g(uvp.u, uvp.v, geom.δx, geom.δy,
                                   δt, GY, Re, γ, i, j)
        end
    end

    apply!(fg, uvp)

    return fg
end

function fill_rhs!(rhs::StaggeredField{CC}, fg::FG,
                   geom::Geometry, δt::Float64)
    for j in 1:(geom.jmax)
        for i in 1:(geom.imax)
            rhs[i, j] = compute_rhs(fg.f, fg.g, geom.δx, geom.δy, δt, i, j)
        end
    end
end

function fill_uv_next!(uvp::FlowFields, fg::FG, geom::Geometry, δt::Float64)
    for j in 1:(geom.jmax)
        for i in 1:(geom.imax - 1)
            uvp.u[i, j] = compute_u_next(fg.f, uvp.p, geom.δx, δt, i, j)
        end
    end
    for j in 1:(geom.jmax - 1)
        for i in 1:(geom.imax)
            uvp.v[i, j] = compute_v_next(fg.g, uvp.p, geom.δy, δt, i, j)
        end
    end
end

include("iter_pressure.jl")
include("interp.jl")
include("domain_bc.jl")

end
