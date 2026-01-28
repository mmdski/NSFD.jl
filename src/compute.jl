# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
