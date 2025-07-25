# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct FG
    Re::Float64
    G::StaggeredVector
    ∇²::Laplacian
    adv::Advection
    values::Matrix{StaggeredVector}
end

function FG(nx::Int, ny::Int, Δx::Float64, Δy::Float64, Re::Float64, G::StaggeredVector,
            γ::Float64)
    ∇² = Laplacian(Δx^2, Δy^2)
    adv = Advection(Δx, Δy, γ)
    fg = Matrix{StaggeredVector}(undef, nx + 2, ny + 2)
    return FG(Re, G, ∇², adv, fg)
end

function FG(state::State)
    (; nx, ny, Δx, Δy) = state.grid_data
    (; Re, G) = state.problem_data
    γ = state.p_iter_data.γ
    return FG(nx, ny, Δx, Δy, Re, G, γ)
end

function Base.getindex(fg::FG, i::Int, j::Int)
    return fg.values[i + 1, j + 1]
end

function Base.setindex!(fg::FG, v::StaggeredVector, i::Int, j::Int)
    return (fg.values[i + 1, j + 1] = v)
end

function compute_fg!(fg::FG, state::State)
    u = state.u
    Δt = state.time_data.Δt
    (; Re, G, ∇², adv) = fg
    for I in state.geometry.fluid_cells
        i, j = Tuple(I)
        fg[i, j] = u[i, j] + Δt * (G + 1.0 / Re * ∇²(u, i, j) - adv(u, u, i, j))
    end
    return fg
end
