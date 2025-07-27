# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct FG
    Re::Float64
    G::StaggeredVector
    ∇²::Laplacian
    adv::Advection
    values::StaggeredVectorField
end

function FG(nx::Int, ny::Int, Δx::Float64, Δy::Float64, Re::Float64, G::StaggeredVector,
            γ::Float64)
    ∇² = Laplacian(Δx^2, Δy^2)
    adv = Advection(Δx, Δy, γ)
    fg = StaggeredVectorField(GridData(nx, ny, Δx, Δy), (0.0, 0.0))
    return FG(Re, G, ∇², adv, fg)
end

function FG(state::State)
    (; nx, ny, Δx, Δy) = state.geometry.grid_data
    (; Re, G) = state.problem_data
    γ = state.p_iter_data.γ
    return FG(nx, ny, Δx, Δy, Re, G, γ)
end

function Base.getindex(fg::FG, i::Int, j::Int)
    return fg.values[i, j]
end

function Base.setindex!(fg::FG, v::StaggeredVector, i::Int, j::Int)
    return (fg.values[i, j] = v)
end

function Base.size(fg::FG)
    return size(fg.values)
end

function Base.size(fg::FG, d::Int)
    return size(fg.values, d)
end

function Base.axes(fg::FG)
    return Base.OneTo.(size(fg))
end

function Base.axes(fg::FG, d::Int)
    return Base.OneTo(size(fg, d))
end

function compute_fg!(fg::FG, state::State)
    u = state.u
    Δt = state.time_data.Δt
    (; Re, G, ∇², adv) = fg
    for I in state.geometry.fluid_cells
        i, j = Tuple(I)
        fg[i, j] = u[i, j] + Δt * (G + 1.0 / Re * ∇²(u, i, j) - adv(u, u, i, j))
    end
    # set boundary conditions for F, G
    for i in 1:(state.geometry.grid_data.nx)
        fg[i, 0] = StaggeredVector(fg[i, 0].x, u[i, 0].y)
        fg[i, end] = StaggeredVector(fg[i, end].x, u[i, end].y)
    end
    for j in 1:(state.geometry.grid_data.ny)
        fg[0, j] = StaggeredVector(u[0, j].x, fg[0, j].y)
        fg[end, j] = StaggeredVector(u[end, j].x, fg[end, j].y)
    end
    return fg
end
