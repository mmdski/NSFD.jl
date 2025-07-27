# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct RHS
    ∇::Divergence
    values::PField
end

function RHS(state::State)
    (; Δx, Δy) = state.geometry.grid_data
    ∇ = Divergence(Δx, Δy)
    values = PField(state.geometry.grid_data, 0.0)
    return RHS(∇, values)
end

function Base.getindex(rhs::RHS, i::Int, j::Int)
    return rhs.values[i, j]
end

function Base.setindex!(rhs::RHS, s::Float64, i::Int, j::Int)
    return (rhs.values[i, j] = s)
end

function compute_rhs!(rhs::RHS, fg::FG, state::State)
    for I in state.geometry.fluid_cells
        i, j = Tuple(I)
        rhs[i, j] = 1.0 / state.time_data.Δt * rhs.∇(fg.values, i, j)
    end
    return rhs
end
