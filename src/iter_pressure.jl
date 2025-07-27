# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct PressureIterator
    data::PressureIterData
    geometry::Geometry
    rit::PField
end

function PressureIterator(state::State)
    rit = PField(state.geometry.grid_data)
    return PressureIterator(state.p_iter_data, state.geometry, rit)
end

function calc_rit!(p_iter::PressureIterator, p::PField, rhs::RHS)
    ∇² = Laplacian(pit.grid.data)
    for I in p_iter.geometry.fluid_cells
        i, j = Tuple(I)
        p_iter.rit[i, j] = ∇²(p, i, j) - rhs[i, j]
    end
end

function calc_norm!(p_iter::PressureIterator)
    s = 0.0
    n = 0
    for I in p_iter.geometry.fluid_cells
        i, j = Tuple(I)
        s = s + p_iter.rit[i, j]^2
        n = n + 1
    end
    return sqrt(s / n)
end

struct PIterOp
    ω::Float64
    Δx²::Float64
    Δy²::Float64
end

function p_iter_op(op::PIterOp, pᵢₜ::PField, rhs::RHS, i::Int, j::Int)
    (; ω, Δx², Δy²) = op
    return (1.0 - ω) * pᵢₜ[i, j] +
           ω / (2.0 / Δx² + 2.0 / Δy²) *
           ((pᵢₜ[i + 1, j] + pᵢₜ[i - 1, j]) / Δx² + (pᵢₜ[i, j + 1] + pᵢₜ[i, j - 1]) / Δy² -
            rhs[i, j])
end

function apply_p_bc!(p_it::PField)
    for i in 1:(p_it.grid.data.nx)
        p_it[i, 0] = p[i, 1]
        p_it[i, end + 1] = p[i, end]
    end
    for j in 1:(p_it.grid.data.nx)
        p[0, j] = p[1, j]
        p[end + 1, j] = p[end, j]
    end
    return p_it
end
