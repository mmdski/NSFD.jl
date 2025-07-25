# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

abstract type AbstractGlobalBCond end

# North boundary conditions
abstract type NBCond <: AbstractGlobalBCond end

struct NOutflowBCond <: NBCond end

function apply!(::NOutflowBCond, state::State)
    for i in 1:(state.grid_data.nx)
        state.u.x[i, end + 1] = state.u.x[i, end]
        state.u.y[i, end] = state.u.y[i, end - 1]
    end
    return state
end

struct NNoSlipBCond <: NBCond
    u_wall::Float64
end

function NNoSlipBCond()
    return NNoSlipBCond(0.0)
end

function apply!(bc::NNoSlipBCond, state::State)
    for i in 1:(state.grid_data.nx)
        state.u.x[i, end + 1] = 2 * bc.u_wall - state.u.x[i, end]
        state.u.y[i, end] = 0
    end
    return state
end

# South boundary conditions
abstract type SBCond <: AbstractGlobalBCond end

struct SOutflowBCond <: SBCond end

function apply!(::SOutflowBCond, state::State)
    for i in 1:(state.grid_data.nx)
        state.u.x[i, 0] = state.u.x[i, 1]
        state.u.y[i, 0] = state.u.y[i, 1]
    end
    return state
end

struct SNoSlipBCond <: SBCond
    u_wall::Float64
end

function SNoSlipBCond()
    return SNoSlipBCond(0.0)
end

function apply!(bc::SNoSlipBCond, state::State)
    for i in 1:(state.grid_data.nx)
        state.u.x[i, 0] = 2 * bc.u_wall - state.u.x[i, 1]
        state.u.y[i, 0] = 0
    end
    return state
end

# East boundary conditions
abstract type EBCond <: AbstractGlobalBCond end

struct EOutflowBCond <: EBCond end

function apply!(::EOutflowBCond, state::State)
    for j in 1:(state.grid_data.ny)
        state.u.x[end, j] = state.u.x[end - 1, j]
        state.u.y[end + 1, j] = state.u.y[end, j]
    end
    return state
end

struct ENoSlipBCond <: EBCond
    v_wall::Float64
end

function ENoSlipBCond()
    return ENoSlipBCond(0.0)
end

function apply!(bc::ENoSlipBCond, state::State)
    for j in 1:(state.grid_data.ny)
        state.u.x[end, j] = 0
        state.u.y[end + 1, j] = 2 * bc.v_wall - state.u.y[end, j]
    end
    return state
end

# West boundary conditions
abstract type WBCond <: AbstractGlobalBCond end

struct WOutflowBCond <: WBCond end

function apply!(::WOutflowBCond, state::State)
    for j in 1:(state.grid_data.ny)
        state.u.x[0, j] = state.u.x[1, j]
        state.u.y[0, j] = state.u.y[1, j]
    end
    return state
end

struct WNoSlipBCond <: WBCond
    v_wall::Float64
end

function WNoSlipBCond()
    return WNoSlipBCond(0.0)
end

function apply!(bc::WNoSlipBCond, state::State)
    for j in 1:(state.grid_data.ny)
        state.u.x[0, j] = 0
        state.u.y[0, j] = 2 * bc.v_wall - state.u.y[1, j]
    end
    return state
end

struct GlobalBCond
    north::NBCond
    south::SBCond
    east::EBCond
    west::WBCond
end

function apply!(bc::GlobalBCond, state::State)
    apply!(bc.north, state)
    apply!(bc.south, state)
    apply!(bc.east, state)
    apply!(bc.west, state)
    return state
end

function Base.show(io::IO, bc::GlobalBCond)
    println(io, "Global Boundary Conditions:")
    println(io, "  North: ", typeof(bc.north))
    println(io, "  South: ", typeof(bc.south))
    println(io, "  East:  ", typeof(bc.east))
    println(io, "  West:  ", typeof(bc.west))
    return nothing
end
