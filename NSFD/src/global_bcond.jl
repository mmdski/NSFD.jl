# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

abstract type AbstractGlobalBCond end

# North boundary conditions
abstract type NBCond <: AbstractGlobalBCond end

struct NOutflowBCond <: NBCond end

function apply!(::NOutflowBCond, fields::StaggeredGridFields)
    for i in 1:(fields.grid_data.nx)
        fields.u[i, end + 1] = fields.u[i, end]
        fields.v[i, end] = fields.v[i, end - 1]
    end
    return fields
end

struct NNoSlipBCond <: NBCond
    u_wall::Float64
end

function NNoSlipBCond()
    return NNoSlipBCond(0.0)
end

function apply!(bc::NNoSlipBCond, fields::StaggeredGridFields)
    for i in 1:(fields.grid_data.nx)
        fields.u[i, end + 1] = 2 * bc.u_wall - fields.u[i, end]
        fields.v[i, end] = 0
    end
    return fields
end

# South boundary conditions
abstract type SBCond <: AbstractGlobalBCond end

struct SOutflowBCond <: SBCond end

function apply!(::SOutflowBCond, fields::StaggeredGridFields)
    for i in 1:(fields.grid_data.nx)
        fields.u[i, 0] = fields.u[i, 1]
        fields.v[i, 0] = fields.v[i, 1]
    end
    return fields
end

struct SNoSlipBCond <: SBCond
    u_wall::Float64
end

function SNoSlipBCond()
    return SNoSlipBCond(0.0)
end

function apply!(bc::SNoSlipBCond, fields::StaggeredGridFields)
    for i in 1:(fields.grid_data.nx)
        fields.u[i, 0] = 2 * bc.u_wall - fields.u[i, 1]
        fields.v[i, 0] = 0
    end
    return fields
end

# East boundary conditions
abstract type EBCond <: AbstractGlobalBCond end

struct EOutflowBCond <: EBCond end

function apply!(::EOutflowBCond, fields::StaggeredGridFields)
    for j in 1:(fields.grid_data.ny)
        fields.u[end, j] = fields.u[end - 1, j]
        fields.v[end + 1, j] = fields.v[end, j]
    end
    return fields
end

struct ENoSlipBCond <: EBCond
    v_wall::Float64
end

function ENoSlipBCond()
    return ENoSlipBCond(0.0)
end

function apply!(bc::ENoSlipBCond, fields::StaggeredGridFields)
    for j in 1:(fields.grid_data.ny)
        fields.u[end, j] = 0
        fields.v[end + 1, j] = 2 * bc.v_wall - fields.v[end, j]
    end
    return fields
end

# West boundary conditions
abstract type WBCond <: AbstractGlobalBCond end

struct WOutflowBCond <: WBCond end

function apply!(::WOutflowBCond, fields::StaggeredGridFields)
    for j in 1:(fields.grid_data.ny)
        fields.u[0, j] = fields.u[1, j]
        fields.v[0, j] = fields.v[1, j]
    end
    return fields
end

struct WNoSlipBCond <: WBCond
    v_wall::Float64
end

function WNoSlipBCond()
    return WNoSlipBCond(0.0)
end

function apply!(bc::WNoSlipBCond, fields::StaggeredGridFields)
    for j in 1:(fields.grid_data.ny)
        fields.u[0, j] = 0
        fields.v[0, j] = 2 * bc.v_wall - fields.v[1, j]
    end
    return fields
end
