# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

function interpolate_u(u::UField, x::Float64, y::Float64)
    Δx = u.grid.data.Δx
    Δy = u.grid.data.Δy
    i = Int(floor(x / Δx)) + 1
    j = Int(floor((y + Δy / 2) / Δy)) + 1
    x₁ = (i - 1.0)Δx
    x₂ = i * Δx
    y₁ = ((j - 1.0) - 0.5)Δy
    y₂ = (j - 0.5)Δy
    u₁ = u[i - 1, j - 1]
    u₂ = u[i, j - 1]
    u₃ = u[i - 1, j]
    u₄ = u[i, j]
    return 1 / (Δx * Δy) *
           ((x₂ - x) * (y₂ - y)u₁ + (x - x₁) * (y₂ - y)u₂ + (x₂ - x) * (y - y₁)u₃ +
            (x - x₁) * (y - y₁)u₄)
end

function interpolate_v(v::VField, x::Float64, y::Float64)
    Δx = v.grid.data.Δx
    Δy = v.grid.data.Δy
    i = Int(floor((x + Δx / 2) / Δx)) + 1
    j = Int(floor(y / Δy)) + 1
    x₁ = ((i - 1) - 0.5)Δx
    x₂ = (i - 0.5)Δx
    y₁ = (j - 1)Δy
    y₂ = j * Δy
    v₁ = v[i - 1, j - 1]
    v₂ = v[i, j - 1]
    v₃ = v[i - 1, j]
    v₄ = v[i, j]
    return 1 / (Δx * Δy) *
           ((x₂ - x) * (y₂ - y)v₁ + (x - x₁) * (y₂ - y)v₂ + (x₂ - x) * (y - y₁)v₃ +
            (x - x₁) * (y - y₁)v₄)
end

abstract type AbstractInterpolatedField end

struct InterpolatedVelocity <: AbstractInterpolatedField
    u::UField
    v::VField
end

function InterpolatedVelocity(fields::StaggeredFields)
    return InterpolatedVelocity(fields.u, fields.v)
end

function (u::InterpolatedVelocity)(x, y)
    return interpolate_u(u.u, x, y), interpolate_v(u.v, x, y)
end
