# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# vector laplacian
struct Laplacian
    Δx²::Float64
    Δy²::Float64
end

function Laplacian(grid_data::GridData)
    return Laplacian(grid_data.Δx^2, grid_data.Δy^2)
end

function (op::Laplacian)(a::StaggeredVectorField, i::Int, j::Int)
    ∂x² = (a[i + 1, j] - 2.0 * a[i, j] + a[i - 1, j]) / op.Δx²
    ∂y² = (a[i, j + 1] - 2.0 * a[i, j] + a[i, j - 1]) / op.Δy²
    return ∂x² + ∂y²
end

function (op::Laplacian)(p::PField, i::Int, j::Int)
    ∂x² = (p[i + 1, j] - 2.0 * p[i, j] + p[i - 1, j]) / op.Δx²
    ∂y² = (p[i, j + 1] - 2.0 * p[i, j] + p[i, j - 1]) / op.Δy²
    return ∂x² + ∂y²
end

"""
    Advection

Computes the conservative form of vector advection:
    ∇ · (u ⊗ a)
where `u` is the advecting field and `a` is the advected field.
"""
struct Advection
    Δx::Float64
    Δy::Float64
    γ::Float64
end

function Advection(Δx, Δy, u::StaggeredVectorField)
    return Advection(Δx, Δy, 0.9, u, u)
end

function Advection(grid_data::GridData, u::StaggeredVectorField, γ)
    return Advection(grid_data.Δx, grid_data.Δy, γ, u, u)
end

function Advection(grid_data::GridData, u::StaggeredVectorField)
    return Advection(grid_data, u, 0.9)
end

function advect_x(a::StaggeredVectorField, u::StaggeredVectorField, Δx::Float64,
                  Δy::Float64, γ::Float64, i::Int, j::Int)
    kᵣᵤ = (u.x[i, j] + u.x[i + 1, j]) / 2.0
    kₗᵤ = (u.x[i - 1, j] + u.x[i, j]) / 2.0
    ∂aₓu∂x = 1.0 / Δx *
             ((kᵣᵤ * (a.x[i, j] + a.x[i + 1, j]) / 2.0 -
               kₗᵤ * (a.x[i - 1, j] + a.x[i, j]) / 2.0) +
              γ * (abs(kᵣᵤ) * (a.x[i, j] - a.x[i + 1, j]) / 2.0 -
                   abs(kₗᵤ) * (a.x[i - 1, j] - a.x[i, j]) / 2.0))

    kᵣᵥ = (u.y[i, j] + u.y[i + 1, j]) / 2.0
    kₗᵥ = (u.y[i, j - 1] + u.y[i + 1, j - 1]) / 2.0
    ∂aₓv∂y = 1.0 / Δy * ((kᵣᵥ * (a.x[i, j] + a.x[i, j + 1]) / 2.0 -
                          kₗᵥ * (a.x[i, j - 1] + a.x[i, j]) / 2.0) +
                         γ * (abs(kᵣᵥ) * (a.x[i, j] - a.x[i, j + 1]) / 2.0 -
                              abs(kₗᵥ) * (a.x[i, j - 1] - a.x[i, j]) / 2.0))
    return ∂aₓu∂x + ∂aₓv∂y
end

function advect_y(a::StaggeredVectorField, u::StaggeredVectorField, Δx::Float64,
                  Δy::Float64, γ::Float64, i::Int, j::Int)
    kᵣᵤ = (u.x[i, j] + u.x[i, j + 1]) / 2.0
    kₗᵤ = (u.x[i - 1, j] + u.x[i - 1, j + 1]) / 2.0
    ∂aᵧu∂x = 1.0 / Δx * ((kᵣᵤ * (a.y[i, j] + a.y[i + 1, j]) / 2.0) -
                         kₗᵤ * (a.y[i - 1, j] + a.y[i, j]) / 2.0 +
                         γ * (abs(kᵣᵤ) * (a.y[i, j] - a.y[i + 1, j]) / 2.0 -
                              abs(kₗᵤ) * (a.y[i - 1, j] - a.y[i, j]) / 2.0))

    kᵣᵥ = (u.y[i, j] + u.y[i, j + 1]) / 2.0
    kₗᵥ = (u.y[i, j - 1] + u.y[i, j]) / 2.0
    ∂aᵧv∂y = 1.0 / Δy * ((kᵣᵥ * (a.y[i, j] + a.y[i, j + 1]) / 2.0 -
                          kₗᵥ * (a.y[i, j - 1] + a.y[i, j]) / 2.0) +
                         γ * (abs(kᵣᵥ) * (a.y[i, j] - a.y[i, j + 1]) / 2.0 -
                              abs(kₗᵥ) * (a.y[i, j - 1] - a.y[i, j]) / 2.0))

    return ∂aᵧu∂x + ∂aᵧv∂y
end

function (adv::Advection)(u::StaggeredVectorField, a::StaggeredVectorField, i::Int, j::Int)
    return StaggeredVector(advect_x(a, u, adv.Δx, adv.Δy, adv.γ, i, j),
                           advect_y(a, u, adv.Δx, adv.Δy, adv.γ, i, j))
end

struct Divergence
    Δx::Float64
    Δy::Float64
end

function (op::Divergence)(u::StaggeredVectorField, i::Int, j::Int)
    return (u[i, j].x - u[i - 1, j].x) / op.Δx + (u[i, j].y - u[i, j - 1].y) / op.Δy
end
