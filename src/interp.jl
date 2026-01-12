# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

function interp(f::StaggeredField{EC}, ::Type{CC}, i::Int, j::Int)
    return StaggeredValue{CC}((f[i, j].value + f[i - 1, j].value) / 2.0)
end

function interp(f::StaggeredField{EC}, ::Type{EN}, i::Int, j::Int)
    return StaggeredValue{EN}((f[i, j].value + f[i, j + 1].value) / 2.0)
end

function interp(f::StaggeredField{CN}, ::Type{CC}, i::Int, j::Int)
    return StaggeredValue{CC}((f[i, j].value + f[i, j - 1].value) / 2.0)
end

function interp(f::StaggeredField{CN}, ::Type{EN}, i::Int, j::Int)
    return StaggeredValue{EN}((f[i, j].value + f[i + 1, j].value) / 2.0)
end

function interp(u::StaggeredField{EC},
                δx::Float64, δy::Float64,
                x::Float64, y::Float64)
    i = Int(floor(x / δx)) + 1
    j = Int(floor((y + δy / 2) / δy)) + 1
    x₁ = (i - 1.0) * δx
    x₂ = i * δx
    y₁ = ((j - 1.0) - 0.5) * δy
    y₂ = (j - 0.5) * δy
    u₁ = u[i - 1, j - 1].value
    u₂ = u[i, j - 1].value
    u₃ = u[i - 1, j].value
    u₄ = u[i, j].value
    return 1 / (δx * δy) *
           ((x₂ - x) * (y₂ - y) * u₁ +
            (x - x₁) * (y₂ - y) * u₂ +
            (x₂ - x) * (y - y₁) * u₃ +
            (x - x₁) * (y - y₁) * u₄)
end

function interp(v::StaggeredField{CN},
                δx::Float64, δy::Float64,
                x::Float64, y::Float64)
    i = Int(floor((x + δx / 2) / δx)) + 1
    j = Int(floor(y / δy)) + 1
    x₁ = ((i - 1.0) - 0.5) * δx
    x₂ = (i - 0.5) * δx
    y₁ = (j - 1.0) * δy
    y₂ = j * δy
    v₁ = v[i - 1, j - 1].value
    v₂ = v[i, j - 1].value
    v₃ = v[i - 1, j].value
    v₄ = v[i, j].value
    return 1 / (δx * δy) *
           ((x₂ - x) * (y₂ - y) * v₁ +
            (x - x₁) * (y₂ - y) * v₂ +
            (x₂ - x) * (y - y₁) * v₃ +
            (x - x₁) * (y - y₁) * v₄)
end
