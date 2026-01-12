# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

function ∂x(f::StaggeredField{EC}, δx::Float64, i::Int, j::Int)
    return StaggeredValue{CC}(((f[i, j] - f[i - 1, j]) / δx).value)
end

function ∂x(p::StaggeredField{CC}, δx::Float64, i::Int, j::Int)
    return StaggeredValue{EC}(((p[i + 1, j] - p[i, i]) / δx).value)
end

function ∂y(g::StaggeredField{CN}, δy::Float64, i::Int, j::Int)
    return StaggeredValue{CC}(((g[i, j] - g[i, j - 1]) / δy).value)
end

function ∂y(p::StaggeredField{CC}, δy::Float64, i::Int, j::Int)
    return StaggeredValue{CN}(((p[i, j + 1] - p[i, i]) / δy).value)
end

function div(f::StaggeredField{EC}, g::StaggeredField{CN}, δx::Float64, δy::Float64, i::Int,
             j::Int)
    return ∂x(f, δx, i, j) + ∂y(g, δy, i, j)
end

function ∂x_cd(f::StaggeredField{EC}, u::StaggeredField{EC}, δx::Float64, i::Int, j::Int)
    return StaggeredValue{EC}((1.0 / δx *
                               (interp(u, CC, i, j) * interp(f, CC, i, j) -
                                interp(u, CC, i, j) * interp(f, CC, i, j))).value)
end

function ∂x_cd(u::StaggeredField{EC}, f::StaggeredField{CN}, δx::Float64, i::Int, j::Int)
    return StaggeredValue{CN}((1.0 / δx *
                               (interp(u, EN, i, j) * interp(f, EN, i, j) -
                                interp(u, EN, i, j) * interp(f, EN, i, j))).value)
end

function ∂y_cd(f::StaggeredField{EC}, v::StaggeredField{CN}, δy::Float64, i::Int, j::Int)
    return StaggeredValue{EC}((1.0 / δy *
                               (interp(v, EN, i, j) * interp(f, EN, i, j) -
                                interp(v, EN, i, j) * interp(f, EN, i, j))).value)
end

function ∂y_cd(f::StaggeredField{CN}, v::StaggeredField{CN}, δy::Float64, i::Int, j::Int)
    return StaggeredValue{CN}((1.0 / δy *
                               (interp(v, CC, i, j) * interp(f, CC, i, j) -
                                interp(v, CC, i, j) * interp(f, CC, i, j))).value)
end

function advect_u_cd(u::StaggeredField{EC}, v::StaggeredField{CN},
                     δx::Float64, δy::Float64,
                     i::Int, j::Int)
    return ∂x_cd(u, u, δx, i, j) + ∂y_cd(u, v, δy, i, j)
end

function advect_v_cd(u::StaggeredField{EC}, v::StaggeredField{CN},
                     δx::Float64, δy::Float64,
                     i::Int, j::Int)
    return ∂x_cd(u, v, δx, i, j) + ∂y_cd(u, v, δy, i, j)
end

function ∂x_dc(u::StaggeredField{EC}, f::StaggeredField{EC}, δx::Float64,
               i::Int, j::Int)
    kᵣ = interp(u, CC, i + 1, j).value
    kₗ = interp(u, CC, i, j).value
    return ((kᵣ - abs(kᵣ)) * f[i + 1, j] + (kᵣ + abs(kᵣ) - kₗ + abs(kₗ))f[i, j] +
            (-kₗ - abs(kₗ)) * f[i - 1, j]) / (2 * δx)
end

function ∂y_dc(f::StaggeredField{EC}, v::StaggeredField{CN}, δy::Float64,
               i::Int, j::Int)
    kᵣ = interp(v, EN, i, j).value
    kₗ = interp(v, EN, i, j - 1).value
    return ((kᵣ - abs(kᵣ)) * f[i, j + 1] + (kᵣ + abs(kᵣ) - kₗ + abs(kₗ))f[i, j] +
            (-kₗ - abs(kₗ)) * f[i, j - 1]) / (2 * δy)
end

function ∂x_dc(u::StaggeredField{EC}, f::StaggeredField{CN}, δx::Float64,
               i::Int, j::Int)
    kᵣ = interp(u, EN, i, j).value
    kₗ = interp(u, EN, i - 1, j).value
    return ((kᵣ - abs(kᵣ)) * f[i + 1, j] + (kᵣ + abs(kᵣ) - kₗ + abs(kₗ))f[i, j] +
            (-kₗ - abs(kₗ)) * f[i - 1, j]) / (2 * δx)
end

function ∂y_dc(f::StaggeredField{CN}, v::StaggeredField{CN}, δy::Float64,
               i::Int, j::Int)
    kᵣ = interp(v, CC, i, j).value
    kₗ = interp(v, CC, i, j).value
    return ((kᵣ - abs(kᵣ)) * f[i, j + 1] + (kᵣ + abs(kᵣ) - kₗ + abs(kₗ))f[i, j] +
            (-kₗ - abs(kₗ)) * f[i, j - 1]) / (2 * δy)
end

function advect_u_dc(u::StaggeredField{EC}, v::StaggeredField{CN},
                     δx::Float64, δy::Float64,
                     i::Int, j::Int)
    return ∂x_dc(u, u, δx, i, j) + ∂y_dc(u, v, δy, i, j)
end

function advect_v_dc(u::StaggeredField{EC}, v::StaggeredField{CN},
                     δx::Float64, δy::Float64,
                     i::Int, j::Int)
    return ∂x_dc(u, v, δx, i, j) + ∂y_dc(u, v, δy, i, j)
end

function advect_u(u::StaggeredField{EC}, v::StaggeredField{CN},
                  δx::Float64, δy::Float64, γ::Float64,
                  i::Int, j::Int)
    return γ * advect_u_cd(u, v, δx, δy, i, j) +
           (1 - γ) * advect_u_dc(u, v, δx, δy, i, j)
end

function advect_v(u::StaggeredField{EC}, v::StaggeredField{CN},
                  δx::Float64, δy::Float64, γ::Float64,
                  i::Int, j::Int)
    return γ * advect_v_cd(u, v, δx, δy, i, j) +
           (1 - γ) * advect_v_dc(u, v, δx, δy, i, j)
end

function ∂x²(f::StaggeredField, δx::Float64, i::Int, j::Int)
    return (f[i + 1, j] - 2 * f[i, j] + f[i - 1, j]) / (δx * δx)
end

function ∂y²(f::StaggeredField, δy::Float64, i::Int, j::Int)
    return (f[i, j + 1] - 2 * f[i, j] + f[i, j - 1]) / (δy * δy)
end

function lap(f::StaggeredField, δx::Float64, δy::Float64, i::Int, j::Int)
    return ∂x²(f, δx, i, j) + ∂y²(f, δy, i, j)
end
