# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

function ∂x(f::StaggeredField{CC}, δx::Float64, i::Int, j::Int)
    return StaggeredValue{EC}((f[i + 1, j] - f[i, j]).value / δx)
end

function ∂y(f::StaggeredField{CC}, δy::Float64, i::Int, j::Int)
    return StaggeredValue{CN}((f[i, j + 1] - f[i, j]).value / δy)
end

function ∂x(f::StaggeredField{EC}, δx::Float64, i::Int, j::Int)
    return StaggeredValue{CC}((f[i, j] - f[i - 1, j]).value / δx)
end

function ∂y(g::StaggeredField{CN}, δy::Float64, i::Int, j::Int)
    return StaggeredValue{CC}((g[i, j] - g[i, j - 1]).value / δy)
end

function ∂x(f::StaggeredField,
            kr::StaggeredValue, kl::StaggeredValue,
            δx::Float64, γ::Float64,
            i::Int, j::Int)
    return 1.0 / δx * ((kr * (f[i, j] + f[i + 1, j]) / 2.0 -
                        kl * (f[i - 1, j] + f[i, j]) / 2.0) +
                       γ * (abs(kr) * (f[i, j] - f[i + 1, j]) / 2.0 -
                            abs(kl) * (f[i - 1, j] - f[i, j]) / 2.0))
end

function ∂x(f::StaggeredField{EC}, u::StaggeredField{EC},
            δx::Float64, γ::Float64,
            i::Int, j::Int)
    kr = interp(u, CC, i + 1, j)
    kl = interp(u, CC, i, j)
    return ∂x(f, kr, kl, δx, γ, i, j)
end

function ∂x(u::StaggeredField{EC}, f::StaggeredField{CN},
            δx::Float64, γ::Float64,
            i::Int, j::Int)
    kr = interp(u, EN, i, j)
    kl = interp(u, EN, i - 1, j)
    return ∂x(f, kr, kl, δx, γ, i, j)
end

function ∂y(f::StaggeredField, ku::StaggeredValue, kd::StaggeredValue,
            δy::Float64, γ::Float64,
            i::Int, j::Int)
    return 1.0 / δy * ((ku * (f[i, j] + f[i, j + 1]) / 2.0 -
                        kd * (f[i, j - 1] + f[i, j]) / 2.0) +
                       γ * (abs(ku) * (f[i, j] - f[i, j + 1]) / 2.0 -
                            abs(kd) * (f[i, j - 1] - f[i, j]) / 2.0))
end

function ∂y(f::StaggeredField{EC}, v::StaggeredField{CN},
            δy::Float64, γ::Float64,
            i::Int, j::Int)
    ku = interp(v, EN, i, j)
    kd = interp(v, EN, i, j - 1)
    return ∂y(f, ku, kd, δy, γ, i, j)
end

function ∂y(f::StaggeredField{CN}, v::StaggeredField{CN},
            δy::Float64, γ::Float64,
            i::Int, j::Int)
    ku = interp(v, CC, i, j + 1)
    kd = interp(v, CC, i, j)
    return ∂y(f, ku, kd, δy, γ, i, j)
end

function ∂x²(f::StaggeredField, δx, i::Int, j::Int)
    return (f[i + 1, j] - 2.0 * f[i, j] + f[i - 1, j]) / δx^2
end

function ∂y²(f::StaggeredField, δy, i::Int, j::Int)
    return (f[i, j + 1] - 2.0 * f[i, j] + f[i, j - 1]) / δy^2
end

function lap(f::StaggeredField, δx::Float64, δy::Float64,
             i::Int, j::Int)
    return ∂x²(f, δx, i, j) + ∂y²(f, δy, i, j)
end

function div(f::StaggeredField{EC}, g::StaggeredField{CN},
             δx::Float64, δy::Float64,
             i::Int, j::Int)
    return ∂x(f, δx, i, j) + ∂y(g, δy, i, j)
end

function advect_u(u::StaggeredField{EC}, v::StaggeredField{CN},
                  δx::Float64, δy::Float64, γ::Float64,
                  i::Int, j::Int)
    return ∂x(u, u, δx, γ, i, j) + ∂y(u, v, δy, γ, i, j)
end

function advect_v(u::StaggeredField{EC}, v::StaggeredField{CN},
                  δx::Float64, δy::Float64, γ::Float64,
                  i::Int, j::Int)
    return ∂x(u, v, δx, γ, i, j) + ∂y(v, v, δy, γ, i, j)
end
