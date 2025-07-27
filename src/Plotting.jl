# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

module Plotting
export plot!

using ..NSFD
import ..NSFD: State

using CairoMakie

import CairoMakie: plot, plot!

function plot!(ax::Axis, grid_data::GridData; grid_lines=false)
    ax.aspect = DataAspect()
    ax.xgridvisible = false
    ax.ygridvisible = false
    (; Δx, Δy, x_length, y_length) = grid_data
    xmin = -Δx
    xmax = x_length + Δx
    ymin = -Δy
    ymax = y_length + Δy
    if grid_lines
        vlines!(ax, 0:Δx:x_length; ymin=ymin, ymax=ymax, color=(:black, 0.1))
        hlines!(ax, 0:Δy:y_length; xmin=xmin, xmax=xmax, color=(:black, 0.1))
    end
    xlims!(xmin, xmax)
    ylims!(ymin, ymax)
    return ax
end

function plot!(ax::Axis, f::NSFD.AbstractStaggeredScalarField; kwargs...)
    (; Δx, Δy, x_length, y_length) = f.grid.data
    return heatmap!(ax, 0:Δx:x_length, 0:Δy:y_length, f.values[2:(end - 1), 2:(end - 1)];
                    kwargs...)
end

function plot!(ax::Axis, u::NSFD.InterpolatedVelocity, grid_data::GridData, nx=10, ny=10;
               kwargs...)
    (; x_length, y_length) = grid_data
    Δx = x_length / nx
    Δy = y_length / ny
    x = range(Δx / 2; stop=x_length - Δx / 2, length=nx)
    y = range(Δy / 2; stop=y_length - Δy / 2, length=ny)
    X = [xi for yi in y, xi in x]
    Y = [yi for yi in y, xi in x]

    F = u.(X, Y)
    U = getindex.(F, 1)
    V = getindex.(F, 2)

    if !haskey(kwargs, :lengthscale)
        max_mag = maximum(sqrt.(U .^ 2 + V .^ 2))
        lengthscale = 2.5 * min(grid_data.Δx, grid_data.Δy) / max_mag
        return arrows2d!(vec(X), vec(Y), vec(U), vec(V); kwargs..., lengthscale=lengthscale)
    else
        return arrows2d!(vec(X), vec(Y), vec(U), vec(V); kwargs...)
    end
end

function plot!(ax::Axis, state::State, nx::Int, ny::Int)
    h = plot!(ax, state.p)
    u = NSFD.InterpolatedVelocity(state.u)
    a = plot!(ax, u, state.geometry.grid_data, nx, ny)
    plot!(ax, state.geometry.grid_data)
    return h, a
end

function plot!(ax::Axis, state::State)
    nx = max(10, round(Int, state.geometry.grid_data.nx / 2))
    aspect = state.geometry.grid_data.y_length / state.geometry.grid_data.x_length
    if aspect ≤ 1.0
        nx = 10
        ny = round(Int, nx * aspect)
    else
        ny = 10
        nx = round(Int, ny / aspect)
    end
    return plot!(ax, state, nx, ny)
end

end # module
