# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

module Plotting
export plot!

using ..NSFD
import ..NSFD: StaggeredFields

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

function plot!(ax::Axis, f::NSFD.AbstractStaggeredField; kwargs...)
    (; Δx, Δy, x_length, y_length) = f.grid.data
    return heatmap!(ax, 0:Δx:x_length, 0:Δy:y_length, f.values[2:(end - 1), 2:(end - 1)];
                    kwargs...)
end

function plot!(ax::Axis, u::NSFD.InterpolatedVelocity, grid_data::GridData, nx=10, ny=10;
               kwargs...)
    (; Δx, Δy, x_length, y_length) = grid_data
    x = range(Δx / 2; stop=x_length - Δx / 2, length=nx)
    y = range(Δy / 2; stop=y_length - Δy / 2, length=ny)
    X = [xi for yi in y, xi in x]
    Y = [yi for yi in y, xi in x]

    F = u.(X, Y)
    U = getindex.(F, 1)
    V = getindex.(F, 2)

    return arrows2d!(vec(X), vec(Y), vec(U), vec(V); kwargs...)
end

function plot!(ax::Axis, f::StaggeredFields, nx::Int, ny::Int)
    h = plot!(ax, f.p)
    u = NSFD.InterpolatedVelocity(f)
    a = plot!(ax, u, f.grid_data, nx, ny)
    plot!(ax, f.grid_data)
    return h, a
end

function plot!(ax::Axis, f::StaggeredFields)
    nx = max(10, round(Int, f.grid_data.nx / 2))
    aspect = f.grid_data.y_length / f.grid_data.x_length
    ny = round(Int, nx * aspect)
    return plot!(ax, f, nx, ny)
end

end # module
