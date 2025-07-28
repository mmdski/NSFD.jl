# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct GridData
    x_length::Float64
    y_length::Float64
    nx::Int64
    ny::Int64
    Δx::Float64
    Δy::Float64
    function GridData(x_length, y_length, nx, ny)
        Δx = x_length / nx
        Δy = y_length / ny
        return new(x_length, y_length, nx, ny, Δx, Δy)
    end
end

function GridData(nx::Int, ny::Int, Δx::Float64, Δy::Float64)
    x_length = nx * Δx
    y_length = ny * Δy
    return GridData(x_length, y_length, nx, ny)
end

struct ProblemData
    Re::Float64 # Reynolds number
    G::StaggeredVector # acceleration due to gravity
end

struct TimeSteppingData
    t::Float64 # current time
    T::Float64 # final time
    Δt::Float64 # time step to produce current state
    τ::Float64 # time step safety factor
end

struct PressureIterData
    iter_max::Int # maximal number of pressure iterations in one time step
    ϵₚ::Float64 # stopping tolerance for pressure iteration
    ω::Float64 # relaxation parameter for SOR iteration
    γ::Float64 # upwind differencing factor
end
