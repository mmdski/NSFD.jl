# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct State
    grid_data::GridData
    problem_data::ProblemData
    time_data::TimeSteppingData
    p_iter_data::PressureIterData
    geometry::Geometry
    u::StaggeredVectorField
    p::PField
end

function State(grid_data::GridData, problem_data::ProblemData, time_data::TimeSteppingData,
               p_iter_data::PressureIterData)
    u = StaggeredVectorField(grid_data)
    p = StaggeredScalarField(grid_data)
    geometry = Geometry(grid_data)
    return State(grid_data, problem_data, time_data, p_iter_data, geometry, u, p)
end

function Base.show(io::IO, state::State)
    (; nx, ny) = state.p.grid.data
    return print(io, "$(typeof(state))(size = ($nx, $ny))")
end
