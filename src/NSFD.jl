# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

module NSFD

# Export
export GridData, StaggeredGrid, meshgrid
export StaggeredVec, StaggeredVecField, UField, VField, PField, set!
export InterpolatedVelocity
export State

# boundary conditions
export apply!,
# North
       NOutflowBCond, NNoSlipBCond,
# South
       SOutflowBCond, SNoSlipBCond,
# East
       EOutflowBCond, ENoSlipBCond,
# West
       WOutflowBCond, WNoSlipBCond

include("grid.jl")
include("staggered_fields.jl")

struct State
    grid_data::GridData
    u::StaggeredVecField
    p::PField
    function State(grid_data::GridData)
        return new(grid_data, StaggeredVecField(grid_data), PField(grid_data))
    end
end

function Base.show(io::IO, state::State)
    (; nx, ny) = state.p.grid.data
    return print(io, "$(typeof(state))(size = ($nx, $ny))")
end

include("interp_fields.jl")
include("global_bcond.jl")

include("time_step.jl")

include("Plotting.jl")

using .Plotting: plot!
export plot!

end # module
