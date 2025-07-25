# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

module NSFD

# Export
export ProblemData, TimeSteppingData, PressureIterData
export GridData, StaggeredGrid, meshgrid
export StaggeredVector, StaggeredVectorField, UField, VField, PField, set!
export InterpolatedVelocity
export State

# operators
export Advection, Laplacian

# boundary conditions
export apply!, GlobalBCond,
# North
       NOutflowBCond, NNoSlipBCond,
# South
       SOutflowBCond, SNoSlipBCond,
# East
       EOutflowBCond, ENoSlipBCond,
# West
       WOutflowBCond, WNoSlipBCond

include("staggered_vector.jl")
include("data.jl")
include("grid.jl")
include("staggered_fields.jl")
include("geometry.jl")
include("operators.jl")
include("state.jl")
include("time_step.jl")
include("global_bcond.jl")
include("fg.jl")
include("interp_fields.jl")
include("Plotting.jl")

using .Plotting: plot!
export plot!

end # module
