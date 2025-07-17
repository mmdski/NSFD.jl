# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

module NSFD

# Export
export GridData, StaggeredGrid, meshgrid
export StaggeredGridFields, set!
export InterpolatedVelocity

include("grid.jl")
include("staggered_fields.jl")
include("interp_fields.jl")

include("Plotting.jl")

using .Plotting: plot!
export plot!

end # module
