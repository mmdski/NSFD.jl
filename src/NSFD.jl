# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

module NSFD

export
       Geometry, CC, EC, CN, EN,
       fill_fg!, fill_rhs!, fill_uv_next!, StaggeredValue,
       StaggeredField, set!, maxabs,
       FlowFields, FG, interp,
       PressureIterator, DomainBC, Domain, East, West, North, South, NoSlip, apply!,
       compute_δt, compute_f, compute_g, compute_u_next, compute_v_next, compute_rhs,
       compute_p_it, compute_p_res

# cell position
abstract type CellPos end

struct CC <: CellPos end  # center, center. location of pressure
struct EC <: CellPos end  # east, center. location of u component of velocity
struct CN <: CellPos end  # center, north. location of v component of velocity
struct EN <: CellPos end  # east, north. used for interpolated u, v values

include("staggered_value.jl")
include("staggered_field.jl")
include("operators.jl")
include("geometry.jl")
include("fg.jl")
include("flow_fields.jl")
include("compute.jl")
include("fill.jl")
include("iter_pressure.jl")
include("interp.jl")
include("domain_bc.jl")

end
