# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct FlowFields
    u::StaggeredField{EC}
    v::StaggeredField{CN}
    p::StaggeredField{CC}
    function FlowFields(imax::Int, jmax::Int)
        return new(StaggeredField(EC, imax, jmax),
                   StaggeredField(CN, imax, jmax),
                   StaggeredField(CC, imax, jmax))
    end
end

function set!(uvp::FlowFields, u::Float64, v::Float64, p::Float64)
    set!(uvp.u, u)
    set!(uvp.v, v)
    set!(uvp.p, p)
    return uvp
end

function fill_uv_next!(uvp::FlowFields, fg::FG, geom::Geometry, δt::Float64)
    for j in 1:(geom.jmax)
        for i in 1:(geom.imax - 1)
            uvp.u[i, j] = compute_u_next(fg.f, uvp.p, geom.δx, δt, i, j)
        end
    end
    for j in 1:(geom.jmax - 1)
        for i in 1:(geom.imax)
            uvp.v[i, j] = compute_v_next(fg.g, uvp.p, geom.δy, δt, i, j)
        end
    end
end
