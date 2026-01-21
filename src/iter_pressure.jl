# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct PressureIterator
    ω::Float64
    ϵₚ::Float64
    iter_max::Int
    p_res::StaggeredField{CC}
    function PressureIterator(geom::Geometry, ω::Float64, ϵₚ::Float64, iter_max::Int)
        return new(ω, ϵₚ, iter_max, StaggeredField(CC, geom.imax, geom.jmax))
    end
end

# pressure boundary condition

function apply!(p::StaggeredField{CC}, imax::Int, jmax::Int)
    for i in 1:imax
        p[i, 0] = p[i, 1]
        p[i, jmax + 1] = p[i, imax]
    end
    for j in 1:jmax
        p[0, j] = p[1, j]
        p[imax + 1, j] = p[imax, j]
    end
end

function (iter_pressure::PressureIterator)(uvp::FlowFields, rhs::StaggeredField{CC},
                                           geom::Geometry)
    it = 0

    p_norm = maxabs(uvp.p)

    while it < iter_pressure.iter_max
        for j in 1:(geom.jmax)
            for i in 1:(geom.imax)
                uvp.p[i, j] = compute_p_it(uvp.p, rhs,
                                           geom.δx, geom.δy, iter_pressure.ω, i, j)
            end
        end

        for j in 1:(geom.jmax)
            for i in 1:(geom.imax)
                iter_pressure.p_res[i, j] = compute_p_res(uvp.p, rhs,
                                                          geom.δx, geom.δy, i, j)
            end
        end

        if maxabs(iter_pressure.p_res) < iter_pressure.ϵₚ * p_norm
            break
        end

        it = it + 1

        apply!(uvp.p, geom.imax, geom.jmax)
    end

    return it
end
