# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

function fill_fg!(fg::FG, uvp::FlowFields, geom::Geometry,
                  δt::Float64, GX::Float64, GY::Float64,
                  Re::Float64, γ::Float64)
    for j in 1:(geom.jmax)
        for i in 1:(geom.imax - 1)
            fg.f[i, j] = compute_f(uvp.u, uvp.v, geom.δx, geom.δy,
                                   δt, GX, Re, γ, i, j)
        end
    end

    for j in 1:(geom.jmax - 1)
        for i in 1:(geom.imax)
            fg.g[i, j] = compute_g(uvp.u, uvp.v, geom.δx, geom.δy,
                                   δt, GY, Re, γ, i, j)
        end
    end

    apply!(fg, uvp)

    return fg
end

function fill_rhs!(rhs::StaggeredField{CC}, fg::FG,
                   geom::Geometry, δt::Float64)
    for j in 1:(geom.jmax)
        for i in 1:(geom.imax)
            rhs[i, j] = compute_rhs(fg.f, fg.g, geom.δx, geom.δy, δt, i, j)
        end
    end
end
