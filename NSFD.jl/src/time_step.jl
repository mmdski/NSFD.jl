# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

function compute_Δt(Δx::Float64, Δy::Float64, Re::Float64, τ::Float64)
    return τ * (Re / 2 * (1 / Δx^2 + 1 / Δy^2)^-1)
end

function compute_Δt(u::StaggeredVectorField, Δx::Float64, Δy::Float64, Re::Float64,
                    τ::Float64)
    u_max_abs = maximum(abs.(@view u.x.values[2:(end - 1), 2:(end - 1)]))
    v_max_abs = maximum(abs.(@view u.y.values[2:(end - 1), 2:(end - 1)]))
    return τ * min(Re / 2 * (1 / Δx^2 + 1 / Δy^2)^-1, Δx / u_max_abs, Δy / v_max_abs)
end
