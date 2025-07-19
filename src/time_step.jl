# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

function comp_Δt(fields::StaggeredGridFields, Re::Float64, τ::Float64)
    u_max_abs = maximum(abs.(@view fields.u.values[2:(end - 1), 2:(end - 1)]))
    v_max_abs = maximum(abs.(@view fields.v.values[2:(end - 1), 2:(end - 1)]))
    (; Δx, Δy) = fields.grid_data
    return τ * min(Re / 2 * (1 / Δx^2 + 1 / Δy^2)^-1, Δx / u_max_abs, Δy / v_max_abs)
end
