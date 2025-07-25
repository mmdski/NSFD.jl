# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct Geometry
    grid_data::GridData
    fluid_cells::Vector{CartesianIndex{2}}
end

function Geometry(grid_data::GridData)
    fluid_indices = vec([CartesianIndex(i, j)
                         for i in 1:(grid_data.nx), j in 1:(grid_data.ny)])
    return Geometry(grid_data, fluid_indices)
end
