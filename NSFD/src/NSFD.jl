# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

module NSFD

export GridData
export PField, UField, VField, interior_shape, grid_data, meshgrid

include("grid.jl")
include("field.jl")

end
