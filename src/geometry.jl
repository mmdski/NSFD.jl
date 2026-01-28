# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct Geometry
    x_length::Float64
    y_length::Float64
    imax::Int
    jmax::Int
    δx::Float64
    δy::Float64
    function Geometry(x_length::Float64, y_length::Float64, imax::Int, jmax::Int)
        return new(x_length, y_length, imax, jmax, x_length / imax, y_length / jmax)
    end
end
