# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

struct FG
    f::StaggeredField{EC}
    g::StaggeredField{CN}
    function FG(imax::Int, jmax::Int)
        return new(StaggeredField(EC, imax, jmax), StaggeredField(CN, imax, jmax))
    end
end
