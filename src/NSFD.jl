# SPDX-License-Identifier: MPL-2.0
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

module NSFD

export N, S, E, W, C, CellPos

abstract type CellDir end
abstract type XDir <: CellDir end
abstract type YDir <: CellDir end

struct N <: YDir end
struct S <: YDir end
struct E <: XDir end
struct W <: XDir end
struct C <: CellDir end

struct CellPos{X<:Union{XDir,C},Y<:Union{YDir,C}} end

function Base.show(io::IO, ::Type{CellPos{X,Y}}) where {X,Y}
    return print(io, X, ",", Y)
end

export StaggeredValue
include("staggered_value.jl")

export StaggeredField, set!
include("staggered_field.jl")

end
