module BubbleSpaces

# Write your package code here.
using Gridap
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.Fields
using Gridap.Helpers
using Gridap.MultiField

export BubbleRefFE

include("BubbleMonomialBase.jl")
include("BubbleReference.jl")
include("AdditionBasis.jl")

end
