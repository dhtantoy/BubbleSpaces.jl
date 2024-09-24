struct AdditionBasis{T1, T2}
    dv1::T1
    dv2::T2
end
Gridap.:+(dv1::MultiField.MultiFieldFEBasisComponent, dv2::MultiField.MultiFieldFEBasisComponent) = AdditionBasis(dv1, dv2)


for op in (:inner, :dot, :*)
    @eval Gridap.$op(a::AdditionBasis, f) = $op(f, a)
    @eval Gridap.$op(f, b::AdditionBasis) = $op(f, b.dv1) + $op(f, b.dv2) 
    @eval Gridap.$op(a::AdditionBasis, b::AdditionBasis) = $op(a.dv1, b) + $op(a.dv2, b)
end
for op in (:gradient,)
    @eval Gridap.$op(a::AdditionBasis) = $op(a.dv1) + $op(a.dv2)
end