struct AdditionBasis{T1, T2}
    dv1::T1
    dv2::T2
end
Gridap.:+(dv1::MultiField.MultiFieldFEBasisComponent, dv2::MultiField.MultiFieldFEBasisComponent) = AdditionBasis(dv1, dv2)

Gridap.dot(f, a::AdditionBasis) = dot(f, a.dv1) + dot(f, a.dv2)
Gridap.dot(a::AdditionBasis, f) = dot(f, a)
for op in (:inner, :dot)
    @eval Gridap.$op(a::AdditionBasis, b::AdditionBasis) = $op(a.dv1, b.dv1) + $op(a.dv2, b.dv2) + $op(a.dv1, b.dv2) + $op(a.dv2, b.dv1)
end
for op in (:gradient,)
    @eval Gridap.$op(a::AdditionBasis) = $op(a.dv1) + $op(a.dv2)
end