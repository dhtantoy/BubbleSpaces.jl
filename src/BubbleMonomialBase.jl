struct BubbleMonomialBasis{T, eT} <: AbstractVector{Monomial}
    monomials::MonomialBasis{2, T}
    vec::Vector{eT}
    function BubbleMonomialBasis(T::Type) 
        eT = eltype(T)
        terms = [
            CartesianIndex(3, 2), # x²y 
            CartesianIndex(2, 3), # xy²
            CartesianIndex(2, 2), # xy
        ]
        monomials = MonomialBasis{2}(T, (2, 2), terms)
        vec = eltype(T)[-27., -27., 27.]
        new{T, eT}(monomials, vec)
    end
end
Base.size(::BubbleMonomialBasis{T}) where T = (num_components(T),)
Base.getindex(::BubbleMonomialBasis, ::Integer) = Monomial()
Base.IndexStyle(::BubbleMonomialBasis) = IndexLinear()
Polynomials.get_order(::BubbleMonomialBasis) = 0
ReferenceFEs.return_type(::BubbleMonomialBasis{T}) where T = T
