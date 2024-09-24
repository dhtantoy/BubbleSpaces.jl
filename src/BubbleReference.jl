struct Bubble <: ReferenceFEName end

function Gridap.ReferenceFEs.return_cache(b::BubbleMonomialBasis, xs::AbstractVector{<:Point{2}})
    cache = return_cache(b.monomials, xs)
    return cache
end

function _lincomb(vals::Matrix{T1}, b::BubbleMonomialBasis{T2}, n_p) where {T1, T2}
    n = num_components(T2)
    vec = b.vec
    @check length(vec) * n == size(vals, 2) "`length(vec) * n` ≠ `size(vals, 2)`"
    @check size(vals, 1) == n_p "`size(vals, 1)` ≠ `n_p`"
    
    r = zeros(T1, n_p, n)
    @inbounds for i = axes(r, 1)
        for j = eachindex(vec)
            for l = axes(r, 2)
                r[i, l] += vec[j] * vals[i, n*(j-1) + l]
            end
        end
    end
    return r
end

function Gridap.evaluate!(cache, b::BubbleMonomialBasis, xs::AbstractVector{<:Point{2}}) 
    vals = evaluate!(cache, b.monomials, xs)
    r = _lincomb(vals, b, length(xs))

    return r
end
function Gridap.ReferenceFEs.return_cache(
    fg::FieldGradientArray{1,<:BubbleMonomialBasis},
    xs::AbstractVector{<:Point{2}})
    cfg = FieldGradientArray{1}(fg.fa.monomials)
    return return_cache(cfg, xs)
end
function Gridap.evaluate!(
    cache,
    fg::FieldGradientArray{1,<:BubbleMonomialBasis},
    xs::AbstractVector{<:Point})
    b = fg.fa
    cfg = FieldGradientArray{1}(b.monomials)
    val = evaluate!(cache, cfg, xs)
    r = _lincomb(val, b, length(xs))

    return r
end

function BubbleRefFE(T::Type, p::Polytope)
    l = LagrangianRefFE(T, p, 0).reffe
    dofs = l.dofs 
    conformity = l.conformity
    face_dofs = l.face_dofs
    ndofs = l.ndofs
    prebasis = BubbleMonomialBasis(T)
    metadata = nothing 
    shapefuncs = compute_shapefuns(dofs, prebasis)

    return GenericRefFE{Bubble}(
        ndofs,
        p,
        prebasis,
        dofs,
        conformity,
        metadata,
        face_dofs,
        shapefuncs
    )
end