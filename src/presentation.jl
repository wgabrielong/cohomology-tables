#=
Presentations of the cohomology ring beyond the singly-generated case.

`identify_ring` in report.jl recognises one shape: a monogenic truncated
polynomial algebra `R[g]/(ord*g, g^(m+1))`. That leaves the torus, `T^3`, the
Wu manifold and most 3-manifolds with no presentation at all, even though their
rings are small and completely determined by the table already computed.

The general construction, standard for a graded-commutative algebra:

  1. GENERATORS are a lift of the indecomposables `H^+ / (H^+)^2`, degree by
     degree. The decomposable part of `H^d` is spanned by the products
     `H^i * H^(d-i)` for `0 < i < d`, which `cup` already gives as coordinate
     vectors. Lifts are chosen from the canonical basis itself, so a generator's
     name is one of the `x2_1` labels already printed in the table.

  2. RELATIONS are the kernel of the surjection from the free
     graded-commutative algebra on those generators, computed degree by degree
     and then reduced to those not already implied by lower degrees.

  3. The relation ideal is generated in degrees `<= top + max(|g_i|)`: every
     monomial above `top` is a relation, and one above the bound factors as
     (something already a relation) * (a generator).

Both steps are the same computation over a field and over `ZZ` once the module
structure is carried along: `H^d` is `R^n` modulo the relations `orders[j]*e_j`,
so those rows are stacked onto every span and every kernel. Over a field the
orders are all zero and the extra rows vanish, recovering plain linear algebra.

Squares of odd-degree generators are the one genuinely ring-dependent point.
Graded commutativity gives `2x^2 = 0` for `|x|` odd, so:

  * `QQ`, `GF(p)` with `p` odd -- 2 is invertible, so `x^2 = 0`; the square is
    not even enumerated.
  * `GF(2)` -- graded commutativity is ordinary commutativity and odd squares
    are generally non-zero (`RP^n` is `F_2[x]/(x^(n+1))`).
  * `ZZ` -- `2x^2 = 0` says `x^2` is 2-torsion, *not* that it vanishes. The
    square is enumerated and `x^2 = 0` or `2x^2 = 0` is discovered, not assumed.

Nothing here is trusted: `_presents(X, P)` re-derives the graded pieces from the
presentation and compares them with the computed ring before it is returned.
A presentation that fails is discarded rather than emitted.
=#

"""
    RingGenerator

One algebra generator: a class of the canonical basis of `H^degree`, chosen as a
lift of an indecomposable. `label` is the basis label already used by the table,
`order` its invariant factor (`0` for a free summand).
"""
struct RingGenerator
  degree::Int
  index::Int
  label::String
  order::Any
end

"""
    RingPresentation

`H^*(K;R)` as generators and relations. `text` is the rendered presentation,
`relations` the rendered relation list; both are derived from `generators` and
the exponent vectors the search produced.
"""
struct RingPresentation
  coeff_ring::Any
  generators::Vector{RingGenerator}
  relations::Vector{String}
  text::String
end

Base.length(P::RingPresentation) = length(P.generators)

### -------------------------------------------------------------------------
### Linear algebra over R, carrying the module structure of each graded piece
### -------------------------------------------------------------------------

# The rows `orders[j] * e_j`. Over a field every order is zero and this is empty,
# which is exactly right: H^d is then a plain vector space.
function _ambient_rows(R, orders)
  rows = Vector{Vector{Any}}()
  n = length(orders)
  for j in 1:n
    is_zero(orders[j]) && continue
    push!(rows, Any[i == j ? R(orders[j]) : zero(R) for i in 1:n])
  end
  return rows
end

_as_matrix(R, rows, ncol) =
  isempty(rows) ? zero_matrix(R, 0, ncol) :
                  matrix(R, length(rows), ncol, reduce(vcat, rows))

# Is the row `v` in the R-span of `rows`?
function _in_span(R, rows, v)
  n = length(v)
  n == 0 && return true
  isempty(rows) && return all(is_zero, v)
  A = _as_matrix(R, rows, n)
  return can_solve(A, matrix(R, 1, n, v); side = :left)
end

# Rows spanning {a : a * V lies in the span of `amb`}, as vectors of length
# nrows(V). This is the kernel of the evaluation map into H^d = R^n / amb.
function _relation_kernel(R, V::Vector{Vector{Any}}, amb, n::Int)
  r = length(V)
  r == 0 && return Vector{Vector{Any}}()
  # H^d = 0: every monomial is a relation.
  n == 0 && return [Any[i == k ? one(R) : zero(R) for i in 1:r] for k in 1:r]
  stacked = _as_matrix(R, vcat(V, amb), n)
  K = kernel(stacked; side = :left)
  out = Vector{Vector{Any}}()
  for i in 1:nrows(K)
    row = Any[K[i, j] for j in 1:r]
    all(is_zero, row) || push!(out, row)
  end
  return out
end

### -------------------------------------------------------------------------
### Step 1 and 2: indecomposable generators
### -------------------------------------------------------------------------

"Coordinate vectors of every product `H^i * H^(d-i)`, `0 < i < d`."
function _decomposable_rows(X::CohomologyRing, d::Int)
  rows = Vector{Vector{Any}}()
  for i in 1:(d - 1)
    j = d - i
    bi, bj = graded_basis(X, i), graded_basis(X, j)
    (isempty(bi.elems) || isempty(bj.elems)) && continue
    for p in eachindex(bi.elems), q in eachindex(bj.elems)
      push!(rows, Any[c for c in cup(X, (i, p), (j, q))])
    end
  end
  return rows
end

"""
    indecomposable_generators(X; max_generators) -> Vector{RingGenerator} or nothing

A generating set for `H^*(K;R)` as an `R`-algebra, taken from the canonical
bases: in each positive degree, the basis classes that are not already in the
span of the decomposables and of the classes chosen before them.

`max_generators` stops the search as soon as that many have been found, and
returns `nothing`. This is a real economy, not just a cap on the answer: the
decomposables in degree `d` cost one `cup` per pair from `H^i x H^(d-i)`, so
abandoning a hopeless space in low degree avoids the quadratic work in high
degree. `Hom_C6_compl_K5_small` has 58 classes in `H^2`, and computing `H^2*H^2`
for degree 4 is 3364 cup products -- 306 seconds, all of it wasted on a space
that was never going to be presentable.
"""
function indecomposable_generators(X::CohomologyRing;
                                   max_generators::Int = typemax(Int))
  R = X.coeff_ring
  gens = RingGenerator[]
  for d in 1:top_degree(X)
    b = graded_basis(X, d)
    n = length(b)
    n == 0 && continue
    check_deadline(X.deadline, "$(X.name): generators in degree $d")
    span = vcat(_decomposable_rows(X, d), _ambient_rows(R, b.orders))
    for j in 1:n
      e = Any[i == j ? one(R) : zero(R) for i in 1:n]
      _in_span(R, span, e) && continue
      length(gens) < max_generators || return nothing
      push!(gens, RingGenerator(d, j, b.labels[j], b.orders[j]))
      push!(span, e)
    end
  end
  return gens
end

### -------------------------------------------------------------------------
### Step 3: monomials of the free graded-commutative algebra
### -------------------------------------------------------------------------

# `x^2 = 0` for odd |x| exactly when 2 is invertible (see the file header).
_odd_squares_possible(R) = !is_unit(R(2))

"""
Exponent vectors of every monomial of total degree `1:bound`, with the exponent
of an odd-degree generator capped at 1 unless the coefficient ring allows an odd
square. Returned grouped by degree: `out[d]` is the list for degree `d`.
"""
function _monomials_by_degree(degs::Vector{Int}, bound::Int, odd_squares::Bool)
  k = length(degs)
  caps = [(!odd_squares && isodd(degs[i])) ? 1 : bound ÷ max(degs[i], 1) for i in 1:k]
  out = [Vector{Vector{Int}}() for _ in 1:bound]
  e = zeros(Int, k)
  function rec(i::Int, used::Int)
    if i > k
      used >= 1 && push!(out[used], copy(e))
      return
    end
    for c in 0:caps[i]
      step = c * degs[i]
      used + step > bound && break
      e[i] = c
      rec(i + 1, used + step)
    end
    e[i] = 0
  end
  rec(1, 0)
  return out
end

### -------------------------------------------------------------------------
### Step 4: evaluating a monomial in the ring
### -------------------------------------------------------------------------

function _evaluate(X::CohomologyRing, gens::Vector{RingGenerator}, e::Vector{Int})
  d = sum(e[i] * gens[i].degree for i in eachindex(e); init = 0)
  d > top_degree(X) && return Any[]
  acc = one(X.A)
  for i in eachindex(e), _ in 1:e[i]
    acc *= graded_basis(X, gens[i].degree).elems[gens[i].index]
  end
  return Any[c for c in class_coords(X, d, acc)]
end

### -------------------------------------------------------------------------
### Step 5 and 6: relations, and dropping those implied by lower degrees
### -------------------------------------------------------------------------

# Monomials are stored as exponent vectors and always evaluated in generator
# order, so a product of two of them has to be reordered to that canonical form
# -- and every time two odd-degree generators swap, the sign flips. Without this
# the implied relations below are wrong: `Sigma_2 x S^1`, whose five generators
# are all in degree 1, came out with a degree-3 piece of dimension 0 instead
# of 1. (Over GF(2) the sign is vacuous, but the formula is the same.)
function _koszul_sign(R, degs::Vector{Int}, a::Vector{Int}, b::Vector{Int})
  swaps = 0
  for i in eachindex(a), j in 1:(i - 1)
    swaps += a[i] * b[j] * degs[i] * degs[j]
  end
  return isodd(swaps) ? -one(R) : one(R)
end

# Relations already forced in degree d by relations of lower degree: each is
# some lower relation multiplied by a monomial, re-expressed in degree d's
# monomial basis.
function _implied_rows(R, mons_by_deg, found, d::Int, index_of, degs::Vector{Int})
  rows = Vector{Vector{Any}}()
  r = length(mons_by_deg[d])
  for (d0, rel, mons0) in found
    d0 < d || continue
    for m in mons_by_deg[d - d0]
      row = Any[zero(R) for _ in 1:r]
      nonzero = false
      for (t, c) in enumerate(rel)
        is_zero(c) && continue
        pos = get(index_of[d], mons0[t] .+ m, 0)
        pos == 0 && continue
        row[pos] += c * _koszul_sign(R, degs, mons0[t], m)
        nonzero = true
      end
      nonzero && any(!is_zero, row) && push!(rows, row)
    end
  end
  return rows
end

"""
    ring_presentation(X; max_generators, max_relations) -> RingPresentation or nothing

`H^*(K;R)` as generators and relations, or `nothing` when there is none worth
printing -- too many generators, too many relations, or a presentation that
failed its own verification.
"""
function ring_presentation(X::CohomologyRing; max_generators::Int = 12,
                           max_relations::Int = 80, max_classes::Int = 500)
  R = X.coeff_ring
  length(graded_basis(X, 0)) == 1 || return nothing
  # A coarse pre-guard, cheap to evaluate; the real economy is the early bail
  # inside the generator search below.
  sum(d -> length(graded_basis(X, d)), 1:top_degree(X); init = 0) <= max_classes ||
    return nothing
  gens = indecomposable_generators(X; max_generators = max_generators)
  (isnothing(gens) || isempty(gens)) && return nothing

  degs = [g.degree for g in gens]
  top = top_degree(X)
  bound = top + maximum(degs)
  mons_by_deg = _monomials_by_degree(degs, bound, _odd_squares_possible(R))
  index_of = [Dict(m => t for (t, m) in enumerate(mons_by_deg[d])) for d in 1:bound]

  found = Tuple{Int,Vector{Any},Vector{Vector{Int}}}[]   # (degree, coeffs, monomials)
  for d in 1:bound
    mons = mons_by_deg[d]
    isempty(mons) && continue
    check_deadline(X.deadline, "$(X.name): relations in degree $d")
    n = d <= top ? length(graded_basis(X, d)) : 0
    amb = d <= top ? _ambient_rows(R, graded_basis(X, d).orders) : Vector{Vector{Any}}()
    V = [_evaluate(X, gens, m) for m in mons]
    kern = _relation_kernel(R, V, amb, n)
    isempty(kern) && continue
    implied = _implied_rows(R, mons_by_deg, found, d, index_of, degs)
    for row in kern
      _in_span(R, implied, row) && continue
      push!(found, (d, row, mons))
      push!(implied, row)
      length(found) <= max_relations || return nothing
    end
  end

  _presents(X, gens, found, mons_by_deg) || return nothing
  # Sort so the presentation is a function of the ring, not of the order the
  # kernel happened to come back in: records must be byte-reproducible.
  rels = sort!([(d, _render(gens, rel, mons)) for (d, rel, mons) in found])
  relations = [r for (_, r) in rels]
  return RingPresentation(R, gens, relations, _render_text(R, gens, relations))
end

### -------------------------------------------------------------------------
### Verification: does the presentation actually present the ring?
### -------------------------------------------------------------------------

# The graded piece of the presented algebra in degree d is
# (monomials of degree d) / (relations of degree d, including those implied).
# Its invariant factors must match those of H^d.
function _presents(X::CohomologyRing, gens, found, mons_by_deg)
  R = X.coeff_ring
  degs = [g.degree for g in gens]
  index_of = [Dict(m => t for (t, m) in enumerate(mons_by_deg[d]))
              for d in 1:length(mons_by_deg)]
  for d in 1:top_degree(X)
    mons = mons_by_deg[d]
    r = length(mons)
    b = graded_basis(X, d)
    rows = _implied_rows(R, mons_by_deg, found, d, index_of, degs)
    for (d0, rel, _) in found
      d0 == d && push!(rows, rel)
    end
    # rank/torsion of R^r / rows must equal that of H^d
    _same_group(R, rows, r, b.orders) || return false
    # and the evaluation map must be onto: every basis class is hit
    V = [_evaluate(X, gens, m) for m in mons]
    amb = _ambient_rows(R, b.orders)
    for j in eachindex(b.elems)
      e = Any[i == j ? one(R) : zero(R) for i in eachindex(b.elems)]
      _in_span(R, vcat(V, amb), e) || return false
    end
  end
  return true
end

# Do R^r / rows and (the group with invariant factors `orders`) agree?
function _same_group(R, rows, r::Int, orders)
  want = sort!([string(m) for m in orders if !is_unit(m)])
  A = _as_matrix(R, rows, r)
  if R === ZZ
    got = String[]
    if r > 0
      S = nrows(A) == 0 ? zero_matrix(ZZ, 0, r) : snf(A)
      diag = [i <= nrows(S) ? S[i, i] : ZZ(0) for i in 1:r]
      got = sort!([string(m) for m in diag if !is_unit(m)])
    end
    return got == want
  end
  # over a field only the dimension can differ
  dim = r - (nrows(A) == 0 ? 0 : rank(A))
  return dim == count(m -> !is_unit(m), orders)
end

### -------------------------------------------------------------------------
### Rendering
### -------------------------------------------------------------------------

function _render_monomial(gens, e::Vector{Int})
  parts = String[]
  for i in eachindex(e)
    e[i] == 0 && continue
    push!(parts, e[i] == 1 ? gens[i].label : "$(gens[i].label)^$(e[i])")
  end
  return isempty(parts) ? "1" : join(parts, "*")
end

function _render(gens, rel::Vector{Any}, mons)
  # A relation and its negative say the same thing; print the one with a
  # positive leading coefficient so `2*x3` does not come out as `-2*x3`.
  # (`is_negative` is meaningless over a finite field, where there is no order;
  # there every coefficient already prints as a non-negative residue.)
  lead = findfirst(!is_zero, rel)
  if !isnothing(lead) && applicable(is_negative, rel[lead]) && is_negative(rel[lead])
    rel = Any[-c for c in rel]
  end
  terms = String[]
  for (t, c) in enumerate(rel)
    is_zero(c) && continue
    m = _render_monomial(gens, mons[t])
    if is_one(c)
      push!(terms, m)
    elseif is_one(-c)
      push!(terms, "-$m")
    else
      push!(terms, "$(c)*$m")
    end
  end
  s = replace(join(terms, " + "), "+ -" => "- ")
  # If it still reads with a leading minus, the negated relation says the same
  # thing and reads better. (Over a finite field there is no order to appeal to,
  # so this is the only normalisation available there.)
  if startswith(s, "-")
    neg = _render_terms(gens, Any[-c for c in rel], mons)
    startswith(neg, "-") || return neg
  end
  return s
end

function _render_terms(gens, rel::Vector{Any}, mons)
  terms = String[]
  for (t, c) in enumerate(rel)
    is_zero(c) && continue
    m = _render_monomial(gens, mons[t])
    if is_one(c);      push!(terms, m)
    elseif is_one(-c); push!(terms, "-$m")
    else               push!(terms, "$(c)*$m")
    end
  end
  return replace(join(terms, " + "), "+ -" => "- ")
end

function _render_text(R, gens, relations)
  sym = ring_symbol(R)
  vars = join([g.label for g in gens], ",")
  degs = join(["|$(g.label)| = $(g.degree)" for g in gens], ", ")
  # An exterior algebra deserves its own notation, and writing it that way makes
  # the answer agree across coefficient rings: over QQ the squares of odd
  # generators are never enumerated (2 is invertible, so x^2 = 0 holds in the
  # free object), while over ZZ they turn up as explicit relations. Both are the
  # same algebra and should print as such.
  if all(g -> isodd(g.degree), gens) &&
     sort(relations) == sort(["$(g.label)^2" for g in gens if !is_zero(g.order)] ∪
                             ["$(g.label)^2" for g in gens if is_zero(g.order)])
    return "Lambda_$sym($(join([g.label for g in gens], ", "))),  $degs"
  end
  all(g -> isodd(g.degree), gens) && isempty(relations) &&
    return "Lambda_$sym($(join([g.label for g in gens], ", "))),  $degs"
  isempty(relations) && return "$sym[$vars] (free graded-commutative),  $degs"
  return "$sym[$vars]/($(join(relations, ", "))),  $degs"
end
