#=
Homology H_*(K; R), in canonical invariant-factor form.

The shape mirrors ring_tables.jl deliberately -- same `Deadline`, same
`_canonical_pieces`, same accessor names -- but there is no ring here. Homology
has no product, so this file is strictly smaller: groups, ranks, torsion, and
nothing else.

Two coefficient regimes, because they have genuinely different costs:

  * **Over a field** (`QQ`, `GF(p)`) there is no torsion, so the groups are
    determined by ranks alone:

        dim H_d = n_d - rank(d_d) - rank(d_{d+1})

    Rank over a field is ordinary linear algebra and stays fast: `RP^5` takes
    about 7 seconds, `CP^3` about 10. This is the case OSCAR cannot do at all,
    and the reason this file exists.

  * **Over `ZZ`** the torsion needs Smith normal form, which suffers the usual
    integer coefficient explosion -- `snf` of `RP^5`'s third boundary matrix
    (2277 x 1174) did not finish in four minutes here. Polymake already has a
    tuned integral routine, so `ZZ` delegates to `Oscar.homology(K, d)` rather
    than recomputing it badly. That call takes no coefficient ring, which is
    exactly why the field case above has to be done by hand.

An earlier version of this file went through `presentation` of the subquotient
`Oscar.homology(C, d)[1]` in every ring. That is correct -- and it is what
`cycle_representative` still uses -- but it does not scale: `RP^5` did not
finish. Note also that the subquotient's `gens` are not a minimal generating
set (over `QQ` the 2-sphere reports 4, 3, 1 generators in degrees 0, 1, 2), so
that route additionally needs Smith normal form to say anything.
=#

"""
Canonical form of a single homology group `H_d`.

* `orders[j]` -- `0` for a free summand `R`, `m` for `R/(m)`
* `labels[j]` -- `a0`, `a1`, `a2_1`, ... (the cohomology side uses `x`)

No representative is stored, here or in any record; see docs/PROVENANCE.md and
[`cycle_representative`](@ref).
"""
struct HomologyGroup
  degree::Int
  orders::Vector{Any}
  labels::Vector{String}
end

Base.length(g::HomologyGroup) = length(g.labels)

"Number of free summands of a homology group, i.e. its Betti number."
free_rank(g::HomologyGroup) = count(is_zero, g.orders)

"The finite invariant factors of a homology group; empty when it is free."
torsion_orders(g::HomologyGroup) = Any[m for m in g.orders if !is_zero(m)]

"""
Everything computed for one space and one coefficient ring, on the homology side.
"""
struct SimplicialHomology
  name::String
  complex::SimplicialComplex
  coeff_ring::Any
  C::Any                          # the SimplicialChainComplex
  groups::Vector{HomologyGroup}   # groups[d+1] is H_d
  deadline::Deadline
end

top_degree(X::SimplicialHomology) = length(X.groups) - 1

"The canonical form of `H_d`."
homology_group(X::SimplicialHomology, d::Int) = X.groups[d + 1]

"""
    homology_ranks(X) -> Vector{Int}

The Betti numbers, degree by degree. Named this rather than `betti_numbers`,
which Oscar already exports.
"""
homology_ranks(X::SimplicialHomology) =
  [free_rank(homology_group(X, d)) for d in 0:top_degree(X)]

### -------------------------------------------------------------------------
### Construction
### -------------------------------------------------------------------------

"""
    simplicial_homology(name, K, R) -> SimplicialHomology

Compute `H_*(K; R)` in canonical invariant-factor form, for `R` any of `ZZ`,
`QQ` or `GF(p)`.

This goes through a chain complex rather than asking OSCAR directly because
`Oscar.homology(K, i)` on a `SimplicialComplex` takes no coefficient ring: it is
Polymake, and integral only. Over `ZZ` the two agree, and the test suite pins
that; over `QQ` and `GF(p)` there is nothing to compare against upstream.

Note this is genuinely different information from the cohomology ring, not a
restatement of it. Over `ZZ` for `RP^4`, `H_* = Z, Z/2, 0, Z/2, 0` while
`H^* = Z, 0, Z/2, 0, Z/2` -- the universal-coefficients degree shift.
"""
function simplicial_homology(name::AbstractString, K::SimplicialComplex, R;
                             time_limit::Real = DEFAULT_TIME_LIMIT,
                             deadline::Deadline = Deadline(time_limit))
  C = SimplicialChainComplex(R, K)
  n = dim(K)
  groups = HomologyGroup[]
  if R === ZZ
    for d in 0:n
      check_deadline(deadline, "$name: building H_$d")
      push!(groups, _integral_group(K, d))
    end
  else
    # One rank per boundary map, reused by the two degrees that need it.
    ranks = zeros(Int, n + 2)                       # ranks[i+1] is rank(d_i)
    for i in 1:n
      check_deadline(deadline, "$name: rank of the degree-$i boundary map")
      ranks[i + 1] = rank(matrix(map(C, 1, (i,))))
    end
    for d in 0:n
      free = Int(f_vector(K)[d + 1]) - ranks[d + 1] - ranks[d + 2]
      push!(groups, _group_from_orders(d, Any[zero(R) for _ in 1:free]))
    end
  end
  return SimplicialHomology(String(name), K, R, C, groups, deadline)
end

_group_from_orders(d::Int, orders::Vector{Any}) =
  HomologyGroup(d, orders, _labels(d, length(orders); letter = "a", unit = false))

"""
`H_d(K; Z)` from Polymake, packed into our canonical order: free summands first,
then torsion. Delegated rather than recomputed because integral Smith normal
form on these boundary matrices is far slower than Polymake's tuned routine.
"""
function _integral_group(K::SimplicialComplex, d::Int)
  v = elementary_divisors(Oscar.homology(K, d))
  orders = vcat(Any[ZZ(0) for _ in 1:count(is_zero, v)],
                Any[m for m in v if !is_zero(m)])
  return _group_from_orders(d, orders)
end

simplicial_homology(e::SpaceEntry, R; kwargs...) =
  simplicial_homology(e.name, build(e), R; kwargs...)

### -------------------------------------------------------------------------
### Representatives -- interactive only, never written to a record
### -------------------------------------------------------------------------

"""
    cycle_representative(X, d, j) -> Vector{Pair{Vector{Int},Any}}

A cycle representing the `j`-th canonical class of `H_d`, as `face => coefficient`
pairs in canonical 0-based vertex labels.

**Deliberately absent from every generated record.** For a closed orientable
`n`-manifold the fundamental class is supported on every `n`-facet with an
orientation sign, so a top-degree cycle representative *is* the facet list:
measured, `S^2` 4 of 4 facets, `T^2` 14 of 14, `RP^3` 40 of 40, `CP^2` 36 of 36.
Serialising that would redistribute the input this repository deliberately does
not ship. Cohomology is safe by the dual accident -- its cocycle supports are
small. See docs/PROVENANCE.md.

Computing one locally is fine, which is what this is for.

Recomputed from the chain complex on demand, via the subquotient
`Oscar.homology(C, d)[1]` and its presentation. That is much slower than the
group computation itself -- it is the route `simplicial_homology` deliberately
avoids -- so expect it to be usable on small complexes only. `j` indexes the
canonical decomposition of that subquotient, which agrees with the recorded
groups as a multiset of invariant factors but need not list them in the same
order.
"""
function cycle_representative(X::SimplicialHomology, d::Int, j::Int)
  H = Oscar.homology(X.C, d)[1]
  k = length(gens(H))
  rel = k == 0 ? zero_matrix(X.coeff_ring, 0, 0) : matrix(map(presentation(H), 1))
  _, Vinv, _, keep = _canonical_pieces(X.coeff_ring, rel, k)
  1 <= j <= length(keep) || throw(BoundsError(keep, j))
  gs = gens(H)
  v = sum((Vinv[keep[j], i] * gs[i] for i in eachindex(gs)); init = zero(H))
  fs = [sort!([w - 1 for w in collect(s)]) for s in faces(X.complex, d)]
  out = Pair{Vector{Int},Any}[]
  for (idx, coeff) in coordinates(repres(v))
    is_zero(coeff) || push!(out, fs[idx] => coeff)
  end
  return out
end

### -------------------------------------------------------------------------
### The independent integral oracle
### -------------------------------------------------------------------------

"""
    integral_invariant_factors(K) -> Vector{Vector{ZZRingElem}}

`H_*(K; Z)` straight from `Oscar.homology`, i.e. Polymake -- a code path with
nothing in common with the chain complex here, which is what makes it usable as
an oracle in the test suite and in [`check_homology`](@ref).
"""
integral_invariant_factors(K::SimplicialComplex) =
  [elementary_divisors(Oscar.homology(K, i)) for i in 0:dim(K)]

"""
    integral_homology_symbols(K) -> Vector{String}

[`integral_invariant_factors`](@ref) rendered as `["Z", "0", "Z/2"]`.
"""
integral_homology_symbols(K::SimplicialComplex) =
  [group_symbol(ZZ, v) for v in integral_invariant_factors(K)]
