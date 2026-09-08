#=
Cup product multiplication tables from Oscar's experimental hypercomplexes.

The Oscar side of this (experimental/DoubleAndHyperComplexes) gives us:

  Oscar.SimplicialCochainComplex(R, K)  the simplicial cochain complex C^*(K;R),
                                        built as a 1-dimensional hypercomplex,
                                        with `mul_cochains` implementing the
                                        Alexander-Whitney cup product
                                          (a u b)(v_0..v_{p+q}) =
                                             a(v_0..v_p) * b(v_p..v_{p+q})
  Oscar.DGAlgCohRing(C)                 the cohomology ring H^*(K;R) of a
                                        cochain complex carrying that DG-algebra
                                        structure
  small_generating_set(A, d)            a small generating set of H^d, pulled
                                        back from a simplified (algorithmically
                                        reduced) model of the cochain complex
  `*` on DGAlgCohRingElem               the induced product on cohomology

What it does *not* give is coordinates: a product comes back as a cohomology
class, not as a combination of named generators. That is what this file adds.

For each degree d we

  1. take `small_generating_set(A, d)` as generators g_1..g_k of H^d;
  2. read off the relation matrix `Rel` between them from `presentation`;
  3. put `Rel` in Smith normal form, `T*Rel*V = S`, which turns H^d into
     the canonical form  (+)_j R/(S_jj);
  4. name the resulting canonical basis and express every cup product in it.

Step 3 is what makes the tables canonical rather than "some representative":
over ZZ it separates free summands from torsion (so RP^4 prints H^2 = Z/2, not
just "one generator"), and it drops generators that were redundant.
=#

### -------------------------------------------------------------------------
### Time budget
### -------------------------------------------------------------------------

"Default wall-clock budget, in seconds, for computing one cohomology ring."
const DEFAULT_TIME_LIMIT = 300.0

struct TimeLimitExceeded <: Exception
  limit::Float64
  stage::String
  elapsed::Float64
end

Base.showerror(io::IO, e::TimeLimitExceeded) =
  print(io, "time limit of $(e.limit)s exceeded during ", e.stage,
            " (after ", round(e.elapsed; digits = 1), "s)")

"""
A wall-clock budget, checked cooperatively.

This is a *cooperative* deadline: it is tested between degrees while the graded
pieces are built and between rows while the multiplication table is filled, so
it stops work at the next boundary rather than pre-empting mid-call. A single
Oscar call that runs long on its own -- `small_generating_set` on a complex with
thousands of cohomology generators is the realistic case -- can therefore
overrun the budget before the next check is reached. Julia cannot interrupt a
blocking call into Singular or Flint, so enforcing a hard limit needs process
isolation, which this repo does not do.
"""
mutable struct Deadline
  start::Float64
  limit::Float64
end

Deadline(limit::Real) = Deadline(time(), Float64(limit))
no_deadline() = Deadline(time(), Inf)
elapsed(d::Deadline) = time() - d.start

function check_deadline(d::Deadline, stage::AbstractString)
  e = elapsed(d)
  e > d.limit && throw(TimeLimitExceeded(d.limit, String(stage), e))
  return nothing
end

### -------------------------------------------------------------------------
### Data
### -------------------------------------------------------------------------

"""
Canonical basis of a single cohomology group `H^d`.

* `elems[j]`  -- the basis class, as an element of the `DGAlgCohRing`
* `orders[j]` -- its order: `0` means a free summand `R`, `m` means `R/(m)`
* `labels[j]` -- the name used in printed tables

`full_orders` and `keep` record the pre-pruning picture: the Smith form has one
diagonal entry per raw generator, and the ones that came out units correspond to
summands that vanish, so they are dropped from the basis. `keep` lists the
surviving positions, which is how a raw coordinate row is cut down to a
canonical one.
"""
struct GradedBasis
  degree::Int
  elems::Vector{Any}
  orders::Vector{Any}
  labels::Vector{String}
  to_coords::Any        # H^d class -> raw coordinates w.r.t. small_generating_set
  V::Any                # canonical coords = raw coords * V
  full_orders::Vector{Any}
  keep::Vector{Int}
end

Base.length(b::GradedBasis) = length(b.elems)

"""
Everything computed for one space and one coefficient ring.
"""
struct CohomologyRing
  name::String
  complex::SimplicialComplex
  coeff_ring::Any
  A::Any                       # the Oscar DGAlgCohRing
  bases::Vector{GradedBasis}   # bases[d+1] is the basis of H^d
  deadline::Deadline           # budget carried into table building and checks
end

top_degree(X::CohomologyRing) = length(X.bases) - 1
graded_basis(X::CohomologyRing, d::Int) = X.bases[d + 1]

"All canonical basis classes, as `(degree, index)` pairs, in degree order."
function basis_indices(X::CohomologyRing)
  return [(d, j) for d in 0:top_degree(X) for j in 1:length(graded_basis(X, d))]
end

### -------------------------------------------------------------------------
### Construction
### -------------------------------------------------------------------------

"""
    simplicial_cohomology_ring(name, K, R) -> CohomologyRing

Compute `H^*(K; R)` with its cup product, in canonical form.
"""
function simplicial_cohomology_ring(name::AbstractString, K::SimplicialComplex, R;
                                   normalize::Bool = true,
                                   time_limit::Real = DEFAULT_TIME_LIMIT,
                                   deadline::Deadline = Deadline(time_limit))
  C = Oscar.SimplicialCochainComplex(R, K)
  A = Oscar.DGAlgCohRing(C)
  bases = GradedBasis[]
  for d in 0:dim(K)
    check_deadline(deadline, "$name: building H^$d")
    push!(bases, _build_basis(A, K, R, d))
  end
  X = CohomologyRing(String(name), K, R, A, bases, deadline)
  return normalize ? normalize_to_powers(X) : X
end

simplicial_cohomology_ring(e::SpaceEntry, R; kwargs...) =
  simplicial_cohomology_ring(e.name, e.build(), R; kwargs...)

"""
    normalize_to_powers(X) -> CohomologyRing

If `H^*` is generated by a single class `t` of degree `e`, re-choose the basis
of each `H^(k*e)` to be `t^k`.

Smith normal form picks a generator of each cyclic summand arbitrarily up to a
unit, so without this the table can come out with entries like `x2 * x4 = -x6`
for `CP^3` -- correct, but a distracting way to display `Z[t]/(t^4)`. Only the
choice of generators changes; the ring does not. A no-op when `H^*` is not
singly generated.
"""
function normalize_to_powers(X::CohomologyRing)
  top = top_degree(X)
  length(graded_basis(X, 0)) == 1 || return X
  positive = [d for d in 1:top if !isempty(graded_basis(X, d).elems)]
  isempty(positive) && return X
  e = minimum(positive)
  m = top ÷ e
  for d in 1:top
    length(graded_basis(X, d)) == ((d % e == 0 && d <= m * e) ? 1 : 0) || return X
  end
  # t^k = u_k * (current basis class of H^(k*e)) for some unit u_k.
  units = Vector{Any}(undef, m)
  for k in 1:m
    c = cup_power(X, (e, 1), k)
    (length(c) == 1 && is_unit(c[1])) || return X
    units[k] = c[1]
  end
  bases = copy(X.bases)
  check_deadline(X.deadline, "$(X.name): normalising generators")
  for k in 2:m                          # k = 1 is already t itself
    bases[k * e + 1] = _rescale(X.A, X.coeff_ring, bases[k * e + 1], units[k])
  end
  return CohomologyRing(X.name, X.complex, X.coeff_ring, X.A, bases, X.deadline)
end

"Replace the single basis class `g` of a degree by `u*g`, `u` a unit."
function _rescale(A, R, b::GradedBasis, u)
  uinv = divexact(one(R), u)
  V = deepcopy(b.V)
  j = only(b.keep)
  for i in 1:nrows(V)
    V[i, j] = V[i, j] * uinv        # coordinates rescale inversely
  end
  elems = Any[A(u) * only(b.elems)]
  return GradedBasis(b.degree, elems, b.orders, b.labels, b.to_coords, V,
                     b.full_orders, b.keep)
end

function _build_basis(A, K::SimplicialComplex, R, d::Int)
  gens_d, to_coords, rel = _raw_generators(A, K, d)
  k = length(gens_d)
  V, Vinv, orders = _smith_form(R, rel, k)

  # Canonical basis element j is  sum_i Vinv[j,i] * g_i, keeping only the
  # summands that did not collapse (order a unit => that summand is zero).
  elems = Any[]
  ords = Any[]
  keep = Int[]
  for j in 1:k
    is_unit(orders[j]) && continue
    e = zero(A)
    for i in 1:k
      c = Vinv[j, i]
      is_zero(c) && continue
      e += A(R(c)) * gens_d[i]
    end
    push!(elems, e)
    push!(ords, orders[j])
    push!(keep, j)
  end
  labels = _labels(d, length(elems))
  return GradedBasis(d, elems, ords, labels, to_coords, V, orders, keep)
end

"""
Generators of `H^d` before canonicalisation, a function turning a class into
its coordinate row with respect to them, and the relations between them.

The normal route is Oscar's `small_generating_set`, which goes through the
simplified cochain complex. A complex concentrated in a single degree (i.e.
`dim(K) == 0`, the point) breaks Oscar's `simplify`, so there we read the
homology of the unsimplified complex instead -- it is already minimal, since
`H^0 = C^0` when there are no differentials.
"""
function _raw_generators(A, K::SimplicialComplex, d::Int)
  if dim(K) == 0
    H = Oscar.graded_part(A, d)
    gens_d = Any[Oscar.DGAlgCohRingElem(A, d, v) for v in gens(H)]
    to_coords = x -> coordinates(Oscar.graded_part(x, d))
    return gens_d, to_coords, matrix(map(presentation(H), 1))
  end

  gens_d = Any[g for g in small_generating_set(A, d)]
  s = Oscar.simplified_cochain_complex(A)
  Hs = Oscar.simplified_graded_part(A, d)
  to_s = Oscar.map_from_original_complex(s)[d]
  to_coords = x -> coordinates(Hs(to_s(repres(Oscar.graded_part(x, d))); check = false))
  return gens_d, to_coords, matrix(map(presentation(Hs), 1))
end

"""
Smith normal form of the relation matrix, as `(V, Vinv, orders)`.

`H^d = R^k / rowspan(rel)`.  With `T*rel*V = S` in Smith form and `T`, `V`
unimodular, `c -> c*V` is an isomorphism carrying `rowspan(rel)` onto
`rowspan(S)`, so `H^d = (+)_j R/(S_jj)`.
"""
function _smith_form(R, rel, k::Int)
  Id = identity_matrix(R, k)
  (k == 0 || nrows(rel) == 0) && return Id, Id, Any[zero(R) for _ in 1:k]
  S, T, V = snf_with_transform(rel)
  @assert T * rel * V == S "unexpected snf_with_transform convention"
  orders = Any[i <= nrows(S) ? S[i, i] : zero(R) for i in 1:k]
  return V, inv(V), orders
end

function _labels(d::Int, k::Int)
  d == 0 && k == 1 && return ["1"]
  k == 1 && return ["x$d"]
  return ["x$(d)_$j" for j in 1:k]
end

### -------------------------------------------------------------------------
### Coordinates and products
### -------------------------------------------------------------------------

"""
    class_coords(X, d, a) -> Vector

Canonical coordinates of the class `a` in `H^d`: the vector `c` with
`a = sum_j c[j] * graded_basis(X, d).elems[j]`, each entry reduced modulo the order of
its summand.
"""
function class_coords(X::CohomologyRing, d::Int, a)
  b = graded_basis(X, d)
  isempty(b.elems) && return Any[]
  raw = b.to_coords(a)
  k = ncols(b.V)
  row = matrix(X.coeff_ring, 1, k, [raw[i] for i in 1:k]) * b.V
  return Any[_reduce(row[1, j], b.full_orders[j]) for j in b.keep]
end

_reduce(c, m) = is_zero(m) ? c : mod(c, m)

"""
    cup(X, (p, i), (q, j)) -> Vector

Canonical coordinates in `H^(p+q)` of the cup product of the `i`-th basis class
of `H^p` with the `j`-th basis class of `H^q`. Empty vector when `p + q` is
above the top degree.
"""
function cup(X::CohomologyRing, (p, i)::Tuple{Int,Int}, (q, j)::Tuple{Int,Int})
  p + q > top_degree(X) && return Any[]
  a = graded_basis(X, p).elems[i]
  b = graded_basis(X, q).elems[j]
  return class_coords(X, p + q, a * b)
end

"""
    cup_power(X, (d, i), k) -> Vector

Canonical coordinates of the `k`-th cup power of a basis class.
"""
function cup_power(X::CohomologyRing, (d, i)::Tuple{Int,Int}, k::Int)
  k * d > top_degree(X) && return Any[]
  a = graded_basis(X, d).elems[i]
  acc = one(X.A)
  for _ in 1:k
    acc = acc * a
  end
  return class_coords(X, k * d, acc)
end
