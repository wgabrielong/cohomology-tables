#=
Pieces shared by the cohomology and the homology paths.

Both reduce to the same problem: a finitely presented module over a PID, given
by a relation matrix, that has to be put in invariant-factor form and named.
Neither side knows about the other -- the cohomology path adds a cup product on
top of this, the homology path adds nothing.

Also holds the wall-clock budget, which both use identically.
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

"""
    _canonical_pieces(R, rel, k) -> (V, Vinv, orders, keep)

[`_smith_form`](@ref) plus the positions that survive: `keep` drops the summands
whose invariant factor is a unit, since `R/(u) = 0`.
"""
function _canonical_pieces(R, rel, k::Int)
  V, Vinv, orders = _smith_form(R, rel, k)
  return V, Vinv, orders, [j for j in 1:k if !is_unit(orders[j])]
end

"""
    _labels(d, k; letter, unit) -> Vector{String}

Names for the `k` canonical basis classes in degree `d`.

`letter` is `"x"` for cohomology and `"a"` for homology, so the two never read
alike. `unit = true` names a lone degree-0 class `"1"`, which is right for a
ring and meaningless for a homology group.
"""
function _labels(d::Int, k::Int; letter::AbstractString = "x", unit::Bool = true)
  unit && d == 0 && k == 1 && return ["1"]
  k == 1 && return ["$letter$d"]
  return ["$(letter)$(d)_$j" for j in 1:k]
end

"Reduce a coordinate modulo the order of its summand; `0` means a free summand."
_reduce(c, m) = is_zero(m) ? c : mod(c, m)
