#=
Generated ring records: a reproducible recipe, a hash of the input complex, and
the computed cohomology ring.

A record contains no facet list. It carries enough to rebuild the complex from
its original source (a Sage constructor invocation, or a file on Lutz's Manifold
Page) plus a SHA-256 over the canonicalised facets, so a rebuild that does not
reproduce the input fails loudly instead of silently computing something else.

The one place faces of the input appear is the support of a cocycle
representative, which names the faces a chosen cochain is non-zero on. That is a
small subset of the complex's faces and does not determine it; the facet list
stays with the upstream source.

The canonicalisation rule is fixed and documented in docs/CANONICALIZATION.md;
`CANONICAL_FORM_ID` versions it.
=#

"Identifier for the canonicalisation rule; bump if the rule ever changes."
const CANONICAL_FORM_ID = "cohomology-tables/canonical-facets/1"

"Identifier for the cup product ring record schema."
const RECORD_SCHEMA_ID = "cohomology-tables/ring-record/1"

"Identifier for the homology record schema."
const HOMOLOGY_RECORD_SCHEMA_ID = "cohomology-tables/homology-record/1"

"""
    canonical_facets(K) -> Vector{Vector{Int}}

The canonical facet list of `K`: 0-based integer labels, vertices ascending
within each facet, facets in lexicographic order.

Complexes reach this point already relabelled to `1:n` by
[`relabel_facets`](@ref) or by the Sage bridge, both of which fix the label
order (numeric for integer labels, `str` order for Sage's arbitrary labels).
This subtracts one and sorts. See docs/CANONICALIZATION.md.
"""
function canonical_facets(K::SimplicialComplex)
  fs = [sort!([v - 1 for v in f]) for f in [sort(collect(s)) for s in facets(K)]]
  return sort!(fs; lt = _lex_less)
end

function _lex_less(a::Vector{Int}, b::Vector{Int})
  for i in 1:min(length(a), length(b))
    a[i] != b[i] && return a[i] < b[i]
  end
  return length(a) < length(b)
end

"""
    canonical_form(K) -> String

The canonical facet list rendered for hashing: one facet per line, vertices
comma-separated in decimal, lines joined with `\\n`, no trailing newline.
"""
canonical_form(K::SimplicialComplex) =
  join((join(f, ",") for f in canonical_facets(K)), "\n")

"""
    facets_sha256(K) -> String

Lowercase hex SHA-256 of [`canonical_form`](@ref), encoded UTF-8. This is the
fingerprint stored in a record and checked when a complex is re-fetched or
regenerated.
"""
facets_sha256(K::SimplicialComplex) = bytes2hex(SHA.sha256(canonical_form(K)))

### -------------------------------------------------------------------------
### A very small JSON writer
### -------------------------------------------------------------------------
#
# Hand-rolled so that the package does not gain a dependency purely to emit a
# few hundred flat records. Only the types used below are supported.

# Every method takes `indent` (most ignore it) so that keyword dispatch from the
# container methods reaches the scalar ones; without it they fall through to the
# stringifying fallback and integers come out quoted.
_json(io::IO, ::Nothing; indent::Int = 0) = print(io, "null")
_json(io::IO, x::Bool; indent::Int = 0) = print(io, x ? "true" : "false")
_json(io::IO, x::Integer; indent::Int = 0) = print(io, x)
_json(io::IO, x::AbstractString; indent::Int = 0) = _json_string(io, x)

function _json_string(io::IO, s::AbstractString)
  print(io, '"')
  for c in s
    if c == '"'
      print(io, "\\\"")
    elseif c == '\\'
      print(io, "\\\\")
    elseif c == '\n'
      print(io, "\\n")
    elseif c == '\r'
      print(io, "\\r")
    elseif c == '\t'
      print(io, "\\t")
    elseif c < ' '
      print(io, "\\u", lpad(string(UInt16(c); base = 16), 4, '0'))
    else
      print(io, c)
    end
  end
  print(io, '"')
end

function _json(io::IO, v::AbstractVector; indent::Int = 0)
  isempty(v) && return print(io, "[]")
  print(io, "[")
  for (i, x) in enumerate(v)
    i > 1 && print(io, ", ")
    _json(io, x; indent = indent)
  end
  print(io, "]")
end

function _json(io::IO, d::AbstractDict; indent::Int = 0)
  isempty(d) && return print(io, "{}")
  pad, inner = " " ^ indent, " " ^ (indent + 2)
  print(io, "{\n")
  ks = sort!(collect(keys(d)); by = String)   # deterministic key order
  for (i, k) in enumerate(ks)
    print(io, inner)
    _json_string(io, String(k))
    print(io, ": ")
    _json(io, d[k]; indent = indent + 2)
    i < length(ks) && print(io, ",")
    print(io, "\n")
  end
  print(io, pad, "}")
end

_json(io::IO, x; indent::Int = 0) = _json_string(io, string(x))

"Serialise a record to JSON text."
function json_string(d::AbstractDict)
  io = IOBuffer()
  _json(io, d; indent = 0)
  print(io, "\n")
  return String(take!(io))
end

### -------------------------------------------------------------------------
### Building a record
### -------------------------------------------------------------------------

"""
    ring_record(entry, X; date_accessed) -> Dict

Assemble the record for one catalogue entry and one computed cohomology ring:
the rebuild recipe, the input fingerprint, the graded pieces with a cocycle
representative for each basis class, and the structure constants of the cup
product.

`date_accessed` is only meaningful for entries fetched from the network; pass
the date the fetch actually happened. It is `nothing` by default rather than
being filled in with today's date, since that would be a guess about when the
cache was populated.
"""
function ring_record(e::SpaceEntry, X::CohomologyRing;
                     date_accessed::Union{Nothing,AbstractString} = nothing,
                     max_basis::Int = 40)
  K = X.complex
  rec = Dict{String,Any}()
  rec["schema"] = RECORD_SCHEMA_ID
  rec["space"] = e.name
  rec["topological_type"] = e.type
  rec["provenance"] = e.provenance
  rec["source"] = _source_block(e, date_accessed)
  rec["complex"] = _complex_block(K)
  rec["coefficients"] = ring_symbol(X.coeff_ring)
  rec["vertex_order_convention"] =
    "Alexander-Whitney cup product with respect to the ascending order of the " *
    "canonical 0-based vertex labels; see docs/CANONICALIZATION.md"
  rec["cohomology"] = [_degree_block(X, d) for d in 0:top_degree(X)]
  n_basis = length(basis_indices(X))
  rec["structure_constants"] = _structure_constants(X; max_basis = max_basis)
  rec["structure_constants_complete"] = n_basis <= max_basis
  n_basis <= max_basis || (rec["structure_constants_omitted_because"] =
    "$n_basis basis classes exceeds max_basis = $max_basis; the full list would " *
    "have $(n_basis^2) entries. Compute individual products with cup(X, (p,i), (q,j)).")
  rec["ring_presentation"] = identify_ring(X)
  rec["reproducible_labelling"] = is_reproducible(e.recipe)
  rec["consistency_checks"] = check_ring(X)
  return rec
end

"The input fingerprint, shared by both record kinds."
_complex_block(K::SimplicialComplex) = Dict{String,Any}(
  "n_vertices" => n_vertices(K),
  "n_facets" => length(facets(K)),
  "dimension" => dim(K),
  "f_vector" => collect(Int, f_vector(K)),
  "canonicalization" => CANONICAL_FORM_ID,
  "facets_sha256" => facets_sha256(K),
)

function _source_block(e::SpaceEntry, date_accessed)
  r = e.recipe
  b = Dict{String,Any}("kind" => r.kind)
  if r.kind == "sage"
    b["invocation"] = "simplicial_complexes.$(r.call)"
    b["sage_version"] = sage_version()
    b["note"] = "regenerate with scripts/regenerate_sage.jl; Sage example " *
                "definitions and vertex labellings can change between releases, " *
                "so the version pin is load-bearing"
    b["reproducible_labelling"] = is_reproducible(e.recipe)
    is_reproducible(e.recipe) || (b["labelling_caveat"] =
      "simplicial_complexes.$(r.call) returns a differently labelled complex on " *
      "each Sage process. The facets_sha256 below is the labelling this record " *
      "was computed from and will not generally be reproduced; the cohomology " *
      "ring is unaffected, being a topological invariant, but the cocycle " *
      "representatives belong to that labelling. regenerate_sage.jl checks " *
      "isomorphism invariants for these entries instead of the hash.")
  elseif r.kind == "lutz"
    b["catalog"] = r.file
    b["label"] = isempty(r.label) ? nothing : r.label
    b["url"] = lutz_url(r.file)
    b["date_accessed"] = date_accessed
    b["note"] = "fetch with scripts/fetch_sources.jl; not redistributed here"
  elseif r.kind == "arxiv"
    b["arxiv_id"] = r.call
    b["url"] = "https://arxiv.org/abs/$(r.call)"
    b["date_accessed"] = date_accessed
  elseif r.kind == "builtin"
    b["invocation"] = r.call
    b["note"] = "constructed by this package; see src/spaces.jl"
  end
  return b
end

function _degree_block(X::CohomologyRing, d::Int)
  b = graded_basis(X, d)
  return Dict{String,Any}(
    "degree" => d,
    "group" => group_symbol(X.coeff_ring, b.orders),
    "invariant_factors" => [string(m) for m in b.orders],
    "basis_labels" => b.labels,
    "cocycle_representatives" =>
      [_cocycle(X, d, j) for j in eachindex(b.elems)],
  )
end

"""
A representative cocycle for one basis class, as a sparse combination of dual
simplices: each term gives the face (canonical 0-based vertex labels) and the
coefficient. An empty list means the class is represented by the zero cochain.
"""
function _cocycle(X::CohomologyRing, d::Int, j::Int)
  v = Oscar.graded_part(graded_basis(X, d).elems[j], d)
  chain = repres(v)
  fs = [sort!([w - 1 for w in collect(s)]) for s in faces(X.complex, d)]
  terms = Dict{String,Any}[]
  for (idx, coeff) in coordinates(chain)
    is_zero(coeff) && continue
    push!(terms, Dict{String,Any}("face" => fs[idx], "coefficient" => string(coeff)))
  end
  return terms
end

"""
Structure constants of the cup product: for each ordered pair of basis classes
whose degrees sum to at most the top degree, the coordinates of their product in
the basis of the target degree. Pairs with a zero product are omitted.

The list is quadratic in the total number of basis classes, so it is only
emitted when that number is at most `max_basis`. `rand2_n25` has `H^2` of rank
475, which would be some 227,000 entries; above the cap the record still carries
every graded piece and cocycle representative, and says the constants were
omitted. Compute them yourself with `cup(X, (p,i), (q,j))`.
"""
function _structure_constants(X::CohomologyRing; max_basis::Int = 40)
  out = Dict{String,Any}[]
  idx = basis_indices(X)
  length(idx) > max_basis && return out
  for a in idx, b in idx
    a[1] + b[1] > top_degree(X) && continue
    c = cup(X, a, b)
    all(is_zero, c) && continue
    push!(out, Dict{String,Any}(
      "a" => Dict{String,Any}("degree" => a[1], "index" => a[2],
                              "label" => class_label(X, a)),
      "b" => Dict{String,Any}("degree" => b[1], "index" => b[2],
                              "label" => class_label(X, b)),
      "product_degree" => a[1] + b[1],
      "coefficients" => [string(x) for x in c],
    ))
  end
  return out
end

### -------------------------------------------------------------------------
### Writing records out
### -------------------------------------------------------------------------

"""
    record_filename(name, R) -> String

A filesystem-safe name for the record of one space over one coefficient ring.
Space names contain `#`, `^`, `(` and `)`, so everything outside `[A-Za-z0-9]`
becomes `_`; a short digest of the original name is appended so that distinct
spaces cannot collide after that flattening.
"""
function record_filename(name::AbstractString, R; kind::AbstractString = "")
  coeff = replace(ring_symbol(R), r"[^A-Za-z0-9]+" => "")
  return "$(record_prefix(name; kind = kind))$(coeff).json"
end

"""
    record_prefix(name) -> String

The filename prefix, up to and including the dot, shared by every record for one
space: `"<slug>_<digest>."`.

Matching records by slug alone is wrong -- `Sigma_4` is a prefix of
`Sigma_4_x_S_1` -- so the eight-hex digest of the untruncated name is part of the
prefix and makes it unambiguous.

`kind` selects a record family: `""` is the cup product ring record (so existing
filenames are unchanged) and `"homology"` gives `<slug>_<digest>.homology.<coeff>.json`.
`record_prefix(name)` with no kind still matches every record for a space.
"""
function record_prefix(name::AbstractString; kind::AbstractString = "")
  slug = strip(replace(name, r"[^A-Za-z0-9]+" => "_"), '_')
  isempty(slug) && (slug = "space")
  base = "$(slug)_$(first(bytes2hex(SHA.sha256(name)), 8))."
  return isempty(kind) ? base : "$(base)$(kind)."
end

"""
    write_record(dir, entry, X; date_accessed) -> String

Write the record for `entry`/`X` into `dir` and return the path.
"""
function write_record(dir::AbstractString, e::SpaceEntry, X::CohomologyRing;
                      date_accessed::Union{Nothing,AbstractString} = nothing,
                      max_basis::Int = 40)
  isdir(dir) || mkpath(dir)
  path = joinpath(dir, record_filename(e.name, X.coeff_ring))
  open(path, "w") do io
    print(io, json_string(ring_record(e, X; date_accessed = date_accessed,
                                      max_basis = max_basis)))
  end
  return path
end

### -------------------------------------------------------------------------
### Homology records
### -------------------------------------------------------------------------

"""
    homology_record(entry, X::SimplicialHomology; date_accessed) -> Dict

The record for one space and one coefficient ring on the homology side: the
rebuild recipe, the input fingerprint, and the groups in invariant-factor form.

Deliberately absent, and not available behind any flag:

* **cycle representatives.** For a closed orientable `n`-manifold the
  fundamental class is supported on every `n`-facet with an orientation sign, so
  a top-degree representative *is* the facet list -- measured `S^2` 4/4, `T^2`
  14/14, `RP^3` 40/40, `CP^2` 36/36. Emitting them would redistribute the input
  this repository does not ship. Use `cycle_representative(X, d, j)` locally
  instead. See docs/PROVENANCE.md.
* **`vertex_order_convention`.** Homology groups do not depend on the vertex
  order, and nothing order-dependent is emitted, so the field would be noise.
  The canonical form still governs `facets_sha256`.
"""
function homology_record(e::SpaceEntry, X::SimplicialHomology;
                         date_accessed::Union{Nothing,AbstractString} = nothing)
  K = X.complex
  rec = Dict{String,Any}()
  rec["schema"] = HOMOLOGY_RECORD_SCHEMA_ID
  rec["space"] = e.name
  rec["topological_type"] = e.type
  rec["provenance"] = e.provenance
  rec["source"] = _source_block(e, date_accessed)
  rec["complex"] = _complex_block(K)
  rec["coefficients"] = ring_symbol(X.coeff_ring)
  rec["homology"] = [_homology_degree_block(X, d) for d in 0:top_degree(X)]
  rec["ranks"] = homology_ranks(X)
  rec["euler_characteristic"] =
    sum((-1)^d * homology_ranks(X)[d + 1] for d in 0:top_degree(X); init = 0)
  rec["reproducible_labelling"] = is_reproducible(e.recipe)
  rec["consistency_checks"] = check_homology(X)
  rec["note"] = "No cycle representatives are stored: for a closed orientable " *
                "manifold the top-degree one is supported on every facet. See " *
                "docs/PROVENANCE.md."
  return rec
end

function _homology_degree_block(X::SimplicialHomology, d::Int)
  g = homology_group(X, d)
  return Dict{String,Any}(
    "degree" => d,
    "group" => group_symbol(X.coeff_ring, g.orders),
    "invariant_factors" => [string(m) for m in g.orders],
    "basis_labels" => g.labels,
    "rank" => free_rank(g),
    "torsion" => [string(m) for m in torsion_orders(g)],
  )
end

"""
    write_homology_record(dir, entry, X; date_accessed) -> String

Write the homology record for `entry`/`X` into `dir` and return the path.
"""
function write_homology_record(dir::AbstractString, e::SpaceEntry, X::SimplicialHomology;
                               date_accessed::Union{Nothing,AbstractString} = nothing)
  isdir(dir) || mkpath(dir)
  path = joinpath(dir, record_filename(e.name, X.coeff_ring; kind = "homology"))
  open(path, "w") do io
    print(io, json_string(homology_record(e, X; date_accessed = date_accessed)))
  end
  return path
end
