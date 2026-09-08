#=
The catalogue: one simplicial model per topological type.

Models come from three places -- a handful built here, SageMath's example
catalogue, and Frank Lutz's Manifold Page (see `sources.jl`). Many spaces are
available from several of them, and many of the library's triangulations are
different models of the same space; the catalogue keeps exactly one entry per
topological type and records the alternatives in the entry's `note`.

Nothing here is stored in the repository. Entries are lazy: `SpaceEntry.build`
is only called when you ask for the complex.
=#

### -------------------------------------------------------------------------
### Models built here
### -------------------------------------------------------------------------

"""
    point_space() -> SimplicialComplex

The one-point space: a single vertex.

This complex is concentrated in cochain degree 0, which trips a bug in OSCAR's
`simplify` for hypercomplexes; `ring_tables.jl` routes 0-dimensional complexes
around it. See the README.
"""
point_space() = simplicial_complex([[1]])

"""
    sphere(n) -> SimplicialComplex

`S^n` as the boundary of the (n+1)-simplex: `n+2` vertices, `n+2` facets, the
vertex-minimal triangulation.
"""
function sphere(n::Int)
  n >= 0 || throw(ArgumentError("sphere(n) needs n >= 0, got $n"))
  return simplicial_complex([[j for j in 1:(n + 2) if j != i] for i in 1:(n + 2)])
end

"""
    kuehnel_real_projective_space(n) -> SimplicialComplex

Kühnel's triangulation of `RP^n` on `2^(n+1) - 1` vertices, for any `n >= 1`.

`S^n` is the boundary of the simplex on `V = {1,...,n+2}`, whose barycentric
subdivision has the proper nonempty subsets of `V` as vertices and the flags
`F_1 < ... < F_{n+1}` as facets. Complementation `F <-> V - F` is a free
simplicial involution on that subdivision and the quotient is `RP^n`; we take
the orbit representative avoiding vertex 1.

Smaller models are known for `n <= 5` and the catalogue prefers them. This is
here for `n >= 6`, where it grows fast: `(n+2)!` permutations are enumerated.
"""
function kuehnel_real_projective_space(n::Int)
  n >= 1 || throw(ArgumentError("kuehnel_real_projective_space(n) needs n >= 1"))
  V = 1:(n + 2)
  rep(F::Vector{Int}) = 1 in F ? sort!(setdiff(collect(V), F)) : F
  ids = Dict{Vector{Int},Int}()
  label(F) = get!(ids, rep(F), length(ids) + 1)
  facets = Set{Vector{Int}}()
  # A flag F_1 < ... < F_{n+1} of proper faces is an ordering of V: F_k is the
  # set of the first k elements.
  for order in _permutations(collect(V))
    push!(facets, sort!([label(sort(order[1:k])) for k in 1:(n + 1)]))
  end
  return simplicial_complex(collect(facets))
end

function _permutations(v::Vector{Int})
  length(v) <= 1 && return [copy(v)]
  out = Vector{Int}[]
  for i in eachindex(v)
    rest = deleteat!(copy(v), i)
    for p in _permutations(rest)
      pushfirst!(p, v[i])
      push!(out, p)
    end
  end
  return out
end

### -------------------------------------------------------------------------
### Catalogue entries
### -------------------------------------------------------------------------

"""
    Recipe

A machine-readable rebuild instruction: enough to reconstruct a model from its
original source without this repository storing any facets.

* `kind = "builtin"` -- `call` is a constructor in this package, e.g. `"sphere(3)"`
* `kind = "sage"`    -- `call` follows `simplicial_complexes.`, e.g. `"K3Surface()"`
* `kind = "lutz"`    -- `file` is a Manifold Page file, `label` picks one complex
                        out of a file holding several (empty if it holds one)
* `kind = "arxiv"`   -- `call` is the arXiv identifier

Recipes are what the generated records store and what
`scripts/fetch_sources.jl` and `scripts/regenerate_sage.jl` replay.
"""
struct Recipe
  kind::String
  call::String
  file::String
  label::String
end

builtin_recipe(call::AbstractString) = Recipe("builtin", call, "", "")
sage_recipe(call::AbstractString) = Recipe("sage", call, "", "")
lutz_recipe(file::AbstractString, label::AbstractString = "") = Recipe("lutz", "", file, label)
arxiv_recipe(id::AbstractString) = Recipe("arxiv", id, "", "")

"""
    is_reproducible(r::Recipe) -> Bool

Whether replaying `r` reproduces the *same labelled* complex, so that a SHA-256
over its canonical form is a meaningful fingerprint.

False for Sage's `SurfaceOfGenus`, which returns a differently labelled
triangulation on each process: same f-vector and homology, different vertex
numbering. Measured over repeated runs -- `SurfaceOfGenus(3)` alternates between
at least two labellings. `SurfaceOfGenus(2)` happened to be stable over three
runs, but it is the same code path, so the whole family is treated as
non-reproducible rather than trusting weak evidence.

For these entries the cohomology ring is still correct -- it is a topological
invariant -- but the cocycle representatives in a record belong to the labelling
that produced them, and `scripts/regenerate_sage.jl` checks isomorphism
invariants instead of the hash.
"""
is_reproducible(r::Recipe) =
  !(r.kind == "sage" && startswith(r.call, "SurfaceOfGenus"))

"""
    build(r::Recipe) -> SimplicialComplex

Execute a recipe. This is the only place a model is actually constructed, so
the record's stated recipe and the complex we compute with cannot drift apart.
"""
function build(r::Recipe)
  r.kind == "sage" && return sage(r.call)
  r.kind == "arxiv" && return bagchi_datta_cp3()
  r.kind == "lutz" && return isempty(r.label) ? lutz(r.file) : lutz(r.file, r.label)
  r.kind == "builtin" && return _build_builtin(r.call)
  throw(ArgumentError("unknown recipe kind \"$(r.kind)\""))
end

function _build_builtin(call::AbstractString)
  call == "point_space()" && return point_space()
  m = match(r"^sphere\((\d+)\)$", call)
  isnothing(m) || return sphere(parse(Int, m.captures[1]))
  m = match(r"^kuehnel_real_projective_space\((\d+)\)$", call)
  isnothing(m) || return kuehnel_real_projective_space(parse(Int, m.captures[1]))
  throw(ArgumentError("unknown built-in constructor \"$call\""))
end

"""
    SpaceEntry

One topological type, one model.

* `name`        -- the catalogue key, e.g. `"K3"`, `"L(3,1)"`, `"RP^4"`
* `type`        -- what the space is, in words
* `provenance`  -- where the model comes from (see [`provenances`](@ref))
* `recipe`      -- how to rebuild it; `build(entry)` executes it
* `coeffs`      -- coefficient rings worth tabulating over
* `note`        -- size, caveats, and other models of the same type
* `heavy`       -- measured not to finish in any reasonable time; the drivers
                   skip these unless asked for by name or with `--all`
"""
struct SpaceEntry
  name::String
  type::String
  provenance::String
  recipe::Recipe
  coeffs::Vector{Any}
  note::String
  heavy::Bool
end

SpaceEntry(name, type, provenance, recipe::Recipe; coeffs = Any[ZZ], note = "", heavy = false) =
  SpaceEntry(name, type, provenance, recipe, collect(Any, coeffs), note, heavy)

"Execute a catalogue entry's recipe."
build(e::SpaceEntry) = build(e.recipe)

is_reproducible(e::SpaceEntry) = is_reproducible(e.recipe)

Base.show(io::IO, e::SpaceEntry) =
  print(io, "SpaceEntry(", repr(e.name), ", ", repr(e.type), ", ", repr(e.provenance), ")")

"The provenance labels used in the catalogue."
provenances() = ["built-in", "Sage", "Lutz/library", "Lutz/further",
                 "Lutz/3-manifolds", "arXiv"]

### --- core spaces ----------------------------------------------------------

function _core_entries()
  E = SpaceEntry[]
  push!(E, SpaceEntry("point", "one-point space", "built-in", builtin_recipe("point_space()");
        note = "every contractible complex has this cohomology ring; other models " *
               "include sage(\"DunceHat()\"), sage(\"RudinBall()\"), " *
               "sage(\"ZieglerBall()\"), sage(\"Simplex(n)\") and Lutz's " *
               "bing, rudin, Abalone, BH_3..BH_5, two_optima, knot"))
  for n in 0:6
    extra = n == 3 ? "; other models: sage(\"BarnetteSphere()\"), " *
                     "sage(\"BrucknerGrunbaumSphere()\"), lutz(\"600_cell\"), " *
                     "lutz(\"trefoil\"), lutz(\"nc_sphere\"), " *
                     "lutz(\"spherical_3manifolds\", \"S^3\")" :
            n == 2 ? "; also sage(\"FareyMap(5)\")" :
            n == 4 ? "; also lutz(\"small_4_manifolds\", \"S^4\")" :
            n == 5 ? "; lutz(\"non_PL_5_sphere\") is a non-PL triangulation of the same space" : ""
    push!(E, SpaceEntry("S^$n", "$n-sphere", "built-in", builtin_recipe("sphere($n)");
          note = "boundary of the $(n+1)-simplex, vertex-minimal" * extra))
  end
  push!(E, SpaceEntry("RP^2", "real projective plane", "Sage", sage_recipe("RealProjectivePlane()"); coeffs = Any[ZZ, GF(2)],
        note = "6 vertices, vertex-minimal; also lutz(\"library\") entries and " *
               "Oscar's real_projective_plane()"))
  push!(E, SpaceEntry("RP^3", "real projective 3-space", "Sage", sage_recipe("RealProjectiveSpace(3)"); coeffs = Any[ZZ, GF(2)],
        note = "11 vertices (Walkup), vertex-minimal; also lutz(\"RP3\"), " *
               "lutz(\"spherical_3manifolds\", \"RP^3\") and lutz(\"Hom_C5_K4\"), " *
               "whose file states Hom(C_5,K_4) = RP^3"))
  push!(E, SpaceEntry("RP^4", "real projective 4-space", "Sage", sage_recipe("RealProjectiveSpace(4)"); coeffs = Any[ZZ, GF(2)],
        note = "16 vertices (Datta), vertex-minimal; also lutz(\"RP4\") and " *
               "lutz(\"small_4_manifolds\", \"RP^4\")"))
  push!(E, SpaceEntry("RP^5", "real projective 5-space", "Lutz/library",
        lutz_recipe("RP5_24"); coeffs = Any[ZZ, GF(2)],
        note = "24 vertices; smaller than Sage's RealProjectiveSpace(5), which " *
               "gives Kühnel's 63-vertex model"))
  push!(E, SpaceEntry("CP^2", "complex projective plane", "Sage", sage_recipe("ComplexProjectivePlane()");
        note = "9 vertices (Kühnel–Banchoff), the unique vertex-minimal model; " *
               "also lutz(\"CP2\") and lutz(\"small_4_manifolds\", \"CP^2\")"))
  push!(E, SpaceEntry("CP^3", "complex projective 3-space", "arXiv", arxiv_recipe("1012.3235");
        note = "18 vertices, 622 facets, Bagchi–Datta arXiv:1012.3235, read from " *
               "the paper's LaTeX source. No explicit triangulation of CP^n is " *
               "known for n >= 4 -- see the README"))
  push!(E, SpaceEntry("HP^2", "quaternionic projective plane", "Sage", sage_recipe("QuaternionicProjectivePlane()");
        note = "15 vertices, 490 facets, dimension 8; the Brehm-Kuehnel " *
               "triangulation with transitive automorphism group, proved to be " *
               "HP^2 by Gorodkov (2019). Also lutz(\"_HP2\")"))
  push!(E, SpaceEntry("Poincare sphere", "Poincaré homology 3-sphere = Sigma(2,3,5)",
        "Sage", sage_recipe("PoincareHomologyThreeSphere()");
        note = "16 vertices; same cohomology ring as S^3. Also lutz(\"poincare\"), " *
               "lutz(\"Poincare_3_sphere\"), lutz(\"spherical_3manifolds\", \"Poincare_sphere\")"))
  push!(E, SpaceEntry("Wu manifold", "SU(3)/SO(3)", "Lutz/library",
        lutz_recipe("SU2_SO3"); coeffs = Any[ZZ, GF(2)],
        note = "13 vertices, vertex-minimal. NB the file's header labels it " *
               "\"SU2/SO3\" and publishes H_* = (Z,0,Z,Z,0,Z), which does not match " *
               "its own facet list -- see the README"))
  push!(E, SpaceEntry("K3", "K3 surface", "Sage", sage_recipe("K3Surface()");
        note = "16 vertices, 288 facets (Casella–Kühnel); H^2 has rank 22, so the " *
               "table is the intersection form. Also lutz(\"K3_16\"), lutz(\"K3_17\")"))
  return E
end

### --- surfaces -------------------------------------------------------------

function _surface_entries()
  E = SpaceEntry[]
  push!(E, SpaceEntry("T^2", "2-torus", "Sage", sage_recipe("Torus()");
        coeffs = Any[ZZ, GF(2)],
        note = "7 vertices, vertex-minimal; exterior algebra on two degree-1 classes"))
  push!(E, SpaceEntry("Klein bottle", "Klein bottle", "Sage", sage_recipe("KleinBottle()");
        coeffs = Any[ZZ, GF(2)], note = "non-orientable genus 2"))
  for g in 2:6
    push!(E, SpaceEntry("Sigma_$g", "orientable surface of genus $g", "Sage", sage_recipe("SurfaceOfGenus($g)"); coeffs = Any[ZZ, GF(2)],
          note = g == 3 ? "also sage(\"FareyMap(7)\")" :
                 g == 6 ? "also sage(\"GenusSix()\") and lutz(\"d2n12g6\")" : ""))
  end
  for k in 3:6
    push!(E, SpaceEntry("N_$k", "non-orientable surface of genus $k", "Sage", sage_recipe("SurfaceOfGenus($k, orientable=False)"); coeffs = Any[ZZ, GF(2)],
          note = "genus 1 is RP^2 and genus 2 is the Klein bottle, both above"))
  end
  # FareyMap(5) and FareyMap(7) are the genus 0 and genus 3 orientable surfaces
  # -- that is, S^2 and Sigma_3, both already here. Only FareyMap(11) reaches a
  # genus this catalogue does not otherwise have.
  push!(E, SpaceEntry("Sigma_26", "orientable surface of genus 26", "Sage",
        sage_recipe("FareyMap(11)"); coeffs = Any[ZZ, GF(2)],
        note = "sage(\"FareyMap(11)\"), the surface of PSL(2, F_11)"))
  push!(E, SpaceEntry("Sigma_15", "orientable surface of genus 15", "Lutz/library",
        lutz_recipe("regular_2_21_23_1"); coeffs = Any[ZZ, GF(2)],
        note = "lutz(\"regular_2_21_23_1\"), 21 vertices, 98 facets"))
  push!(E, SpaceEntry("rand2_n25", "random 2-complex, H_2 = Z^475", "Lutz/library",
        lutz_recipe("rand2_n25_p0.328"),
        note = "25 vertices, 751 facets; H^2 has rank 475, so the table is skipped " *
               "at the default max_table"))
  push!(E, SpaceEntry("PG64", "orientable surface of genus 620", "Lutz/library",
        lutz_recipe("PG64_n2017_o1_g620", "facets"); heavy = true,
        note = "2017 vertices, 6510 facets; H^1 has rank 1240 and it does not finish"))
  push!(E, SpaceEntry("PG128", "non-orientable surface of genus 2542", "Lutz/library",
        lutz_recipe("PG128_PG128P7"); heavy = true,
        note = "127 vertices but H^1 has rank 2541; does not finish"))
  push!(E, SpaceEntry("AG_5_3", "orientable surface of genus 9680", "Lutz/library",
        lutz_recipe("AG_5_3_n29647_o1_g9680", "facets"); heavy = true,
        note = "29647 vertices, 98010 facets; H^1 has rank 19360 and it does not finish"))
  return E
end

### --- Moore spaces ---------------------------------------------------------

"""
Moore spaces `M(Z/q, n)`: one non-trivial reduced homology group, `Z/q` in
degree `n`.

Sage's `MooreSpace(q)` builds `M(Z/q, 1)` only. Higher ones are its iterated
suspension, since `M(Z/q, n) = S^(n-1) M(Z/q, 1)`, which Sage will do for us.
Tabulated over `GF(p)` for the smallest prime `p` dividing `q` -- that is where
the torsion shows; over `ZZ` the cup product vanishes in positive degrees.

`M(Z/2, 1)` is not listed: it is `RP^2`, already in the catalogue.
"""
function _moore_entries()
  spec = [(3, 1, 3), (4, 1, 2), (5, 2, 5), (7, 3, 7), (8, 4, 2), (9, 2, 3)]
  return [SpaceEntry("M(Z/$q,$n)", "Moore space with H_$n = Z/$q", "Sage",
                     sage_recipe(n == 1 ? "MooreSpace($q)" :
                                          "MooreSpace($q).suspension($(n - 1))");
                     coeffs = Any[ZZ, GF(p)],
                     note = n == 1 ? "" :
                            "the $(n - 1)-fold suspension of M(Z/$q,1)")
          for (q, n, p) in spec]
end

### --- 3-manifolds ----------------------------------------------------------

# The labels of Lutz's geometric 3-manifold catalogues, hardcoded so that
# listing the catalogue needs no network. The complexes themselves are fetched
# on demand. Entries duplicating something already in the catalogue are dropped.

const _SPHERICAL_3MANIFOLDS = [
  "L_3_1", "L_4_1", "L_5_1", "L_5_2", "L_6_1", "L_7_1", "L_7_2", "L_8_1",
  "L_8_3", "L_9_1", "L_9_2", "L_10_1", "L_10_3", "cube_space", "P_3", "P_4",
  "P_5", "P_6", "P_7", "P_8", "P_9", "P_10", "octahedron_space",
  "truncated_cube_space"]   # "S^3", "RP^3", "Poincare_sphere" are already above

const _FLAT_3MANIFOLDS = ["T^3", "G2", "G3", "G4", "G5", "G6", "KxS^1", "B2", "B3", "B4"]

const _NIL_3MANIFOLDS = ["Nil_Oo1_1", "Nil_Oo1_2", "Nil_Oo1_3", "Nil_Oo1_4", "Nil_Oo1_5"]

const _S2xR_SPACES = ["S^2twistS^1", "S^2xS^1", "RP^2xS^1", "RP^3#RP^3"]

const _H2xR_SPACES = ["o_g2xS1", "o_g3xS1", "o_g4xS1", "o_g5xS1", "n_g3xS1",
                      "n_g4xS1", "n_g5xS1", "n_g6xS1", "n_g7xS1", "n_g8xS1",
                      "n_g9xS1", "n_g10xS1"]

const _HYPERBOLIC_3MANIFOLDS = [
  "or_0.94270736", "or_0.98136883", "or_1.01494161", "or_1.26370924",
  "or_1.28448530", "or_1.39850888", "or_1.41406104_a", "or_1.41406104_b",
  "or_1.42361190", "or_1.44069901", "or_1.46377664", "or_1.52947733",
  "or_1.54356891_a", "or_1.54356891_b", "or_1.58316666_a", "or_1.58316666_b",
  "or_1.58864664_a", "or_1.58864664_b", "or_1.64960972", "or_1.75712603"]

const _HOMOLOGY_3SPHERES = ["Sigma_2_3_7", "Sigma_2_5_7", "Sigma_3_4_5",
                            "Sigma_3_4_7", "Sigma_3_5_7", "Sigma_4_5_7"]

const _CONNECTED_SUMS_3D = vcat(
  ["(S^2xS^1)#$k" for k in 2:20], ["(S^2twistS^1)#$k" for k in 2:20],
  ["(S^2xS^1)#RP^3", "(S^2twistS^1)#RP^3"],
  ["(S^2xS^1)#$k#RP^3" for k in 2:5], ["(S^2twistS^1)#$k#RP^3" for k in 2:5],
  ["(S^2xS^1)#L_3_1", "(S^2twistS^1)#L_3_1",
   "(S^2xS^1)#2#L_3_1", "(S^2twistS^1)#2#L_3_1",
   "L_3_1#L_3_1", "L_3_1#-L_3_1"])

function _three_manifold_entries()
  E = SpaceEntry[]
  groups = [
    ("spherical_3manifolds", _SPHERICAL_3MANIFOLDS, "spherical 3-manifold"),
    ("flat_3manifolds", _FLAT_3MANIFOLDS, "flat 3-manifold"),
    ("nil_3manifolds", _NIL_3MANIFOLDS, "Nil 3-manifold"),
    ("S2xR_spaces", _S2xR_SPACES, "S^2 x R 3-manifold"),
    ("H2xR_spaces", _H2xR_SPACES, "H^2 x R 3-manifold"),
    ("hyperbolic_3manifolds", _HYPERBOLIC_3MANIFOLDS, "hyperbolic 3-manifold"),
    ("homology_3spheres", _HOMOLOGY_3SPHERES, "Brieskorn homology 3-sphere"),
    ("connected_sums_3d", _CONNECTED_SUMS_3D, "connected sum of 3-manifolds"),
  ]
  for (file, labels, kind) in groups
    for l in labels
      push!(E, SpaceEntry(_pretty_3manifold(l), kind, "Lutz/3-manifolds",
            lutz_recipe(file, l); coeffs = Any[ZZ],
            note = "lutz(\"$file\", \"$l\")"))
    end
  end
  push!(E, SpaceEntry("Weber-Seifert space", "hyperbolic dodecahedral space",
        "Lutz/library", lutz_recipe("hyperbolic_dodecahedral_space");
        coeffs = Any[ZZ, GF(5)], note = "21 vertices; H_1 = (Z/5)^3"))
  return E
end

"`L_3_1` -> `L(3,1)`, `o_g2xS1` -> `Sigma_2 x S^1`, and so on."
function _pretty_3manifold(l::AbstractString)
  m = match(r"^L_(\d+)_(\d+)$", l)
  isnothing(m) || return "L($(m.captures[1]),$(m.captures[2]))"
  m = match(r"^Sigma_(\d+)_(\d+)_(\d+)$", l)
  isnothing(m) || return "Sigma($(m.captures[1]),$(m.captures[2]),$(m.captures[3]))"
  m = match(r"^o_g(\d+)xS1$", l)
  isnothing(m) || return "Sigma_$(m.captures[1]) x S^1"
  m = match(r"^n_g(\d+)xS1$", l)
  isnothing(m) || return "N_$(m.captures[1]) x S^1"
  m = match(r"^or_([\d.]+)(_[ab])?$", l)
  isnothing(m) || return "hyperbolic vol $(m.captures[1])$(something(m.captures[2], ""))"
  return l
end

### --- 4-manifolds and up ---------------------------------------------------

function _high_dimensional_entries()
  E = SpaceEntry[]
  for l in ["S^3xS^1", "S^2xS^2", "CP^2#CP^2", "CP^2#-CP^2", "S^3twistS^1",
            "(S^2xS^2)#(S^2xS^2)", "CP^2#(S^2xS^2)"]
    push!(E, SpaceEntry(l, "4-manifold", "Lutz/further",
          lutz_recipe("small_4_manifolds", l); coeffs = Any[ZZ],
          note = "lutz(\"small_4_manifolds\", \"$l\"); vertex-minimal or close to it"))
  end
  push!(E, SpaceEntry("RP^4#K3", "4-manifold", "Lutz/library",
        lutz_recipe("RP4_K3_17"); coeffs = Any[ZZ, GF(2)], note = "28 vertices, 460 facets"))
  push!(E, SpaceEntry("RP^4#11(S^2xS^2)", "4-manifold", "Lutz/library",
        lutz_recipe("RP4_11S2xS2"); coeffs = Any[ZZ, GF(2)], note = "31 vertices, 518 facets"))
  push!(E, SpaceEntry("S^3xS^2", "5-manifold", "Lutz/further",
        lutz_recipe("S3xS2", "S3xS2_12^a"); coeffs = Any[ZZ],
        note = "12 vertices, vertex-minimal; one of 25 combinatorially distinct " *
               "models in that file -- lutz_labels(\"S3xS2\") lists them"))
  push!(E, SpaceEntry("S^2 x Poincare sphere", "5-manifold", "Lutz/library",
        lutz_recipe("S2xpoincare"); coeffs = Any[ZZ],
        note = "64 vertices, 3600 facets; a few minutes"))
  return E
end

### --- torsion and combinatorial complexes ----------------------------------

function _other_entries()
  E = SpaceEntry[]
  for k in (4, 8, 16, 32)
    push!(E, SpaceEntry("HMT_$k", "Hadamard matrix torsion 2-complex", "Lutz/library",
          lutz_recipe("HMT_$k"); coeffs = Any[ZZ, GF(2)],
          note = "H_1 is a product of cyclic 2-groups of increasing order"))
  end
  # Hom_C5_K4 is not listed: its own file says Hom(C_5,K_4) = RP^3, already here.
  for (nm, kind, coeffs) in [("Hom_C6_compl_K5_small", "Hom complex, 4-dimensional", Any[ZZ]),
                             ("Hom_n9_655_compl_K4", "Hom complex, 3-dimensional", Any[ZZ]),
                             ("Hom_C6_compl_K5", "Hom complex, 4-dimensional", Any[ZZ]),
                             ("Hom_C5_K5", "Hom complex, 5-dimensional", Any[ZZ])]
    big = nm in ("Hom_n9_655_compl_K4", "Hom_C6_compl_K5", "Hom_C5_K5")
    push!(E, SpaceEntry(nm, kind, "Lutz/library", lutz_recipe(nm); coeffs = coeffs,
          heavy = big, note = big ? "19k-95k facets; does not finish" : ""))
  end
  for (n, i) in ((3, 3), (4, 4))
    push!(E, SpaceEntry("Chessboard($n,$i)", "chessboard complex", "Sage", sage_recipe("ChessboardComplex($n,$i)")))
  end
  for n in (5, 6, 7)
    push!(E, SpaceEntry("Matching($n)", "matching complex on $n vertices", "Sage", sage_recipe("MatchingComplex($n)")))
  end
  push!(E, SpaceEntry("NotIConnected(5,2)", "not-2-connected graphs on 5 vertices",
        "Sage", sage_recipe("NotIConnectedGraphs(5,2)")))
  push!(E, SpaceEntry("SumComplex(5,[0,1,3])", "sum complex of Linial–Meshulam–Rosenthal",
        "Sage", sage_recipe("SumComplex(5,[0,1,3])")))
  return E
end

### --- the catalogue --------------------------------------------------------

"""
    catalogue(; heavy = true) -> Vector{SpaceEntry}

Every space in the collection, one entry per topological type. Pass
`heavy = false` to drop the six entries measured not to finish (`PG64`,
`PG128`, `AG_5_3` and the three big Hom complexes).

Building it needs no network and no Sage: entries are lazy, and only
`SpaceEntry.build` (or [`space`](@ref)) reaches out. `catalogue_table` prints
it; `catalogue_space(name)` builds one.
"""
function catalogue(; heavy::Bool = true)
  E = vcat(_core_entries(), _surface_entries(), _moore_entries(),
           _three_manifold_entries(), _high_dimensional_entries(), _other_entries())
  names = [e.name for e in E]
  allunique(names) || throw(ErrorException(
    "duplicate catalogue names: $(join(unique(filter(n -> count(==(n), names) > 1, names)), ", "))"))
  return heavy ? E : filter(e -> !e.heavy, E)
end

"The catalogue keys, in catalogue order."
catalogue_names() = [e.name for e in catalogue()]

"""
    find_space(name) -> SpaceEntry

Look a catalogue entry up by name.
"""
function find_space(name::AbstractString)
  E = catalogue()
  i = findfirst(e -> e.name == name, E)
  isnothing(i) && throw(ArgumentError(
    "no catalogue entry named \"$name\"; see catalogue_names()"))
  return E[i]
end

"""
    catalogue_space(name) -> SimplicialComplex

Build the model for a catalogue entry, e.g. `catalogue_space("K3")`, `catalogue_space("L(3,1)")`.
This is where the network call or the Sage process happens.
"""
catalogue_space(name::AbstractString) = build(find_space(name))

"""
    catalogue_table(io = stdout; provenance = nothing)

Print the catalogue: name, topological type and provenance. Pass `provenance`
to restrict, e.g. `catalogue_table(; provenance = "Sage")`.
"""
function catalogue_table(io::IO = stdout; provenance::Union{Nothing,AbstractString} = nothing)
  E = catalogue()
  isnothing(provenance) || (E = filter(e -> e.provenance == provenance, E))
  w1 = maximum(length(e.name) for e in E)
  w2 = maximum(length(e.type) for e in E)
  println(io, rpad("space", w1), "  ", rpad("topological type", w2), "  provenance")
  println(io, "-" ^ (w1 + w2 + 16))
  for e in E
    println(io, rpad(e.name, w1), "  ", rpad(e.type, w2), "  ", e.provenance)
  end
  println(io, "\n", length(E), " entries")
  return nothing
end
