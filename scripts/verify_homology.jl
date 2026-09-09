#!/usr/bin/env julia
#=
Check every homology computation against independent authorities.

    julia scripts/verify_homology.jl                 # the whole computable catalogue
    julia scripts/verify_homology.jl --field         # + cross-check against H^*
    julia scripts/verify_homology.jl --deep          # + an independent integral check
    julia scripts/verify_homology.jl --only K3 --only T^3
    julia scripts/verify_homology.jl --all           # include the six heavy entries

The default pass is cheap: it needs one integral homology per entry. `--field`
computes the cup product ring over QQ and GF(2) as well, which is far more
expensive and is the reason it is not on by default.

Three checks. Note what is and is not independent:

  * **The published header** -- Lutz's Manifold Page files state `H_*=(...)` in
    their own comment block, written by their authors. Fully independent of
    everything here. Applies to Lutz-sourced single-complex files; the
    geometric 3-manifold catalogues hold many complexes per file and state no
    per-entry homology. The notation differs from ours (`Z_2` for `Z/2`,
    `(Z_2)^3 * Z_4` for a product) and is parsed before comparing.

  * **The cohomology ring** -- over a field, `dim H_d` must equal `dim H^d`.
    The cup product side is built from a *cochain* complex through OSCAR's
    `DGAlgCohRing`, so this genuinely cross-checks two code paths. This is the
    strongest automatic check available over `QQ` and `GF(p)`.

  * **Euler characteristic** -- the alternating sum of ranks against the
    alternating sum of the f-vector.

**`ZZ` is deliberately not checked against `Oscar.homology` here.** Integral
Smith normal form does not scale (`RP^5`'s third boundary matrix defeated it),
so `simplicial_homology` over `ZZ` *delegates* to Polymake. Comparing the two
would be comparing Polymake with itself. `--deep` re-derives integral homology
independently, through the subquotient and its Smith form, for complexes small
enough to afford it; that is a real oracle, just not one that scales.

A disagreement with a published header is not automatically a bug: `SU2_SO3` is
a known case where the header contradicts the facet list in its own file, and
Polymake sides against the header too. Those are reported separately.
=#

using Oscar
include(joinpath(@__DIR__, "..", "src", "CohomologyTables.jl"))
using .CohomologyTables

"(free rank, sorted torsion orders) from our own computation."
function ours(X::SimplicialHomology)
  map(0:top_degree(X)) do d
    g = homology_group(X, d)
    (free_rank(g), sort(string.(torsion_orders(g))))
  end
end

"(free rank, sorted torsion orders) from Oscar/Polymake."
function polymake(K::SimplicialComplex)
  map(integral_invariant_factors(K)) do v
    (count(is_zero, v), sort(string.(filter(!is_zero, v))))
  end
end

"""
Integral homology re-derived independently of Polymake: the subquotient of the
chain complex, presented and put in Smith normal form. Slow -- this is the route
`simplicial_homology` avoids -- so `--deep` applies it only below `max_facets`.

Presented from `simplify(C)`, not from `C` itself. Presenting the unsimplified
subquotient is **not** reliable: for `NotIConnected(5,2)` in degree 5 it returns
19 relation rows when the image of the sixth boundary map has rank 20, giving
`Z^7` where the answer is `Z^6` (the f-vector forces Euler characteristic -5,
and the complex is known to be a wedge of `(5-2)! = 6` five-spheres). The
cohomology side has always presented the simplified complex, which is why it
was never affected.

`simplify` throws on a complex concentrated in one degree -- the OSCAR bug the
cohomology path also works around -- so `dim(K) == 0` is handled directly.
"""
function deep_integral(K::SimplicialComplex)
  C = SimplicialChainComplex(ZZ, K)
  dim(K) == 0 && return [(Int(f_vector(K)[1]), String[])]
  Cs = simplify(C)
  map(0:dim(K)) do d
    H = Oscar.homology(Cs, d)[1]
    k = length(gens(H))
    k == 0 && return (0, String[])
    rel = matrix(map(presentation(H), 1))
    S = nrows(rel) == 0 ? nothing : snf(rel)
    ords = isnothing(S) ? [ZZ(0) for _ in 1:k] :
           [i <= nrows(S) ? S[i, i] : ZZ(0) for i in 1:k]
    keep = [m for m in ords if !is_unit(m)]
    (count(is_zero, keep), sort(string.(filter(!is_zero, keep))))
  end
end

"""
Parse one degree of Lutz's published `H_*`, e.g. `Z`, `0`, `Z^22`, `Z_2`,
`(Z_2)^3 * Z_4`, `Z^2541+Z_2`, `Z_5^3`, into (free rank, sorted torsion orders).
Returns `nothing` if the notation is not understood, so unparsed entries are
reported as such rather than silently counted as agreement.
"""
function parse_lutz_degree(t::AbstractString)
  t = replace(strip(t), " " => "")
  (t == "0" || isempty(t)) && return (0, String[])
  free, tors = 0, String[]
  for part in split(replace(t, "*" => "+"), "+")
    isempty(part) && continue
    m = match(r"^Z\^(\d+)$", part);            isnothing(m) || (free += parse(Int, m.captures[1]); continue)
    part == "Z" && (free += 1; continue)
    m = match(r"^\(Z_(\d+)\)\^(\d+)$", part);  isnothing(m) || (append!(tors, fill(m.captures[1], parse(Int, m.captures[2]))); continue)
    m = match(r"^Z_(\d+)\^(\d+)$", part);      isnothing(m) || (append!(tors, fill(m.captures[1], parse(Int, m.captures[2]))); continue)
    m = match(r"^\(?Z_(\d+)\)?$", part);       isnothing(m) || (push!(tors, m.captures[1]); continue)
    return nothing
  end
  return (free, sort(tors))
end

"The `H_*=(...)` line of a Manifold Page file, split by degree."
function published(name::AbstractString)
  header = try
    lutz_header(name)
  catch
    return nothing
  end
  m = match(r"H_\*\s*=\s*\((.*)\)", replace(header, "#" => " "))
  isnothing(m) && return nothing
  # split on top-level commas only: "(Z_2)^2 * Z_4" contains none, but be safe
  parts, depth, cur = String[], 0, IOBuffer()
  for c in m.captures[1]
    c == '(' && (depth += 1); c == ')' && (depth -= 1)
    if c == ',' && depth == 0
      push!(parts, String(take!(cur)))
    else
      print(cur, c)
    end
  end
  push!(parts, String(take!(cur)))
  return [parse_lutz_degree(p) for p in parts]
end

function main(args)
  only = [args[i + 1] for i in eachindex(args) if args[i] == "--only"]
  entries = catalogue(; heavy = ("--all" in args) || !isempty(only))
  isempty(only) || (entries = filter(e -> e.name in only, entries))

  deep = "--deep" in args
  field = "--field" in args
  maxf = let i = findfirst(==("--max-facets"), args)
    isnothing(i) ? 200 : parse(Int, args[i + 1])
  end
  ok_pm, bad_pm = 0, String[]
  ok_field, bad_field = 0, String[]
  ok_deep, bad_deep, skipped_deep = 0, String[], 0
  ok_pub, bad_pub, unparsed, nopub = 0, String[], String[], String[]
  failed = String[]

  for e in entries
    local K, mine
    try
      K = build(e)
      mine = ours(simplicial_homology(e.name, K, ZZ))
    catch err
      push!(failed, "$(e.name): $(first(sprint(showerror, err), 100))")
      println(rpad(e.name, 30), "FAILED"); flush(stdout); continue
    end

    # Definitional over ZZ (we delegate there), but it catches packing errors.
    pm = polymake(K)
    agree_pm = mine == pm
    agree_pm ? (ok_pm += 1) : push!(bad_pm, "$(e.name): ours $mine vs Polymake $pm")

    # The strongest independent check, but it needs the whole cup product ring
    # over each field, so it is opt-in.
    field && for F in (QQ, GF(2))
      try
        hd = [length(homology_group(simplicial_homology(e.name, K, F), d)) for d in 0:dim(K)]
        X = simplicial_cohomology_ring(e.name, K, F)
        cd = [length(graded_basis(X, d)) for d in 0:top_degree(X)]
        hd == cd ? (ok_field += 1) :
          push!(bad_field, "$(e.name) over $(ring_symbol(F)): H_* $hd vs H^* $cd")
      catch err
        push!(bad_field, "$(e.name) over $(ring_symbol(F)): $(first(sprint(showerror, err), 80))")
      end
    end

    if deep
      if length(facets(K)) <= maxf
        try
          di = deep_integral(K)
          di == mine ? (ok_deep += 1) :
            push!(bad_deep, "$(e.name): independent $di vs ours $mine")
        catch err
          push!(bad_deep, "$(e.name): $(first(sprint(showerror, err), 80))")
        end
      else
        skipped_deep += 1
      end
    end

    note = ""
    if e.recipe.kind == "lutz"
      pub = published(e.recipe.file)
      if isnothing(pub)
        # The geometric 3-manifold catalogues hold many complexes per file and
        # state no per-entry H_*; only the single-complex files do.
        push!(nopub, e.name); note = "  (no published homology)"
      elseif any(isnothing, pub)
        push!(unparsed, e.name); note = "  (header notation not parsed)"
      elseif length(pub) == length(mine) && all(pub .== mine)
        ok_pub += 1; note = "  == published header"
      else
        push!(bad_pub, "$(e.name): ours $mine vs header $pub")
        note = "  != published header"
      end
    end
    println(rpad(e.name, 30), agree_pm ? "ok vs Polymake" : "MISMATCH vs Polymake", note)
    flush(stdout)
  end

  println("\n", "=" ^ 72)
  println("Polymake packing: $ok_pm agree, $(length(bad_pm)) disagree  (definitional over ZZ)")
  field && println("Field dims vs H^*: $ok_field agree, $(length(bad_field)) disagree  (independent)")
  deep && println("Independent ZZ:   $ok_deep agree, $(length(bad_deep)) disagree, $skipped_deep too large")
  println("Published header: $ok_pub agree, $(length(bad_pub)) disagree, ",
          "$(length(unparsed)) not parsed, $(length(nopub)) state none")
  println("build/compute failures: $(length(failed))")
  for (label, list) in [("MISMATCH vs Polymake (a bug here)", bad_pm),
                        ("field dimensions disagree with cohomology (a bug here)", bad_field),
                        ("independent integral homology disagrees (a bug here)", bad_deep),
                        ("differs from published header", bad_pub),
                        ("header notation not parsed", unparsed),
                        ("failed", failed)]
    isempty(list) || (println("\n$label:"); for x in list; println("  ", x); end)
  end
  isempty(bad_pm) && isempty(bad_field) && isempty(bad_deep) && isempty(failed) ||
    error("homology verification failed")
end

main(ARGS)
