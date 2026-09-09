#!/usr/bin/env julia
#=
Generate the homology records under records/.

    julia scripts/export_homology.jl                      # the whole catalogue
    julia scripts/export_homology.jl --only K3 --only "L(3,1)"
    julia scripts/export_homology.jl --coeffs GF7
    julia scripts/export_homology.jl --time-limit 120 --all

One record per (space, coefficient ring), over ZZ, QQ, GF(2), GF(3), GF(5) and
GF(7) by default. Files are named `<slug>_<digest>.homology.<coeff>.json`,
alongside the cup product ring records, which are left untouched.

Records carry the groups in invariant-factor form and **no cycle
representatives** -- for a closed orientable manifold the top-degree one is
supported on every facet, so emitting it would redistribute the input this
repository deliberately does not ship. See docs/PROVENANCE.md.

This does not touch MANIFEST.tsv: the manifest anchors the *inputs*, and these
records reference exactly the complexes it already lists.
=#

using Oscar
include(joinpath(@__DIR__, "..", "src", "CohomologyTables.jl"))
using .CohomologyTables

const DEFAULT_GRID = Any[ZZ, QQ, GF(2), GF(3), GF(5), GF(7)]

function parse_coeffs(s)
  s == "ZZ" && return Any[ZZ]
  s == "QQ" && return Any[QQ]
  m = match(r"^GF(\d+)$", s)
  isnothing(m) || return Any[GF(parse(Int, m.captures[1]))]
  error("unrecognised --coeffs $s (try ZZ, QQ, GF2, GF7)")
end

function date_for(r::Recipe)
  r.kind == "lutz" && return lutz_cached_at(r.file)
  r.kind == "arxiv" && return arxiv_cached_at(r.call)
  return nothing
end

function main(args)
  outdir = let i = findfirst(==("--out"), args)
    isnothing(i) ? joinpath(@__DIR__, "..", "records") : args[i + 1]
  end
  only = [args[i + 1] for i in eachindex(args) if args[i] == "--only"]
  entries = catalogue(; heavy = ("--all" in args) || !isempty(only))
  isempty(only) || (entries = filter(e -> e.name in only, entries))
  c = findfirst(==("--coeffs"), args)
  coeffs = isnothing(c) ? DEFAULT_GRID : parse_coeffs(args[c + 1])
  j = findfirst(==("--time-limit"), args)
  tl = isnothing(j) ? DEFAULT_TIME_LIMIT : parse(Float64, args[j + 1])

  mkpath(outdir)
  written, skipped = 0, String[]
  for e in entries
    local K
    try
      K = build(e)
    catch err
      push!(skipped, "$(e.name): could not build -- $(sprint(showerror, err))")
      continue
    end
    date = date_for(e.recipe)
    for R in coeffs
      try
        X = simplicial_homology(e.name, K, R; time_limit = tl)
        path = write_homology_record(outdir, e, X; date_accessed = date)
        written += 1
        println("wrote ", basename(path))
      catch err
        push!(skipped, "$(e.name) over $(ring_symbol(R)): $(sprint(showerror, err))")
      end
      flush(stdout)
    end
  end

  println("\n$written homology records in $outdir")
  if !isempty(skipped)
    println("$(length(skipped)) skipped:")
    for s in skipped
      println("  ", s)
    end
  end
end

main(ARGS)
