#!/usr/bin/env julia
#=
Rebuild every Sage-sourced complex named in records/MANIFEST.tsv and verify it
against the hash recorded there.

    julia scripts/regenerate_sage.jl
    julia scripts/regenerate_sage.jl --manifest records/MANIFEST.tsv

Sage's example definitions and vertex labellings can change between releases, so
the recorded version string is load-bearing. A version difference is reported;
a hash mismatch is a hard failure, not a warning -- the record was computed from
a different complex than the one your Sage just built.

Sage library code is GPL-2.0-or-later. This repository invokes the constructors
and does not copy any facet list out of Sage; see docs/PROVENANCE.md.
=#

using Oscar
include(joinpath(@__DIR__, "..", "src", "CohomologyTables.jl"))
using .CohomologyTables

function read_manifest(path)
  lines = readlines(path)
  isempty(lines) && error("empty manifest: $path")
  header = split(lines[1], '\t')
  return [Dict(zip(header, split(l, '\t'))) for l in lines[2:end] if !isempty(strip(l))]
end

"""
For an entry whose Sage constructor does not reproduce its vertex labelling,
compare what *is* invariant under relabelling -- the f-vector and the integral
homology -- against the manifest. A difference here is a real failure: it means
Sage is now building a different space, not merely a different numbering.

Read from the manifest rather than by scanning `records/`: a space has more than
one record (a ring record per coefficient ring, and a homology record per ring),
so picking one by filename prefix was ambiguous.
"""
function check_invariants(K, row)
  want_f = get(row, "f_vector", "")
  isempty(want_f) && return (false, "no recorded f-vector")
  join(collect(Int, f_vector(K)), ",") == want_f || return (false, "f-vector differs")
  want_h = get(row, "integral_homology", "")
  isempty(want_h) && return (true, "f-vector matches; no recorded homology")
  join(integral_homology_symbols(K), ",") == want_h ||
    return (false, "integral homology differs")
  return (true, "f-vector and homology match")
end

function main(args)
  i = findfirst(==("--manifest"), args)
  path = isnothing(i) ? joinpath(@__DIR__, "..", "records", "MANIFEST.tsv") : args[i + 1]
  isfile(path) || error("no manifest at $path -- run scripts/export_records.jl first")
  sage_available() || error("no working `sage` on the PATH (set SAGE to override)")
  here = sage_version()
  rows = filter(r -> r["kind"] == "sage", read_manifest(path))
  println("verifying $(length(rows)) Sage models with: $here")
  ok, bad, drift = 0, String[], String[]
  for r in rows
    name, want = r["space"], r["facets_sha256"]
    r["sage_version"] == here || push!(drift, "$name recorded under $(r["sage_version"])")
    try
      K = sage(r["call"])
      got = facets_sha256(K)
      if got == want
        ok += 1
        println("  ok        ", rpad(name, 24), "simplicial_complexes.", r["call"])
      elseif get(r, "reproducible", "true") == "false"
        # This constructor does not reproduce its labelling; check what is
        # actually invariant instead, and still fail if that differs.
        inv_ok, why = check_invariants(K, r)
        inv_ok ? (ok += 1) : push!(bad, "$name: $why")
        println(inv_ok ? "  ok (iso)  " : "  MISMATCH  ", rpad(name, 24),
                inv_ok ? "relabelled by Sage; $why" : why)
      else
        push!(bad, "$name (simplicial_complexes.$(r["call"])): expected $want, got $got")
        println("  MISMATCH  ", name)
      end
    catch err
      push!(bad, "$name: $(sprint(showerror, err))")
      println("  FAILED    ", name)
    end
    flush(stdout)
  end
  println("\n$ok verified, $(length(bad)) failed")
  if !isempty(drift)
    println("$(length(drift)) recorded under a different Sage version:")
    for d in drift
      println("  ", d)
    end
  end
  if !isempty(bad)
    for b in bad
      println("  ", b)
    end
    error("Sage regeneration failed -- the recorded complexes are not what this " *
          "Sage builds; do not trust the affected records")
  end
end

main(ARGS)
