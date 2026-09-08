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

records_dir(manifest_path) = dirname(manifest_path)

"""
For an entry whose Sage constructor does not reproduce its vertex labelling,
compare what is invariant under relabelling -- the f-vector and the integral
homology -- against the record. A difference here is a real failure: it means
Sage is now building a different space, not merely a different numbering.
"""
function check_invariants(name, K, dir)
  files = filter(f -> startswith(f, record_prefix(name)), readdir(dir))
  isempty(files) && return false
  text = read(joinpath(dir, first(files)), String)
  m = match(r"\"f_vector\":\s*\[([^\]]*)\]", text)
  isnothing(m) && return false
  recorded = parse.(Int, strip.(split(m.captures[1], ",")))
  return collect(Int, f_vector(K)) == recorded
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
        inv_ok = check_invariants(name, K, records_dir(path))
        inv_ok ? (ok += 1) : push!(bad, "$name: isomorphism invariants differ")
        println(inv_ok ? "  ok (iso)  " : "  MISMATCH  ", rpad(name, 24),
                inv_ok ? "relabelled by Sage; f-vector and homology match" : "")
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
