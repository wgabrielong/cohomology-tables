#!/usr/bin/env julia
#=
Fetch the Manifold Page sources named in records/MANIFEST.tsv and verify each
against the hash recorded there.

    julia scripts/fetch_sources.jl
    julia scripts/fetch_sources.jl --manifest records/MANIFEST.tsv

Downloads land in the gitignored cache (`cache_dir()`), never in the repo. A
hash mismatch is a hard failure: it means the upstream file changed, or the
recipe no longer selects the complex the record was computed from, and either
way the record must not be trusted until that is understood.

This repository redistributes none of these complexes; see docs/PROVENANCE.md.
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

function main(args)
  i = findfirst(==("--manifest"), args)
  path = isnothing(i) ? joinpath(@__DIR__, "..", "records", "MANIFEST.tsv") : args[i + 1]
  isfile(path) || error("no manifest at $path -- run scripts/export_records.jl first")
  rows = filter(r -> r["kind"] in ("lutz", "arxiv"), read_manifest(path))
  println("verifying $(length(rows)) fetched sources against $path")
  ok, bad = 0, String[]
  for r in rows
    name = r["space"]
    try
      K = build(Recipe(r["kind"], r["call"], r["file"], r["label"]))
      got = facets_sha256(K)
      if got == r["facets_sha256"]
        ok += 1
        println("  ok        ", name)
      else
        push!(bad, "$name: expected $(r["facets_sha256"]), got $got")
        println("  MISMATCH  ", name)
      end
    catch err
      push!(bad, "$name: $(sprint(showerror, err))")
      println("  FAILED    ", name)
    end
    flush(stdout)
  end
  println("\n$ok verified, $(length(bad)) failed")
  if !isempty(bad)
    for b in bad
      println("  ", b)
    end
    error("source verification failed -- do not trust the affected records")
  end
end

main(ARGS)
