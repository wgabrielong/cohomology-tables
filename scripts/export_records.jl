#!/usr/bin/env julia
#=
Generate the ring records under records/.

    julia scripts/export_records.jl                       # everything computable
    julia scripts/export_records.jl --only K3 --only S^3
    julia scripts/export_records.jl --coeffs GF7 --out records
    julia scripts/export_records.jl --manifest-only    # rebuild MANIFEST.tsv alone
    julia scripts/export_records.jl --time-limit 120

Each record holds a rebuild recipe, a SHA-256 over the canonicalised facets, the
graded pieces with a cocycle representative per basis class, and the structure
constants of the cup product. It holds no facets: see docs/CANONICALIZATION.md
and docs/PROVENANCE.md.

Alongside the records this writes records/MANIFEST.tsv, the table that
scripts/fetch_sources.jl and scripts/regenerate_sage.jl verify against.
=#

using Oscar
include(joinpath(@__DIR__, "..", "src", "CohomologyTables.jl"))
using .CohomologyTables

const MANIFEST_HEADER = ["space", "kind", "call", "file", "label",
                         "sage_version", "url", "date_accessed", "reproducible",
                         "facets_sha256"]

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

function manifest_row(e::SpaceEntry, K, date)
  r = e.recipe
  url = r.kind == "lutz" ? lutz_url(r.file) :
        r.kind == "arxiv" ? "https://arxiv.org/abs/$(r.call)" : ""
  return [e.name, r.kind, r.call, r.file, r.label,
          r.kind == "sage" ? sage_version() : "",
          url, something(date, ""), string(is_reproducible(e)), facets_sha256(K)]
end

function main(args)
  outdir = let i = findfirst(==("--out"), args)
    isnothing(i) ? joinpath(@__DIR__, "..", "records") : args[i + 1]
  end
  only = [args[i + 1] for i in eachindex(args) if args[i] == "--only"]
  entries = catalogue(; heavy = ("--all" in args) || !isempty(only))
  isempty(only) || (entries = filter(e -> e.name in only, entries))
  c = findfirst(==("--coeffs"), args)
  coeffs = isnothing(c) ? nothing : parse_coeffs(args[c + 1])
  j = findfirst(==("--time-limit"), args)
  tl = isnothing(j) ? DEFAULT_TIME_LIMIT : parse(Float64, args[j + 1])

  mkpath(outdir)
  rows = Vector{String}[]
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
    push!(rows, manifest_row(e, K, date))
    "--manifest-only" in args && continue
    for R in (isnothing(coeffs) ? e.coeffs : coeffs)
      try
        X = simplicial_cohomology_ring(e.name, K, R; time_limit = tl)
        path = write_record(outdir, e, X; date_accessed = date)
        written += 1
        println("wrote ", basename(path))
      catch err
        push!(skipped, "$(e.name) over $(ring_symbol(R)): $(sprint(showerror, err))")
      end
      flush(stdout)
    end
  end

  # Merge rather than overwrite: a `--only` run must not drop the rows for every
  # other space, since the verification scripts read this file as the whole
  # picture.
  path = joinpath(outdir, "MANIFEST.tsv")
  merged = Dict{String,Vector{String}}()
  if isfile(path)
    old = readlines(path)
    if !isempty(old) && split(old[1], '\t') == MANIFEST_HEADER
      for l in old[2:end]
        isempty(strip(l)) && continue
        cols = String.(split(l, '\t'))
        merged[cols[1]] = cols
      end
    end
  end
  for r in rows
    merged[r[1]] = r
  end
  open(path, "w") do io
    println(io, join(MANIFEST_HEADER, "\t"))
    for k in sort!(collect(keys(merged)))
      println(io, join(merged[k], "\t"))
    end
  end

  println("\n$written records + MANIFEST.tsv in $outdir")
  if !isempty(skipped)
    println("$(length(skipped)) skipped:")
    for s in skipped
      println("  ", s)
    end
  end
end

main(ARGS)
