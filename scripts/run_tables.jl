#!/usr/bin/env julia
#=
Print cohomology rings and cup product multiplication tables for the catalogue.

    julia scripts/run_tables.jl --list                 # just show the catalogue
    julia scripts/run_tables.jl                        # compute everything
    julia scripts/run_tables.jl --only K3 --only "L(3,1)"
    julia scripts/run_tables.jl --provenance Sage      # one source only
    julia scripts/run_tables.jl --coeffs GF7           # override coefficients
    julia scripts/run_tables.jl --time-limit 120 --out tables.txt
    julia scripts/run_tables.jl --no-homology            # cup products only

Nothing is stored in the repo: models are built by Sage or downloaded from the
Manifold Page on first use and cached outside it (see `cache_dir()`).

The --time-limit budget is cooperative -- one long call into OSCAR can overrun
it. scripts/run_isolated.jl enforces the same budget as a hard wall-clock kill.
=#

using Oscar
include(joinpath(@__DIR__, "..", "src", "CohomologyTables.jl"))
using .CohomologyTables

function parse_coeffs(s)
  s == "ZZ" && return Any[ZZ]
  s == "QQ" && return Any[QQ]
  m = match(r"^GF(\d+)$", s)
  isnothing(m) || return Any[GF(parse(Int, m.captures[1]))]
  error("unrecognised --coeffs $s (try ZZ, QQ, GF2, GF7)")
end

function run_all(io::IO, entries; time_limit = DEFAULT_TIME_LIMIT, coeffs = nothing,
                 homology::Bool = true)
  for e in entries
    for R in (isnothing(coeffs) ? e.coeffs : coeffs)
      t0 = time()
      # Homology first and in its own try: it is much cheaper than the cup
      # product, so a ring that times out still leaves the groups printed.
      if homology
        try
          H = simplicial_homology(e, R; time_limit = time_limit)
          print_homology(io, H; note = "$(e.type) -- $(e.provenance)")
        catch err
          println(io, "$(e.name) homology over $(ring_symbol(R)) SKIPPED: ",
                  sprint(showerror, err))
          println(io)
        end
      end
      try
        X = simplicial_cohomology_ring(e, R; time_limit = time_limit)
        print_report(io, X; note = "$(e.type) -- $(e.provenance)" *
                                   (isempty(e.note) ? "" : "; $(e.note)"))
        println(io, "(computed in $(round(time() - t0; digits = 2)) s)")
      catch err
        println(io, "=" ^ 72)
        println(io, "$(e.name)   coefficients in $(ring_symbol(R))")
        println(io, "=" ^ 72)
        println(io, "SKIPPED: ", sprint(showerror, err))
      end
      println(io)
      flush(io)
    end
  end
end

function main(args)
  if "--list" in args
    p = findfirst(==("--provenance"), args)
    catalogue_table(; provenance = isnothing(p) ? nothing : args[p + 1])
    return
  end
  only = [args[i + 1] for i in eachindex(args) if args[i] == "--only"]
  # `--all` includes the six entries measured not to finish; naming one
  # explicitly with --only always works.
  entries = catalogue(; heavy = ("--all" in args) || !isempty(only))
  isempty(only) || (entries = filter(e -> e.name in only, entries))
  p = findfirst(==("--provenance"), args)
  isnothing(p) || (entries = filter(e -> e.provenance == args[p + 1], entries))
  c = findfirst(==("--coeffs"), args)
  coeffs = isnothing(c) ? nothing : parse_coeffs(args[c + 1])
  hom = !("--no-homology" in args)
  j = findfirst(==("--time-limit"), args)
  tl = isnothing(j) ? DEFAULT_TIME_LIMIT : parse(Float64, args[j + 1])
  i = findfirst(==("--out"), args)
  if isnothing(i)
    run_all(stdout, entries; time_limit = tl, coeffs = coeffs, homology = hom)
  else
    open(args[i + 1], "w") do io
      run_all(io, entries; time_limit = tl, coeffs = coeffs, homology = hom)
    end
    @info "wrote $(args[i + 1])"
  end
end

main(ARGS)
