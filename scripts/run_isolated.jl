#!/usr/bin/env julia
#=
Hard wall-clock gate: run each space in its own process and kill it at the limit.

    julia scripts/run_isolated.jl                       # 300s per space
    julia scripts/run_isolated.jl --time-limit 120
    julia scripts/run_isolated.jl --provenance Sage --out tables.txt

The `--time-limit` inside `run_tables.jl` is cooperative -- it is only tested
between degrees and between table rows, so one long call into OSCAR (and from
there into Singular or Flint, which Julia cannot interrupt) can blow past it.
This driver gives up on interrupting and kills the process instead, which is the
only way to actually bound the wall clock. The cost is one Julia + OSCAR startup
per space, roughly 35 seconds, so prefer `run_tables.jl` unless you are sweeping
spaces whose cost you do not already know.
=#

using Oscar
include(joinpath(@__DIR__, "..", "src", "CohomologyTables.jl"))
using .CohomologyTables

const WORKER = joinpath(@__DIR__, "run_tables.jl")

function run_isolated(io::IO, entries; time_limit::Real = DEFAULT_TIME_LIMIT, coeffs = nothing)
  for e in entries
    println(io, "-" ^ 72)
    println(io, "### $(e.name)  (hard limit $(time_limit)s)")
    flush(io)
    args = String[WORKER, "--only", e.name, "--time-limit", string(time_limit)]
    isnothing(coeffs) || append!(args, ["--coeffs", coeffs])
    t0 = time()
    proc = run(pipeline(`$(Base.julia_cmd()) $args`; stdout = io, stderr = io); wait = false)
    killer = Timer(Float64(time_limit)) do _
      process_running(proc) && kill(proc, Base.SIGKILL)
    end
    wait(proc)
    close(killer)
    success(proc) || println(io, "KILLED after $(round(time() - t0; digits = 1))s ",
                                 "(hard limit $(time_limit)s)")
    flush(io)
  end
end

function main(args)
  only = [args[i + 1] for i in eachindex(args) if args[i] == "--only"]
  # `--all` includes the six entries measured not to finish; naming one
  # explicitly with --only always works.
  entries = catalogue(; heavy = ("--all" in args) || !isempty(only))
  isempty(only) || (entries = filter(e -> e.name in only, entries))
  p = findfirst(==("--provenance"), args)
  isnothing(p) || (entries = filter(e -> e.provenance == args[p + 1], entries))
  c = findfirst(==("--coeffs"), args)
  coeffs = isnothing(c) ? nothing : args[c + 1]
  j = findfirst(==("--time-limit"), args)
  tl = isnothing(j) ? DEFAULT_TIME_LIMIT : parse(Float64, args[j + 1])
  i = findfirst(==("--out"), args)
  if isnothing(i)
    run_isolated(stdout, entries; time_limit = tl, coeffs = coeffs)
  else
    open(args[i + 1], "w") do io
      run_isolated(io, entries; time_limit = tl, coeffs = coeffs)
    end
    @info "wrote $(args[i + 1])"
  end
end

main(ARGS)
