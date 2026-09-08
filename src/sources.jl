#=
Where the simplicial complexes come from.

Nothing is stored in this repository. Every model is either constructed on the
spot, asked of a local SageMath, or downloaded from Frank Lutz's Manifold Page
and cached outside the repo (see `cache_dir`).

Three source functions:

  sage("Torus()")            -> SimplicialComplex   via a local `sage` process
  lutz("K3_16")              -> SimplicialComplex   one file = one complex
  lutz("flat_3manifolds", "T^3")
                             -> SimplicialComplex   one file = many complexes
=#

### -------------------------------------------------------------------------
### Cache
### -------------------------------------------------------------------------

"""
    cache_dir() -> String

Where downloaded triangulations are kept. Defaults to a scratchspace in the
Julia depot, never inside this repository; override with the
`COHOMOLOGY_TABLES_CACHE` environment variable. Safe to delete at any time.
"""
function cache_dir()
  dir = get(ENV, "COHOMOLOGY_TABLES_CACHE",
            joinpath(first(DEPOT_PATH), "scratchspaces", "cohomology-tables"))
  isdir(dir) || mkpath(dir)
  return dir
end

"""
    clear_cache()

Delete everything `cache_dir` holds. The next call re-downloads.
"""
function clear_cache()
  dir = cache_dir()
  isdir(dir) && rm(dir; recursive = true)
  return nothing
end

const LUTZ_STELLAR = "https://www3.math.tu-berlin.de/IfM/Nachrufe/Frank_Lutz/stellar"
const LUTZ_LIBRARY = "$LUTZ_STELLAR/library_of_triangulations"

"""
    fetch_text(url; name) -> String

Download `url` once, cache it under [`cache_dir`](@ref), and return the body
with non-ASCII bytes dropped (several of Lutz's files carry Latin-1 bytes in
author names, and every field we parse is ASCII).
"""
function fetch_text(url::AbstractString; name::AbstractString = basename(url))
  path = joinpath(cache_dir(), name)
  if !isfile(path)
    try
      Downloads.download(url, path)
    catch e
      rm(path; force = true)
      throw(ErrorException("could not download $url -- $(sprint(showerror, e))"))
    end
  end
  return String(filter(<=(0x7f), read(path)))
end

### -------------------------------------------------------------------------
### Parsing the GAP-style facet lists
### -------------------------------------------------------------------------

"""
    relabel_facets(facets) -> Vector{Vector{Int}}

Renumber vertices to `1:n`, preserving their order. Necessary because the
sources disagree: some files are 0-based, `contractible_non_5_ball` numbers 5013
vertices up to 7076, and Sage hands out tuples.

The relabelling is order-preserving, which matters: the Alexander-Whitney cup
product is defined relative to a total order on the vertices. Any order gives an
isomorphic ring, but a deterministic one makes results reproducible.
"""
function relabel_facets(facets::Vector{Vector{Int}})
  isempty(facets) && return facets
  labels = sort!(unique!(reduce(vcat, deepcopy(facets))))
  index = Dict(v => i for (i, v) in enumerate(labels))
  return [sort!([index[v] for v in f]) for f in facets]
end

"""
    parse_gap_blocks(text) -> Vector{Pair{String,Vector{Vector{Int}}}}

Pull every `label => facets` pair out of a Manifold Page file.

The pages use three layouts and this handles all of them:

  * `facets:=[ [1,2,3], ... ];`       -- one complex per file (the library)
  * `##  T^3   f = (15,105,180,90)`   -- many per file (geometric 3-manifolds);
    `facets:=[...]` or `facets=[...]`    the label comes from the `##` line above
  * `CP^2=[[1,2,3,4,5], ... ]`        -- many per file (small_4_manifolds)

Bracket nesting is tracked by hand rather than by regex, since the lists run to
thousands of entries across many lines.
"""
function parse_gap_blocks(text::AbstractString)
  out = Pair{String,Vector{Vector{Int}}}[]
  s = String(text)
  # An assignment `<label> :=` or `<label> =` introducing a bracketed list.
  # Matched rather than searching for "[[" because the library files space
  # their lists out as `facets:=[ [ 1, 2, 4 ], ... ];`.
  for m in eachmatch(r"([A-Za-z_(][A-Za-z0-9_^#()\.\-]*)\s*:?=\s*\[", s)
    open_at = m.offset + length(m.match) - 1        # index of the opening '['
    label = m.captures[1]
    if label == "facets"
      # The label lives on the nearest preceding `##  <label>   f = (...)` line.
      # Requiring the `f =` keeps the files' `###...###` banners out of it.
      before = s[firstindex(s):(open_at - 1)]
      hits = collect(eachmatch(r"(?m)^##\s+(\S+)\s+f\s*=", before))
      label = isempty(hits) ? "facets" : last(hits).captures[1]
    end
    # Scan to the bracket that closes the outer list.
    depth, k = 0, open_at
    while k <= lastindex(s)
      c = s[k]
      c == '[' && (depth += 1)
      c == ']' && (depth -= 1; depth == 0 && break)
      k = nextind(s, k)
    end
    depth == 0 || throw(ArgumentError("unbalanced brackets after `$label`"))
    body = s[open_at:k]
    facets = [parse.(Int, [x.match for x in eachmatch(r"-?\d+", inner.captures[1])])
              for inner in eachmatch(r"\[([^\[\]]*)\]", body)]
    filter!(!isempty, facets)
    isempty(facets) || push!(out, label => relabel_facets(facets))
  end
  return out
end

### -------------------------------------------------------------------------
### Lutz's Manifold Page
### -------------------------------------------------------------------------

"""
    lutz(name) -> SimplicialComplex
    lutz(file, label) -> SimplicialComplex

A triangulation from Frank Lutz's Manifold Page,
<https://www3.math.tu-berlin.de/IfM/Nachrufe/Frank_Lutz/stellar/>, downloaded on
first use and cached outside the repo.

The one-argument form takes a file holding a single complex, looking first in
the Benedetti-Lutz Library of Triangulations and then in the Further Examples at
the top level -- `lutz("K3_16")`, `lutz("RP3")`.

The two-argument form picks one complex out of a file holding many, by the label
in its `##` header -- `lutz("flat_3manifolds", "T^3")`,
`lutz("spherical_3manifolds", "L_3_1")`. [`lutz_labels`](@ref) lists what a
given file contains.
"""
function lutz(name::AbstractString)
  blocks = _lutz_blocks(name)
  isempty(blocks) && throw(ArgumentError("no facet list found in $name.txt"))
  length(blocks) == 1 || throw(ArgumentError(
    "$name.txt holds $(length(blocks)) complexes; use lutz(\"$name\", label) " *
    "with one of: $(join(first.(blocks), ", "))"))
  return simplicial_complex(last(only(blocks)))
end

function lutz(file::AbstractString, label::AbstractString)
  blocks = _lutz_blocks(file)
  i = findfirst(p -> first(p) == label, blocks)
  isnothing(i) && throw(ArgumentError(
    "$file.txt has no complex labelled \"$label\"; available: " *
    join(first.(blocks), ", ")))
  return simplicial_complex(last(blocks[i]))
end

"""
    lutz_labels(file) -> Vector{String}

The labels of the complexes in one Manifold Page file.
"""
lutz_labels(file::AbstractString) = first.(_lutz_blocks(file))

function _lutz_blocks(name::AbstractString)
  for base in (LUTZ_LIBRARY, LUTZ_STELLAR)
    try
      return parse_gap_blocks(fetch_text("$base/$name.txt"; name = "lutz_$name.txt"))
    catch
      continue
    end
  end
  throw(ArgumentError("could not fetch \"$name.txt\" from the Manifold Page"))
end

"""
    lutz_header(name) -> String

The `###`/`##` comment header of a Manifold Page file: description, f-vector,
integral homology and references, as its authors wrote them. Useful as an
independent check on what we compute.
"""
function lutz_header(name::AbstractString)
  for base in (LUTZ_LIBRARY, LUTZ_STELLAR)
    try
      text = fetch_text("$base/$name.txt"; name = "lutz_$name.txt")
      lines = split(text, '\n')
      stop = something(findfirst(l -> occursin(r"\[\[", l), lines), length(lines) + 1)
      return join(lines[1:(stop - 1)], '\n')
    catch
      continue
    end
  end
  throw(ArgumentError("could not fetch \"$name.txt\" from the Manifold Page"))
end

### -------------------------------------------------------------------------
### SageMath
### -------------------------------------------------------------------------

"""
    sage(expr; command = "sage") -> SimplicialComplex

Build one of SageMath's example complexes, e.g. `sage("Torus()")`,
`sage("RealProjectiveSpace(4)")`, `sage("MooreSpace(5)")`. The argument is
whatever follows `simplicial_complexes.` in Sage's
[catalogue](https://doc.sagemath.org/html/en/reference/topology/sage/topology/simplicial_complex_examples.html).

Requires `sage` on the `PATH`; set the `SAGE` environment variable to point
elsewhere. Sage's vertices can be tuples or integers mod n, so the relabelling
to `1:n` happens on the Sage side, ordered by `str` for reproducibility.
"""
function sage(expr::AbstractString; command::AbstractString = get(ENV, "SAGE", "sage"))
  script = """
  K = simplicial_complexes.$expr
  V = sorted(K.vertices(), key=str)
  idx = {v: i + 1 for i, v in enumerate(V)}
  print([sorted(idx[v] for v in f) for f in K.facets()])
  """
  out = try
    read(`$command -c $script`, String)
  catch e
    throw(ErrorException(
      "could not run `$command` for simplicial_complexes.$expr -- " *
      "is SageMath installed and on the PATH? ($(sprint(showerror, e)))"))
  end
  m = match(r"\[\[.*\]\]", out)
  isnothing(m) && throw(ErrorException(
    "unexpected output from Sage for simplicial_complexes.$expr:\n$out"))
  facets = [parse.(Int, [x.match for x in eachmatch(r"\d+", inner.captures[1])])
            for inner in eachmatch(r"\[([^\[\]]*)\]", m.match)]
  return simplicial_complex(relabel_facets(filter!(!isempty, facets)))
end

"Whether a working `sage` is on the PATH."
function sage_available(; command::AbstractString = get(ENV, "SAGE", "sage"))
  return try
    success(`$command --version`)
  catch
    false
  end
end

### -------------------------------------------------------------------------
### arXiv
### -------------------------------------------------------------------------

"""
    bagchi_datta_cp3() -> SimplicialComplex

The 18-vertex, 622-facet triangulation of `CP^3` of Bagchi and Datta,
[arXiv:1012.3235](https://arxiv.org/abs/1012.3235), read out of the paper's own
LaTeX source, where the facets are listed as words in `a_1..a_9`, `b_1..b_9`.
Any triangulation of `CP^3` needs at least 17 vertices (Arnoux-Marin), so this
is close to vertex-minimal.

Neither SageMath nor the Manifold Page carries `CP^3`, and this repository
stores no facet lists, so the paper is the source. The e-print is gzipped and
is decompressed with `gzip -dc`; on a system without `gzip` this is the one
model that will not build.
"""
function bagchi_datta_cp3()
  path = joinpath(cache_dir(), "arxiv_1012.3235.tex")
  if !isfile(path)
    gz = joinpath(cache_dir(), "arxiv_1012.3235.gz")
    Downloads.download("https://arxiv.org/e-print/1012.3235", gz)
    try
      open(path, "w") do io
        run(pipeline(`gzip -dc $gz`; stdout = io))
      end
    catch e
      rm(path; force = true)
      throw(ErrorException(
        "could not decompress the arXiv source for CP^3 -- is `gzip` available? " *
        "($(sprint(showerror, e)))"))
    finally
      rm(gz; force = true)
    end
  end
  text = String(filter(<=(0x7f), read(path)))
  start = findfirst("complete list of the 622 facets", text)
  isnothing(start) && throw(ErrorException(
    "the facet list is not where expected in arXiv:1012.3235; the source may have changed"))
  words = [m.captures[1] for m in eachmatch(r"((?:[ab]_\d+){7})\s*\$", text[first(start):end])]
  facets = map(words) do w
    sort!([(c == "a" ? 0 : 9) + parse(Int, k)
           for (c, k) in (m.captures for m in eachmatch(r"([ab])_(\d+)", w))])
  end
  unique!(facets)
  length(facets) == 622 || throw(ErrorException(
    "expected 622 facets for CP^3, parsed $(length(facets))"))
  return simplicial_complex(facets)
end

### -------------------------------------------------------------------------
### Provenance helpers used by the generated records
### -------------------------------------------------------------------------

"""
    lutz_url(file) -> String

The URL a Manifold Page file is fetched from. Files live either in the Library
of Triangulations or at the top level of the stellar directory; this reports
whichever one actually resolved, having fetched it if necessary.
"""
function lutz_url(file::AbstractString)
  for base in (LUTZ_LIBRARY, LUTZ_STELLAR)
    url = "$base/$file.txt"
    try
      fetch_text(url; name = "lutz_$file.txt")
      return url
    catch
      continue
    end
  end
  throw(ArgumentError("could not locate \"$file.txt\" on the Manifold Page"))
end

"""
    sage_version(; command) -> String

The banner string of the local SageMath, e.g. `"SageMath version 10.3, Release
Date: 2024-03-19"`.

Recorded with every Sage-sourced model: Sage's example definitions and vertex
labellings can change between releases, so a record is only reproducible against
a stated version. Returns `"unavailable"` if Sage cannot be run, so that record
generation on a machine without Sage still says so honestly rather than failing.
"""
function sage_version(; command::AbstractString = get(ENV, "SAGE", "sage"))
  return try
    strip(first(split(read(`$command --version`, String), '\n')))
  catch
    "unavailable"
  end
end

"""
    lutz_cached_at(file) -> Union{String,Nothing}

The date the cached copy of a Manifold Page file was written, as `YYYY-MM-DD`,
or `nothing` if it is not cached. This is an observation of when the download
happened, not a guess: it reads the cache file's modification time.
"""
function lutz_cached_at(file::AbstractString)
  path = joinpath(cache_dir(), "lutz_$file.txt")
  isfile(path) || return nothing
  return Dates.format(Dates.unix2datetime(mtime(path)), "yyyy-mm-dd")
end

"""
    arxiv_cached_at(id) -> Union{String,Nothing}

As [`lutz_cached_at`](@ref), for a cached arXiv source.
"""
function arxiv_cached_at(id::AbstractString)
  path = joinpath(cache_dir(), "arxiv_$id.tex")
  isfile(path) || return nothing
  return Dates.format(Dates.unix2datetime(mtime(path)), "yyyy-mm-dd")
end
