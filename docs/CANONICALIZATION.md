# Canonical form of a simplicial complex

Records identify their input complex by a SHA-256 over a canonical rendering of
its facets, never by storing the facets. This document fixes that rendering.

The rule is versioned by `CANONICAL_FORM_ID`, currently
`cohomology-tables/canonical-facets/1`. If the rule ever changes, that string
changes with it and old records stay interpretable.

## The rule

Given a simplicial complex with vertex labels `L` and facets given as subsets of
`L`:

1. **Order the labels.** Sort `L` under the total order for its source (below)
   and assign `0, 1, 2, …` in that order. Every label is replaced by its index,
   so the canonical form is always 0-based.
2. **Order within each facet.** Sort each facet's indices ascending.
3. **Order the facets.** Sort the facets lexicographically as sequences of
   integers: compare element by element; if one is a prefix of the other, the
   shorter comes first.
4. **Render.** Each facet becomes its indices in decimal, comma-separated, with
   no spaces. Facets are joined with `\n` (U+000A). There is **no trailing
   newline**.
5. **Hash.** SHA-256 of that string encoded as UTF-8, lowercase hex.

## The label order, per source

Step 1 needs a total order on labels that may not be integers, so it is fixed
per source and each is deterministic:

| source | labels | order |
|---|---|---|
| SageMath | arbitrary hashable Python objects — integers, tuples such as `(2, 4)`, elements of `Z/n` | ascending by Python `str(label)`, compared as Python strings (code point order) |
| Lutz's Manifold Page | integers, sometimes 0-based, sometimes with gaps | ascending as integers |
| built here | integers `1..n` | ascending as integers |

Sage's labels are the reason this has to be stated at all: `MatchingComplex(n)`
has tuples for vertices and `SumComplex(n, A)` has residues mod `n`, neither of
which has a canonical integer order of its own. Sorting by `str` is applied in
the Sage process itself, before anything crosses into Julia.

## Why the order matters beyond the hash

The Alexander–Whitney cup product is defined relative to a total order on the
vertices: on dual basis simplices `sigma` and `tau` it is non-zero only when
`max(sigma) == min(tau)`. Different orders give isomorphic rings but different
cochain-level representatives, so a record's cocycle representatives are only
meaningful together with the order that produced them. Every record therefore
carries a `vertex_order_convention` field naming this document, and the
representatives it lists are expressed in canonical 0-based labels.

## When the hash is not reproducible

The scheme assumes replaying a recipe rebuilds the *same labelled* complex. That
holds for every source here except Sage's `SurfaceOfGenus`, which returns a
differently numbered triangulation on each process: same f-vector, same
homology, different vertex labels. Measured over repeated runs,
`SurfaceOfGenus(3)` alternates between at least two labellings.
`SurfaceOfGenus(2)` was stable over three runs, but it is the same code path, so
the whole family is treated as non-reproducible rather than trusting weak
evidence.

`is_reproducible(recipe)` reports this, `records/MANIFEST.tsv` carries it in the
`reproducible` column, and affected records set `reproducible_labelling: false`
and explain it in `source.labelling_caveat`. For those entries
`scripts/regenerate_sage.jl` compares f-vector and homology, which are invariant
under relabelling, instead of the hash — and still fails if *those* differ,
since that would mean Sage is building a different space rather than a different
numbering.

The cohomology ring is unaffected: it is a topological invariant. What is tied
to the labelling is the cocycle representatives, which belong to the particular
numbering the record was computed from.

## Worked example

The boundary of the tetrahedron, `S^2`, as built by `sphere(2)`:

```
0,1,2
0,1,3
0,2,3
1,2,3
```

with SHA-256

```
41d49c975f160871b5e3210fa4e0d84f9b4abfb79c57c63a692ab77ead820200
```

Reproduce it with:

```julia
julia> K = catalogue_space("S^2");
julia> canonical_form(K)     # the string above
julia> facets_sha256(K)      # the digest above
```
