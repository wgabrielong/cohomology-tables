# Provenance of the input complexes

This repository **redistributes no simplicial complexes**. Every model is either
constructed by code here, obtained by invoking a SageMath constructor on your
machine, or downloaded at run time into a gitignored cache. What is stored here
is the *recipe* to rebuild each complex and a SHA-256 of its canonical form (see
[CANONICALIZATION.md](CANONICALIZATION.md)), so a rebuild that does not
reproduce the input fails loudly.

To be precise about what "no complexes" means: no record contains a facet list.
Records do name the individual faces in the support of each cocycle
representative — the faces a chosen cochain is non-zero on. That is a small
subset of the faces, only in dimensions carrying cohomology, and it does not
determine the complex. If you consider even that too close to the upstream data
for your purposes, drop the `cocycle_representatives` field when regenerating.

The generated records under `records/` are the author's own work and are
released under CC0-1.0. The package source is GPL-3.0-or-later. Neither licence
is asserted over any upstream complex.

---

## Source (a) — Frank H. Lutz's Manifold Page

<https://www3.math.tu-berlin.de/IfM/Nachrufe/Frank_Lutz/stellar/>

Used for: the Library of Triangulations, the Further Examples page, and the
geometric 3-manifold catalogues (spherical, flat, Nil, `S^2 x R`, `H^2 x R`,
hyperbolic, Brieskorn homology spheres, connected sums), plus `small_4_manifolds`
and `S3xS2`.

The page carries material by **Thom Sulanke**, and the Library of Triangulations
is joint work with **Bruno Benedetti**.

**Licence status: no terms stated.** The page states no licence for the
triangulation data, and Frank Lutz has died, so no licence can now be sought
from him. The author of this repository therefore holds no licence to grant over
this material and **applies no SPDX identifier to it**. Nothing from this source
is committed here; `scripts/fetch_sources.jl` downloads it on demand into
`cache_dir()`, which is outside the repository and gitignored. If you
redistribute a cache you have populated, that is your decision to justify, not
one this repository makes for you.

### What to cite

| catalog | cite |
|---|---|
| Library of Triangulations (`lutz("K3_16")` and friends) | B. Benedetti and F. H. Lutz, *Random discrete Morse theory and a new library of triangulations*, Exp. Math. **23**, 66–94 (2014) |
| `sphere_bundles`, `connected_sums_3d`, and the 3-manifold f-vector data | F. H. Lutz, T. Sulanke and E. Swartz, *f-vectors of 3-manifolds*, Electron. J. Combin. **16**(2), #R13 (2009) |
| `hyperbolic_dodecahedral_space` | F. H. Lutz, T. Sulanke and E. Swartz, *f-vectors of 3-manifolds*, Electron. J. Combin. **16**(2), #R13 (2009); and F. H. Lutz, *Triangulated manifolds with few vertices: geometric 3-manifolds*, preprint, 48 pp., 2003, [arXiv:math/0311116](https://arxiv.org/abs/math/0311116) |
| `non_PL_5_sphere`, `Poincare_3_sphere` | A. Björner and F. H. Lutz, *Simplicial manifolds, bistellar flips and a 16-vertex triangulation of the Poincaré homology 3-sphere*, Exp. Math. **9**, 275–289 (2000) |
| `HMT_4`, `HMT_8`, `HMT_16`, `HMT_32` | D. Lofano and F. H. Lutz, *Hadamard matrix torsion*, arXiv:2109.13052 (2021) |
| `Abalone`, `BH_*`, `d2_n8_*torsion`, `d2_n9_5torsion` | B. Benedetti, C. Lai, D. Lofano and F. H. Lutz, *Random simple-homotopy theory*, arXiv:2107.09862 (2021) |
| `PG64_*`, `PG128_*`, `AG_5_3_*`, `non_4_2_colorable` | F. H. Lutz and J. M. Møller, *Chromatic numbers of simplicial manifolds*, Beitr. Algebra Geom. **61**, 419–453 (2020) |
| geometric 3-manifold catalogues (`flat_3manifolds`, `spherical_3manifolds`, `nil_3manifolds`, `S2xR_spaces`, `H2xR_spaces`, `hyperbolic_3manifolds`, `homology_3spheres`) | F. H. Lutz, T. Sulanke and E. Swartz, *f-vectors of 3-manifolds*, Electron. J. Combin. **16**(2), #R13 (2009); and F. H. Lutz, *Triangulated manifolds with few vertices: geometric 3-manifolds*, preprint, 48 pp., 2003, [arXiv:math/0311116](https://arxiv.org/abs/math/0311116) |
| `small_4_manifolds` | no publication is given for these; cite the page itself, <https://www3.math.tu-berlin.de/IfM/Nachrufe/Frank_Lutz/stellar/further-examples.html> |
| `S3xS2` | F. H. Lutz, T. Sulanke and E. Swartz, *f-vectors of 3-manifolds*, Electron. J. Combin. **16**(2), #R13 (2009) |

`lutz_header(name)` prints any file's own header, which is where its authors put
the description, f-vector, integral homology and references. That header is the
authoritative attribution for a given file; the table above is a convenience.

---

## Source (b) — SageMath

`sage.topology.simplicial_complex_examples`, via
<https://doc.sagemath.org/html/en/reference/topology/sage/topology/simplicial_complex_examples.html>

Used for: spheres and projective spaces, surfaces, Moore spaces, `K3Surface`,
`QuaternionicProjectivePlane`, `PoincareHomologyThreeSphere`, the chessboard,
matching and sum complexes, and `FareyMap`.

**Licence status: GPL-2.0-or-later** (SageMath library code). That is
upgrade-compatible with GPL-3.0-or-later, so there is no conflict with this
package's Oscar.jl dependency.

**Facet lists hardcoded inside that Sage module are GPL'd source text and are
not copied into this repository.** Complexes are obtained only by invoking the
constructors in a local Sage process (`src/sources.jl`, function `sage`), and
what is stored here is the invocation string, the Sage version it was recorded
under, and a hash. `scripts/regenerate_sage.jl` replays the invocations and
verifies the hashes.

Sage's example definitions and vertex labellings can change between releases, so
the recorded version string is load-bearing and a hash mismatch is treated as an
error rather than a warning.

One family is non-reproducible even within a single Sage version:
`SurfaceOfGenus` returns a differently labelled complex on each process. Those
records say so (`reproducible_labelling: false`) and are verified against
isomorphism invariants instead. See
[CANONICALIZATION.md](CANONICALIZATION.md#when-the-hash-is-not-reproducible).

### Upstream attribution inside Sage

Several Sage examples are themselves due to the authors above, or to Kühnel, and
should carry that attribution as well as Sage's:

| Sage constructor | originally due to |
|---|---|
| `ComplexProjectivePlane()` | Kühnel and Banchoff |
| `K3Surface()` | Casella and Kühnel; labelling from Spreer and Kühnel |
| `RealProjectiveSpace(3)` | Walkup; explicit facets given by Lutz |
| `RealProjectiveSpace(4)` | Datta |
| `RealProjectiveSpace(n)` for `n >= 5` | Kühnel |
| `QuaternionicProjectivePlane()` | Brehm and Kühnel, whose description Sage implements; Gorodkov proved it triangulates `HP^2` |
| `PoincareHomologyThreeSphere()` | Björner and Lutz |
| `BarnetteSphere()` | Barnette |
| `BrucknerGrunbaumSphere()` | Brückner and Grünbaum |
| `RudinBall()` | Rudin |
| `ZieglerBall()` | Ziegler |
| `DunceHat()` | Zeeman; the triangulation is Hachimori's |

For `QuaternionicProjectivePlane()` specifically, Sage builds the one of the
three triangulations in Brehm–Kühnel that has a transitive automorphism group
(isomorphic to `A_5`), and it was only later proved to be `HP^2`:

* U. Brehm and W. Kühnel, *15-vertex triangulations of an 8-manifold*,
  Math. Ann. **294** (1992), no. 1, 167–193.
* D. Gorodkov, *A 15-vertex triangulation of the quaternionic projective plane*,
  Discrete Comput. Geom. **62** (2019), 348–373,
  [doi:10.1007/s00454-018-00055-w](https://doi.org/10.1007/s00454-018-00055-w).

---

## Source (c) — arXiv

`CP^3` is read out of the LaTeX source of

* B. Bagchi and B. Datta, *A triangulation of CP^3 as symmetric cube of S^2*,
  [arXiv:1012.3235](https://arxiv.org/abs/1012.3235),

which lists its 622 facets in the paper text. The source is downloaded on demand
into the same gitignored cache and is not committed here. Cite the paper.

**Licence status: arXiv non-exclusive distribution licence**
(<http://arxiv.org/licenses/nonexclusive-distrib/1.0/>). Under it the submitter
grants arXiv "a perpetual, non-exclusive license to distribute this article".
That is a licence to *arXiv*, not to the public: it confers no onward
redistribution right on third parties. This repository therefore does not
redistribute the source or the facet list, and downloads the e-print on demand
into the gitignored cache like everything else.

---

## Fields still to be established

None. Every citation above is confirmed, `REUSE.toml` carries the copyright line
`2026 Wern Juin Gabriel Ong`, and `CITATION.cff` describes this repository
itself.

No dates, names or email addresses have been invented anywhere in this
repository. Where a date appears in a record's `date_accessed` it is read from
the modification time of the cached download, i.e. it is an observation of when
that fetch happened on the machine that generated the record.
