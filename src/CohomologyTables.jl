#=
Cup product multiplication tables for spaces with a simplicial model, built on
OSCAR's experimental hypercomplex machinery (`DoubleAndHyperComplexes`).

    julia> include("src/CohomologyTables.jl"); using .CohomologyTables
    julia> using Oscar

    julia> X = simplicial_cohomology_ring("K3", space("K3"), GF(7))
    julia> print_report(X)      # the 22 x 22 intersection form of K3 mod 7

    julia> K = simplicial_complex([[1,2,3], [2,3,4], [1,3,4], [1,2,4]])   # your own
    julia> print_report(simplicial_cohomology_ring("my sphere", K, QQ))

The four pieces:

  sources.jl      where complexes come from: a local SageMath, Frank Lutz's
                  Manifold Page (downloaded and cached outside the repo), arXiv
  spaces.jl       the catalogue -- one model per topological type
  ring_tables.jl  the computation: canonical bases per degree, cup products as
                  coordinates, and the wall-clock budget
  report.jl       printing, ring identification, consistency checks

`scripts/run_tables.jl` is the batch driver. See README.md for the full space
list and the limitations worth knowing about.
=#
module CohomologyTables

using Oscar
using Downloads
using Dates
using SHA

include("sources.jl")
include("spaces.jl")
include("ring_tables.jl")
include("report.jl")
include("records.jl")

export
  # sources
  sage, sage_available, sage_version, lutz, lutz_labels, lutz_header, lutz_url,
  lutz_cached_at, arxiv_cached_at, bagchi_datta_cp3,
  parse_gap_blocks, relabel_facets, fetch_text, cache_dir, clear_cache,
  # the catalogue
  SpaceEntry, Recipe, catalogue, catalogue_names, catalogue_table, find_space,
  catalogue_space, provenances, build, is_reproducible, builtin_recipe, sage_recipe,
  lutz_recipe,
  arxiv_recipe, point_space, sphere, kuehnel_real_projective_space,
  # the ring
  CohomologyRing, GradedBasis, simplicial_cohomology_ring, normalize_to_powers,
  graded_basis, basis_indices, top_degree,
  class_coords, cup, cup_power,
  Deadline, TimeLimitExceeded, DEFAULT_TIME_LIMIT, check_deadline, no_deadline,
  # reporting
  cup_product_table, identify_ring, check_ring, print_report,
  format_combination, group_symbol, ring_symbol, class_label,
  # generated records
  canonical_facets, canonical_form, facets_sha256, ring_record, json_string,
  record_filename, record_prefix, write_record, CANONICAL_FORM_ID, RECORD_SCHEMA_ID

end # module
