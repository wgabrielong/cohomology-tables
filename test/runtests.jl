#=
Checks against the textbook answers.

    julia test/runtests.jl              # offline core only
    julia test/runtests.jl --online     # + SageMath and the Manifold Page

Everything here is a theorem, not a regression baseline: if one of these fails,
either a model is wrong or the cup product machinery is.

The offline set uses only models built in this repo, so it runs without Sage and
without a network. `--online` adds the catalogue entries that need one or both;
those are also the ones that check our results against the homology their
sources publish.
=#

using Test
using Oscar
include(joinpath(@__DIR__, "..", "src", "CohomologyTables.jl"))
using .CohomologyTables

const ONLINE = ("--online" in ARGS) || get(ENV, "COHOMOLOGY_TABLES_ONLINE", "") == "1"

"The list of cohomology groups, as strings: `[\"Z\", \"0\", \"Z/2\"]`."
groups(X) = [group_symbol(X.coeff_ring, graded_basis(X, d).orders) for d in 0:top_degree(X)]

"Dimensions of `H^*` over a field."
dims(X) = [length(graded_basis(X, d)) for d in 0:top_degree(X)]

@testset "CohomologyTables" begin

  @testset "the point" begin
    X = simplicial_cohomology_ring("point", point_space(), ZZ)
    @test groups(X) == ["Z"]
    @test identify_ring(X) == "Z"
    @test isempty(check_ring(X))
  end

  @testset "S^$n" for n in 1:4
    X = simplicial_cohomology_ring("S^$n", sphere(n), ZZ)
    @test groups(X) == ["Z"; fill("0", n - 1); "Z"]
    @test identify_ring(X) == "Z[x$n]/(x$n^2),  |x$n| = $n"
    @test isempty(check_ring(X))
  end

  @testset "RP^n from the built-in Kuehnel construction" for n in 2:3
    K = kuehnel_real_projective_space(n)
    @test n_vertices(K) == 2^(n + 1) - 1
    X = simplicial_cohomology_ring("RP^$n", K, GF(2))
    # H^*(RP^n; F_2) = F_2[x]/(x^(n+1)), |x| = 1
    @test dims(X) == fill(1, n + 1)
    @test identify_ring(X) == "F_2[x1]/(x1^$(n + 1)),  |x1| = 1"
    @test isempty(check_ring(X))
  end

  @testset "facet parsing" begin
    blocks = parse_gap_blocks("""
      ##  T^3   f = (4,6,4), g_2 = 0.
      facets:=[ [ 1, 2, 3 ], [ 1, 2, 4 ] ];
      CP^2=[[0,1,2],[0,1,3]]
      """)
    @test first.(blocks) == ["T^3", "CP^2"]
    @test last(blocks[1]) == [[1, 2, 3], [1, 2, 4]]
    @test last(blocks[2]) == [[1, 2, 3], [1, 2, 4]]      # 0-based input relabelled
    @test relabel_facets([[5, 9], [9, 20]]) == [[1, 2], [2, 3]]
  end

  @testset "canonical form and hashing" begin
    K = sphere(2)
    # 0-based, vertices ascending inside facets, facets lexicographic, no
    # trailing newline. See docs/CANONICALIZATION.md.
    @test canonical_facets(K) == [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    @test canonical_form(K) == "0,1,2\n0,1,3\n0,2,3\n1,2,3"
    @test !endswith(canonical_form(K), "\n")
    @test facets_sha256(K) ==
          "41d49c975f160871b5e3210fa4e0d84f9b4abfb79c57c63a692ab77ead820200"
    # the hash is over the complex, so relabelling the input does not change it
    shifted = simplicial_complex(relabel_facets(
      [[v + 100 for v in f] for f in [[1, 2, 3], [1, 2, 4], [1, 3, 4], [2, 3, 4]]]))
    @test facets_sha256(shifted) == facets_sha256(K)
  end

  @testset "records" begin
    e = find_space("S^2")
    X = simplicial_cohomology_ring(e.name, build(e), ZZ)
    rec = ring_record(e, X)
    @test rec["schema"] == RECORD_SCHEMA_ID
    @test rec["complex"]["canonicalization"] == CANONICAL_FORM_ID
    @test rec["complex"]["facets_sha256"] == facets_sha256(build(e))
    @test rec["source"]["kind"] == "builtin"
    @test rec["ring_presentation"] == "Z[x2]/(x2^2),  |x2| = 2"
    @test isempty(rec["consistency_checks"])
    @test rec["structure_constants_complete"]
    # a record carries no facets, only the hash
    text = json_string(rec)
    @test !occursin("\"facets\"", text)
    @test occursin(facets_sha256(build(e)), text)
    # filenames stay distinct after the slug flattening
    names = catalogue_names()
    @test allunique([record_filename(n, ZZ) for n in names])
    @test !occursin(r"[^A-Za-z0-9._]", record_filename("(S^2xS^1)#2#RP^3", ZZ))
  end

  @testset "recipes rebuild what the catalogue says" begin
    for n in ["S^3", "point", "S^0"]
      e = find_space(n)
      @test e.recipe.kind == "builtin"
      @test build(e.recipe) isa SimplicialComplex
    end
    @test find_space("K3").recipe == sage_recipe("K3Surface()")
    @test find_space("T^3").recipe == lutz_recipe("flat_3manifolds", "T^3")
    @test find_space("CP^3").recipe == arxiv_recipe("1012.3235")
    @test_throws ArgumentError build(Recipe("nonsense", "", "", ""))
  end

  @testset "reproducibility of recipes" begin
    # Sage's SurfaceOfGenus relabels on every process, so its hash is not a
    # usable fingerprint; everything else here is reproducible.
    @test !is_reproducible(sage_recipe("SurfaceOfGenus(3)"))
    @test !is_reproducible(sage_recipe("SurfaceOfGenus(4, orientable=False)"))
    @test is_reproducible(sage_recipe("K3Surface()"))
    @test is_reproducible(lutz_recipe("K3_16"))
    @test is_reproducible(builtin_recipe("sphere(3)"))
    @test !is_reproducible(find_space("Sigma_3"))
    @test is_reproducible(find_space("K3"))
    # record_prefix must not confuse a name with one it is a prefix of
    @test record_prefix("Sigma_4") != record_prefix("Sigma_4 x S^1")
    @test !startswith(record_prefix("Sigma_4 x S^1"), record_prefix("Sigma_4"))
  end

  @testset "catalogue integrity" begin
    C = catalogue()
    @test length(C) > 200
    @test allunique([e.name for e in C])
    @test all(e -> e.provenance in provenances(), C)
    @test all(e -> !isempty(e.type), C)
    @test all(e -> !isempty(e.coeffs), C)
    @test find_space("K3").name == "K3"
    @test_throws ArgumentError find_space("no such space")
    # one model per topological type: the well-known duplicates are gone
    names = [e.name for e in C]
    @test !("K3_16" in names) && !("K3_17" in names)     # both are "K3"
    @test !("BarnetteSphere" in names)                    # that is "S^3"
    @test !("dunce_hat" in names)                         # that is "point"
    @test !("Hom_C5_K4" in names)                         # that is "RP^3"
    @test !("d2n12g6" in names)                           # that is "Sigma_6"
    @test !any(startswith(n, "FareyMap") for n in names)  # genus 0/3/26 instead
    @test !("M(Z/2,1)" in names)                          # that is "RP^2"
    # the Moore family, one entry per (q, n)
    @test [n for n in names if startswith(n, "M(Z/")] ==
          ["M(Z/3,1)", "M(Z/4,1)", "M(Z/5,2)", "M(Z/7,3)", "M(Z/8,4)", "M(Z/9,2)"]
  end

  @testset "time budget" begin
    @test simplicial_cohomology_ring("S^4", sphere(4), ZZ; time_limit = 600) isa CohomologyRing
    @test_throws TimeLimitExceeded simplicial_cohomology_ring(
      "S^4", sphere(4), ZZ; deadline = Deadline(time() - 1000, 1.0))
  end

  if !ONLINE
    @info "offline run: pass --online to also check Sage and Manifold Page models"
  else
    @testset "Sage models" begin
      @test sage_available()
      # RP^n over F_2 is the truncated polynomial ring F_2[x]/(x^(n+1))
      for n in 2:4
        X = simplicial_cohomology_ring("RP^$n", catalogue_space("RP^$n"), GF(2))
        @test dims(X) == fill(1, n + 1)
        @test identify_ring(X) == "F_2[x1]/(x1^$(n + 1)),  |x1| = 1"
      end
      # H^*(RP^n; Z)
      @test groups(simplicial_cohomology_ring("RP^3", catalogue_space("RP^3"), ZZ)) ==
            ["Z", "0", "Z/2", "Z"]
      @test groups(simplicial_cohomology_ring("RP^4", catalogue_space("RP^4"), ZZ)) ==
            ["Z", "0", "Z/2", "0", "Z/2"]
      # CP^2 = Z[y]/(y^3), |y| = 2
      XC = simplicial_cohomology_ring("CP^2", catalogue_space("CP^2"), ZZ)
      @test identify_ring(XC) == "Z[x2]/(x2^3),  |x2| = 2"
      # torus: exterior algebra on two degree-1 classes
      XT = simplicial_cohomology_ring("T^2", catalogue_space("T^2"), ZZ)
      @test groups(XT) == ["Z", "Z^2", "Z"]
      @test cup(XT, (1, 1), (1, 1)) == [zero(ZZ)]
      @test cup(XT, (1, 1), (1, 2)) == [-c for c in cup(XT, (1, 2), (1, 1))]
      @test isnothing(identify_ring(XT))
      # K3: H^2 has rank 22
      @test dims(simplicial_cohomology_ring("K3", catalogue_space("K3"), GF(7))) == [1, 0, 22, 0, 1]
      # M(Z/5,2), the suspension of Sage's MooreSpace(5): the Z/5 in H_2 shows
      # over F_5 in degrees 2 and 3 (universal coefficients), and not over F_7.
      KM = catalogue_space("M(Z/5,2)")
      @test dims(simplicial_cohomology_ring("M", KM, GF(5))) == [1, 0, 1, 1]
      @test dims(simplicial_cohomology_ring("M", KM, GF(7))) == [1, 0, 0, 0]
    end

    @testset "Manifold Page models" begin
      # CP^3 = Z[y]/(y^4) -- read out of arXiv:1012.3235
      XC = simplicial_cohomology_ring("CP^3", catalogue_space("CP^3"), ZZ)
      @test identify_ring(XC) == "Z[x2]/(x2^4),  |x2| = 2"
      # RP^5 over F_2
      @test identify_ring(simplicial_cohomology_ring("RP^5", catalogue_space("RP^5"), GF(2))) ==
            "F_2[x1]/(x1^6),  |x1| = 1"
      # lens space L(3,1): H^* = (Z, 0, Z/3, Z)
      @test groups(simplicial_cohomology_ring("L(3,1)", catalogue_space("L(3,1)"), ZZ)) ==
            ["Z", "0", "Z/3", "Z"]
      # the 3-torus
      @test groups(simplicial_cohomology_ring("T^3", catalogue_space("T^3"), ZZ)) ==
            ["Z", "Z^3", "Z^3", "Z"]
      # S^2 x S^2
      @test groups(simplicial_cohomology_ring("S^2xS^2", catalogue_space("S^2xS^2"), ZZ)) ==
            ["Z", "0", "Z^2", "0", "Z"]
      # The Wu manifold's file header disagrees with its own facet list: it
      # claims H_* = (Z,0,Z,Z,0,Z), but the facets give (Z,0,Z/2,0,0,Z).
      KW = catalogue_space("Wu manifold")
      @test f_vector(KW) == [13, 78, 286, 533, 468, 156]      # the published f-vector
      @test [elementary_divisors(Oscar.homology(KW, i)) for i in 0:5] ==
            [ZZRingElem[0], ZZRingElem[], ZZRingElem[2],
             ZZRingElem[], ZZRingElem[], ZZRingElem[0]]
      @test dims(simplicial_cohomology_ring("Wu", KW, GF(2))) == [1, 0, 1, 1, 0, 1]
      @test dims(simplicial_cohomology_ring("Wu", KW, GF(7))) == [1, 0, 0, 0, 0, 1]
    end
  end

end
