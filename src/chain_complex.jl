#=
The simplicial chain complex C_*(K; R).

OSCAR ships the *cochain* complex (experimental/DoubleAndHyperComplexes/src/
Objects/simplicial_complex.jl) because that is what carries the cup product.
This is its mirror image, and exists because `Oscar.homology(K, i)` on a
`SimplicialComplex` takes no coefficient ring -- it is Polymake, integral only --
so homology over QQ or GF(p) has nowhere else to come from.

Three differences from OSCAR's cochain version, and no others:

  * direction `:chain`, so maps run i -> i-1 rather than i -> i+1;
  * the map is the boundary matrix *untransposed* (the cochain version
    transposes it to get the coboundary);
  * maps exist for 1 <= i <= dim(K) rather than 0 <= i < dim(K).

The chain factory is identical: C_i is free on the i-faces. Should this ever be
upstreamed, the diff against OSCAR's file is those three lines.

Note the `Oscar.` qualifications throughout: `AbsHyperComplex`, `HyperComplex`,
`HyperComplexChainFactory`, `HyperComplexMapFactory`, `can_compute`,
`chain_factory` and `underlying_complex` are not in Oscar's export list.
=#

### Production of the chains: C_i is free on the i-faces.
struct SimplicialChainFactory{ChainType} <: Oscar.HyperComplexChainFactory{ChainType}
  R::Ring
  K::SimplicialComplex

  function SimplicialChainFactory(R::Ring, K::SimplicialComplex)
    return new{FreeMod{elem_type(R)}}(R, K)
  end
end

function (fac::SimplicialChainFactory)(self::Oscar.AbsHyperComplex, i::Tuple)
  return FreeMod(fac.R, f_vector(fac.K)[only(i) + 1])
end

function Oscar.can_compute(fac::SimplicialChainFactory, self::Oscar.AbsHyperComplex, i::Tuple)
  return 0 <= only(i) <= dim(fac.K)
end

### Production of the morphisms: d_i : C_i -> C_{i-1}.
struct SimplicialChainMapFactory{MorphismType} <: Oscar.HyperComplexMapFactory{MorphismType}
  function SimplicialChainMapFactory()
    return new{FreeModuleHom}()
  end
end

function (::SimplicialChainMapFactory)(self::Oscar.AbsHyperComplex, p::Int, I::Tuple)
  fac = Oscar.chain_factory(self)
  i = only(I)
  return hom(self[i], self[i - 1],
             matrix(fac.R, Polymake.topaz.boundary_matrix(fac.K.pm_simplicialcomplex, i)))
end

function Oscar.can_compute(::SimplicialChainMapFactory, self::Oscar.AbsHyperComplex,
                           p::Int, I::Tuple)
  p == 1 || return false
  return 1 <= only(I) <= dim(Oscar.chain_factory(self).K)
end

@doc raw"""
    SimplicialChainComplex(R::Ring, K::SimplicialComplex)

The simplicial chain complex of `K` with coefficients in `R`, as a
1-dimensional hypercomplex concentrated in degrees `0 .. dim(K)`.

`Oscar.homology(C, i)` then gives `H_i(K; R)` as a subquotient. Beware that its
generators are not a minimal generating set -- see `simplicial_homology`.
"""
@attributes mutable struct SimplicialChainComplex{ChainType, MorphismType} <:
                           Oscar.AbsHyperComplex{ChainType, MorphismType}
  internal_complex::Oscar.HyperComplex{ChainType, MorphismType}

  function SimplicialChainComplex(R::Ring, K::SimplicialComplex)
    internal_complex = Oscar.HyperComplex(
      1, SimplicialChainFactory(R, K), SimplicialChainMapFactory(), [:chain];
      upper_bounds = Union{Int, Nothing}[dim(K)],
      lower_bounds = Union{Int, Nothing}[0])
    return new{FreeMod{elem_type(R)}, FreeModuleHom}(internal_complex)
  end
end

Oscar.underlying_complex(C::SimplicialChainComplex) = C.internal_complex

"The simplicial complex a chain complex was built from."
Oscar.simplicial_complex(C::SimplicialChainComplex) = Oscar.chain_factory(C).K

"The ring a chain complex is defined over."
Oscar.base_ring(C::SimplicialChainComplex) = Oscar.chain_factory(C).R
