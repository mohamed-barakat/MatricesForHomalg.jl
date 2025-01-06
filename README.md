<!-- BEGIN HEADER -->
# MatricesForHomalg

### Matrices for the homalg project

<!--[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://homalg-project.github.io/MatricesForHomalg.jl/stable/)-->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://homalg-project.github.io/MatricesForHomalg.jl/dev/)
[![Build Status](https://github.com/homalg-project/MatricesForHomalg.jl/actions/workflows/Tests.yml/badge.svg?branch=main)](https://github.com/homalg-project/MatricesForHomalg.jl/actions/workflows/Tests.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/homalg-project/MatricesForHomalg.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/homalg-project/MatricesForHomalg.jl)

<!-- END HEADER -->

This package implements the homalg matrix interface in Julia. The
corresponding GAP package [MatricesForHomalg][MatricesForHomalg] is
used by the projects:

* [homalg project][homalg project],
* [CAP project][CAP project],
* [CategoricalTowers][CategoricalTowers],
* [HigherHomologicalAlgebra][HigherHomologicalAlgebra].

We are gradually creating Julia-versions of the various GAP packages
in these projects and many of them will rely on this package.

The current version only supports matrices over $\mathbb{Z}$ and
$\mathbb{Q}$. Matrices over [Bézout domains][Bézout domains] will be
supported next. Future versions should cover the full functionality of
the GAP package [MatricesForHomalg][MatricesForHomalg] by supporting
matrices over so-called computable rings [[BR][BR], [BL][BL]].

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).

[MatricesForHomalg]: https://homalg-project.github.io/pkg/MatricesForHomalg/
[CAP project]: https://homalg-project.github.io/prj/CAP_project/
[homalg project]: https://homalg-project.github.io/prj/homalg_project/
[CategoricalTowers]: https://homalg-project.github.io/prj/CategoricalTowers/
[HigherHomologicalAlgebra]: https://homalg-project.github.io/prj/HigherHomologicalAlgebra/
[Bézout domains]: https://ncatlab.org/nlab/show/B%C3%A9zout+domain
[BR]: https://arxiv.org/abs/math/0701146
[BL]: https://arxiv.org/abs/1003.1943
