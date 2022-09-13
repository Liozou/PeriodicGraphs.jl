# PeriodicGraphs

[![Build Status](https://ci.appveyor.com/api/projects/status/github/Liozou/PeriodicGraphs.jl?svg=true)](https://ci.appveyor.com/project/Liozou/PeriodicGraphs-jl)
[![codecov](https://codecov.io/gh/Liozou/PeriodicGraphs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/Liozou/PeriodicGraphs.jl)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://liozou.github.io/PeriodicGraphs.jl/)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This module provides the new `PeriodicGraph{N}` type to manipulate `N`-dimensional periodic
graphs, extending the `AbstractGraph` API from [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl/).

Various optimized algorithms are implemented as well, including:

- Neighborhood exploration and coordination sequence computation
- Separation by dimension of the connected components
- Simple ring statistics

and more, see [the documentation](https://liozou.github.io/PeriodicGraphs.jl/).

See also:

- [PeriodicGraphEmbeddings.jl](https://github.com/Liozou/PeriodicGraphEmbeddings.jl)
  for a dependent package specialized on the manipulation of periodic graph embeddings,
  including symmetry computations for 3D graphs
- [CrystalNets.jl](https://github.com/coudertlab/CrystalNets.jl) for a dependent package
  specialized on crystal nets, including a periodic graph canonization algorithm for 3D,
  2D and 1D periodic graphs.
