# PeriodicGraphs.jl

A Julia package for the manipulation of periodic graphs.

## Usages

`PeriodicGraphs.jl` defines types and algorithms useful for the manipulation of `N`-periodic undirected graphs, with `N` known at compile-time. It extends `Graphs.jl` by providing the new `PeriodicGraph` type which implements the `AbstractGraph` API. It also provides an API for the manipulation of symmetries.

This package serves as support for the `PeriodicGraphEmbeddings.jl` library, which focuses on euclidean embeddings of periodic graphs, and for `CrystalNets.jl` for the manipulation and identification of crystal nets, which are 2 or 3-periodic graphs.

## Package installation

The installation follows the usual procedure. Start by downloading and installing [Julia](https://julialang.org/) (v1.6 or higher for `PeriodicGraphs.jl`). Then, either

- open the Julia REPL and enter the package manager by typing `]`, then install `PeriodicGraphs.jl` by entering:
  ```julia
  pkg> add PeriodicGraphs
  ```
- alternatively, you can do it from a shell by executing:
  ```bash
  julia -e 'import Pkg; Pkg.add("PeriodicGraphs")
  ```

To use the package, open a REPL and enter

```julia
julia> using PeriodicGraphs
```
