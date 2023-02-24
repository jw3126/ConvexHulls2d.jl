# ConvexHulls2d

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jw3126.github.io/ConvexHulls2d.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jw3126.github.io/ConvexHulls2d.jl/dev/)
[![Build Status](https://github.com/jw3126/ConvexHulls2d.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jw3126/ConvexHulls2d.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jw3126/ConvexHulls2d.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jw3126/ConvexHulls2d.jl)

Zero dependency convex hulls in 2d.

# Alternatives
There are many packages in julia capapable of computing convex hulls. Some examples are:
* [LazySet](https://github.com/JuliaReach/LazySets.jl)
* [QHull](https://github.com/JuliaPolyhedra/QHull.jl)
* [MiniQHull](https://github.com/gridap/MiniQhull.jl)
* [LibGEOS](https://github.com/JuliaGeo/LibGEOS.jl)

All of these packages offer much more than this package, however they are more heavy depependencies.
Also this package is faster then most of them for small to medium numbers of points.
