# WIP 2D sampling

This is a small application for exploring different 2D sampling strategies. It is mostly exploratory
and implementations here have not been optimized for performance. Currently Poisson Disk Sampling
is doing a bunch of avoidable mallocs.

Currently implemented:

1. Uniform random sampling
2. Stratified sampling
3. Poisson Disk Sampling [Fast Poisson Disk Sampling in Arbitrary Dimensions, Bridson]

## TODO

- Visualization: Allow switching between different samplings during runtime.
- Low discrepancy sequences
- Higher dimensional sampling
- Variants of sampling according to probability map

## Compilation and dependencies

This is just a single file with dependencies for visualization. You can compile it
by simply:

~~~c
gcc sampling.c -o sampling -lglfw3 -lopengl32 -lglew32 -lnanovg
~~~

Dependency list:

- glfw
- glew
- nanovg