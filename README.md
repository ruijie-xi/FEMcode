# FEMcode

MATLAB finite element routines for research and learning on two-dimensional triangular meshes.

[简体中文](README.zh-CN.md)

FEMcode provides mesh topology, finite element spaces, basis evaluation, sparse assembly, boundary treatment, interpolation, error measurement, and plotting. It is a research codebase with a runnable Poisson example and focused regression checks; comprehensive validation of all element types is still incomplete.

## Requirements and quick start

- MATLAB. A minimum supported release has not been established; GNU Octave compatibility has not been verified.
- Optional: external DistMesh for `mesh_circle_triangle` and `mesh_square_triangle_unstructured`. It is not bundled and is not needed for the structured Poisson example.

Open MATLAB in the repository root and run:

```matlab
addpath(pwd);
addpath(fullfile(pwd,'poisson'));
solver_2D_Poisson_d;
```

The example uses P1 elements on a structured triangular mesh with 100 subdivisions per coordinate direction. It solves

$$
-\Delta u=-2e^{x+y}\quad\text{in }(-1,1)^2,
\qquad u=e^{x+y}\quad\text{on the boundary}.
$$

The exact solution is `exp(x+y)`. The program prints mesh size, elapsed time, and L2, H1, and sampled maximum errors. Problem data are in `poisson/fun_*.m`.

## Implemented spaces

These identifiers have implementation branches; inclusion does not imply comprehensive numerical verification.

| Family | Identifiers |
| --- | --- |
| Piecewise polynomial and enriched scalar spaces | `P0`, `P1`, `P2`, `bubbleP1` |
| Discontinuous spaces | `P1dc`, `DGP1`, `DG-P2-quad` |
| Crouzeix–Raviart | `CR` |
| Raviart–Thomas | `RT0` |
| Bernardi–Raugel and Mardal–Tai–Winther | `BR`, `MTW` |
| Nédélec | `Ned1-1`, `Ned1-2`, `Ned2-1` |
| Research/legacy variants | `CR-P0`, `CR-RT0`, `BR0`, `BR-RT0`, `RT0-MTW` |

Scalar spaces can be stacked through `dim`. Intrinsically vector-valued spaces such as RT0 and Nédélec use `dim = 2` and share degrees of freedom across vector-basis components. Consult the implementation before using research variants. Complete quadrilateral and three-dimensional workflows are not provided.

## Core workflow

```matlab
[T,P] = mesh_square_triangle(0,1,0,1,16,16);
mesh = mesh_topo(T,P);
space = mesh_FE('P1',1,mesh,3,0,[]);
A = assemble_matrix_2D(@(x) ones(1,size(x,2)),[1,1], ...
    mesh,space,space,[1,1,1,0,1,0;1,1,0,1,0,1],3);
```

This constructs a mesh, a scalar P1 space, and a Laplace stiffness matrix. The Poisson example shows load assembly, boundary treatment, and solution.

### Mesh and degrees of freedom

`P` is a `2 x N_node` coordinate array; `T` is a `3 x N_elem` vertex-index array. `mesh_topo` makes triangle orientation consistent and builds connectivity.

| Field | Meaning |
| --- | --- |
| `mesh.E` | Two adjacent element slots followed by two endpoint indices for each oriented edge |
| `mesh.TE` | Signed global edge indices; local order `(1,2)`, `(2,3)`, `(3,1)` |
| `mesh.ET` | Local edge numbers in the adjacent elements, corresponding to `E(1:2,:)` |
| `mesh.tau`, `mesh.N` | Unit tangent and its clockwise-rotated unit normal |
| `mesh.s`, `mesh.hmax` | Edge lengths and maximum triangle edge length |

On boundary edges, `E(2,:)` holds the element and `E(1,:)` holds a nonpositive boundary label. `mesh_topo(T,P)` labels all boundary edges zero. Optional functions in `mesh_topo(T,P,f0,f1,...)` assign labels `0,-1,...` when both endpoints satisfy the function's zero test. Unmatched boundary edges remain zero. Labels select edges; they do not impose boundary conditions by themselves.

The space constructor takes six arguments:

```matlab
space = mesh_FE(basis_type,dim,mesh,Gauss_type_2D,Gauss_type_1D,bndy_number);
```

- `Gauss_type_2D`: triangle quadrature point count (`1`, `3`, `4`, or `9`), also used to cache basis values and first derivatives.
- `Gauss_type_1D`: edge quadrature point count (`1`–`5`); `0` skips edge caching.
- `bndy_number`: a boundary label for edge caching; `[]` includes all edges.

`space.N_node` means **total degrees of freedom**, despite its name. `N_lb` is the local basis count; `T` maps local to global degrees of freedom; `basis_type` stores the character identifier. `space.P` contains coordinates or element-dependent degree-of-freedom metadata. Vector-basis components may share global indices. Use `get_Dbndynodes(mesh,space,label)` to retrieve boundary degrees of freedom; not every space provides a `Dbndynodes` field.

### Assembly and cache conventions

`assemble_matrix_2D` assembles terms of the form

$$
c\int_\Omega f\,(\partial_x^i\partial_y^j u_k)
(\partial_x^m\partial_y^n v_l)\,dx.
$$

Each row of `mat_info` is `[k,l,i,j,m,n]`; `coe_num` supplies the scalar coefficients. Trial and test spaces are separate arguments. `assemble_vector_2D` assembles loads. The `*_FE_*` assemblers accept finite element coefficient functions; see their source headers for argument conventions.

`treat_Dirichlet(A,b,mesh,space,label,boundary_fun)` replaces matrix rows to prescribe degrees of freedom. This generally does not preserve symmetry. The Poisson example uses MATLAB's direct backslash solver.

**Quadrature must match the basis cache.** Assemblers and cached error routines must use the triangle rule used to construct the space. Before evaluating an existing solution with a different error rule, rebuild the cache:

```matlab
[space.basis,space.basis_dx,space.basis_dy] = ...
    generate_basis_data_triangle(mesh,space,9);
% Subsequent cached triangle evaluations must use rule 9.
```

## Errors and verification

- `compute_error`: `L2`, `H1`, `H10` (H1 seminorm), `L_inf`, and `Hdiv`. The exact-solution callback is `u_fun(points,dx,dy)`, with one row per component.
- `compute_norm`: `L2`, `H1`, `H1-semi`, `L_inf`, and `Hdiv-semi`.
- `Hdiv` error is the full norm `sqrt(||e||_L2^2 + ||div(e)||_L2^2)`; `Hdiv-semi` is `||div(u_h)||_L2`. Both require two-component fields. Elementwise derivatives give broken quantities for nonconforming fields.
- `L_inf` samples quadrature points; it is not an exact continuous maximum.

Run the focused regression checks from the repository root:

```matlab
addpath(pwd);
addpath(fullfile(pwd,'test'));
test_Hdiv;
```

These check divergence cancellation, nonzero divergence, exact interpolation, and constant-error L2 contributions with vector P1 fields. Other scripts under `test/` are exploratory, not a comprehensive regression suite.

A separate local refinement experiment for the Poisson problem above, retaining three-point assembly and using nine-point error integration, gave:

| Subdivisions per direction | L2 error | L2 order | H1 error | H1 order |
| ---: | ---: | ---: | ---: | ---: |
| 8 | 2.08012e-02 | — | 3.71880e-01 | — |
| 16 | 5.18008e-03 | 2.0056 | 1.85299e-01 | 1.0050 |
| 32 | 1.29373e-03 | 2.0014 | 9.25683e-02 | 1.0013 |
| 64 | 3.23353e-04 | 2.0004 | 4.62740e-02 | 1.0003 |
| 128 | 8.08331e-05 | 2.0001 | 2.31357e-02 | 1.0001 |

Orders are `log2(E_N/E_2N)`. This experiment is not yet shipped as an automated convergence test. The default example uses three-point error integration, so its error values differ. These results validate this smooth Poisson case, not every element or boundary treatment. Stokes, magnetic diffusion, and ALE examples are not included; no CI workflow is configured.

## Source map

| Files | Purpose |
| --- | --- |
| `mesh_*.m`, `orientation_consistence.m` | Mesh generation, topology, and spaces |
| `basis_*.m`, `generate_basis_data_*.m` | Basis evaluation and caching |
| `generate_Gauss_*.m`, `local_*.m`, `assemble_*.m` | Quadrature and assembly |
| `compute_dof.m`, `interpolate.m`, `FE_function_*.m` | Degrees of freedom, interpolation, and evaluation |
| `treat_Dirichlet*.m`, `get_Dbndynodes.m` | Boundary treatment |
| `compute_error.m`, `compute_norm.m`, `plot_*.m` | Diagnostics and plotting |
| `poisson/`, `test/` | Example problem and checks |

## Origin and license

The initial code grew out of Xiaoming He's finite element programming course at Missouri University of Science and Technology, followed by modifications and extensions. The original course reference is retained in the [Chinese README](README.zh-CN.md#来源与许可证).

The repository includes the [GNU General Public License, version 3](LICENSE.md).
