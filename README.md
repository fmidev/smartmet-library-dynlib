# smartmet-library-dynlib

C++ bindings for a subset of [dynlib](https://git.app.uib.no/Clemens.Spensberger/dynlib),
Clemens Spensberger's dynamic-meteorology feature-detection library
(University of Bergen, MIT-licensed).

The upstream library is a Fortran core with Python (f2py) bindings. This
package vendors the numerical kernels, adds an `ISO_C_BINDING` shim, and
exposes the detection entry points through a small C++ API in the `Fmi::Dynlib`
namespace.

## Current scope

- Front detection, Jenkner et al. (2010) max-gradient method: `detectFrontsMaxGrad`

Further detection algorithms already present in the vendored Fortran
(front max-curv, jet axis, convergence / deformation / vorticity lines,
cyclone-by-contour, blocking indicator, Rossby wave breaking, precipitation
blobs) will be exposed incrementally via the same pattern: add a bind(c)
subroutine in `dynlib_wrapper.f90`, a declaration in `DynlibC.h`, a C++
overload in `Dynlib.h`.

## Build

Requires `gcc-gfortran`, `lapack-devel`, and `smartmet-library-macgyver-devel`.

```bash
make                # libsmartmet-dynlib.so
make test           # run unit tests against the local build
make install        # PREFIX=/usr by default
make rpm
```

## Third-party sources

`third_party/dynlib/` contains a vendored subset of the upstream dynlib
Fortran sources (commit `cc4fc9e`, fetched 2026-04-17). Files are
unmodified; see `third_party/dynlib/UPSTREAM` for the refresh procedure.
The upstream MIT licence is preserved at `third_party/dynlib/LICENSE`.

## Runtime dependencies

- `liblapack.so.3` — called from `diag.f90` (3×3 symmetric eigendecomp)
- `libgfortran.so.5` — Fortran runtime

## Array-order convention

All C++ inputs are `Fmi::Matrix<double>(nx, ny)` in standard row-major
order (x-index innermost). The Fortran shim transposes internally
because upstream dynlib subroutines declare fields as `dat(nz, ny, nx)`.
The transpose cost is negligible relative to the detection work.
