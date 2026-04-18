# smartmet-library-dynlib

C++ bindings for a subset of [dynlib](https://git.app.uib.no/Clemens.Spensberger/dynlib),
Clemens Spensberger's dynamic-meteorology feature-detection library
(University of Bergen, MIT-licensed).

The upstream library is a Fortran core with Python (f2py) bindings. This
package vendors the numerical kernels, adds an `ISO_C_BINDING` shim, and
exposes the detection entry points through a small C++ API in the `Fmi::Dynlib`
namespace.

## Current scope

All entry points live in `Fmi::Dynlib`.

| Detector                           | Function(s)                                              |
|------------------------------------|----------------------------------------------------------|
| Jet axes (zero-shear)              | `detectJetAxes`, `detectJetAxesFFThres` (30 m/s cutoff)  |
| Blocking indicator (Masato 2012)   | `blockingIndicator`                                      |
| Trough axes / vorticity lines      | `detectTroughAxes`, `detectVorticityLines`               |
| Convergence / deformation lines    | `detectConvergenceLines`, `detectDeformationLines`       |
| Cyclones (Wernli & Schwierz)       | `detectCyclonesByContour`                                |
| Rossby wave breaking               | `detectRossbyWaveBreakingGradRev`                        |
| Precipitation blobs                | `detectPrecipitationBlobs`                               |
| Fronts (gradient / curvature)      | `detectFrontsMaxGrad`, `detectFrontsMaxCurv`             |

### Caveat on front detection

`detectFrontsMaxGrad` and `detectFrontsMaxCurv` produce gradient-based
line output that does **not** match meteorologist-drawn frontal
analysis closely enough for operational chart-quality use. They remain
in the library as inputs to downstream diagnostics and for research
workflows. Do not advertise them as a ready-to-use product. See
`DESIGN_NOTES.md` for the validation experiment and the conclusion
that chart-quality fronts should come from a precomputed external
source (ML-trained or a tuned classical pipeline with post-processing).

The other detectors above use objective physical criteria and produce
output that is meaningful to serve directly.

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
