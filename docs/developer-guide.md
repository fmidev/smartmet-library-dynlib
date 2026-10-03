# dynlib developer guide

This guide is for developers who change `smartmet-library-dynlib`. The library wraps a
vendored subset of Clemens Spensberger's dynlib (Fortran) behind a small C++ API, for
detecting meteorological features from gridded fields: jet axes, troughs, convergence and
deformation lines, fronts, cyclones, Rossby wave breaking, blocking and precipitation blobs.
The WMS plugin's `WeatherObjectsLayer` uses it.

[README.md](../README.md) lists the detectors and their caveats, and
[DESIGN_NOTES.md](../DESIGN_NOTES.md) records the design decisions and the front-detection
experiment.

## Contents

1. [Building and testing](#1-building-and-testing)
2. [Layers of the code](#2-layers-of-the-code)
3. [Arrays and grid spacing](#3-arrays-and-grid-spacing)
4. [Configuration and threads](#4-configuration-and-threads)
5. [Updating the Fortran code](#5-updating-the-fortran-code)
6. [Known pitfalls](#6-known-pitfalls)

---

## 1. Building and testing

```bash
make           # needs gcc-gfortran and lapack-devel
make test
```

`test/DynlibTest.cpp` (Boost.Test) runs each detector on synthetic fields and checks the C++ ↔ Fortran
boundary (shapes, output arrays, missing values). Check real output through the WMS
`WeatherObjectsLayer` products.

## 2. Layers of the code

| Layer | Files | Role |
|-------|-------|------|
| C++ API | `dynlib/Dynlib.h`, `Dynlib.cpp` | `Fmi::Dynlib::detect…()` functions taking `Fmi::Matrix<double>` and returning lines, cyclones, blobs or fields. |
| C ABI | `dynlib/DynlibC.h` | One `extern "C"` entry point per detector. |
| Fortran shim | `third_party/dynlib/dynlib_wrapper.f90` | `ISO_C_BINDING` wrappers that transpose the arrays, set the configuration, call the upstream routines and fill the outputs. |
| Upstream | the other `third_party/dynlib/*.f90` | Unmodified upstream sources. |

## 3. Arrays and grid spacing

* C++ inputs are `Fmi::Matrix<double>(nx, ny)` with x varying fastest. Upstream declares
  fields as `(nz, ny, nx)`, so the shim transposes them.
* The detectors take `dx`, `dy` as **double** grid spacings in metres, per grid point: the
  distance between the neighbours `i-1` and `i+1`, as dynlib's centred differences expect.
  `latLonDoubleGridSpacing()` computes them for regular lat/lon grids.
* Line outputs come back as an offset array and a point array, in fixed-size buffers that
  the shim initialises to NaN. The C++ side decodes them (`decodeLines()`), stopping at the
  first NaN or invalid offset, so the unused part of a buffer is never read as data.

## 4. Configuration and threads

The upstream detectors read their thresholds and smoothing from Fortran **module
variables** (`config.f90`). The shim's `apply_config()` resets them to the defaults and
applies the caller's overrides at the start of every call.

Module variables are process-wide, so the C++ API serialises all calls into Fortran with
one mutex: each detection runs with the options it asked for, and concurrent detections
run one at a time. Code that calls the C functions of `DynlibC.h` directly must serialise
them itself.

## 5. Updating the Fortran code

`third_party/dynlib/UPSTREAM` records the vendored commit, the subset of files and the
refresh procedure. Keep the upstream files unmodified and make all adaptations in
`dynlib_wrapper.f90`, so that a refresh stays a copy. After a refresh, rebuild and run the
tests, since upstream may change the argument lists the shim relies on.

## 6. Known pitfalls

* **Detections do not run in parallel** (§4), and direct callers of `DynlibC.h` must
  serialise their calls.
* **Front detection is not chart quality** (see the README); do not offer it as a finished
  product.
* **Grid spacings are double spacings in metres** (§3). Passing the plain spacing doubles
  every gradient, and passing degrees makes the thresholds meaningless.
* **The library needs the Fortran runtime and LAPACK** at run time (`libgfortran`,
  `liblapack`).
