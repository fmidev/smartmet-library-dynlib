# smartmet-library-dynlib — design notes

Written overnight 2026-04-17 → 2026-04-18 for morning review.

## Status at hand-off

### smartmet-library-dynlib (new repo, `~/hub/smartmet-library-dynlib/`)

- Compiles cleanly with `make`; 10 unit tests pass with `make test`.
- Staged locally at `$HOME/smartmet-dynlib-staging/` (lib64 + include/smartmet/dynlib).
- Vendored dynlib commit `cc4fc9e` (Clemens Spensberger, UiB, MIT) —
  only the Fortran files needed for detection, not the Python layer.
- Fortran is compiled with `gfortran -O2 -fPIC -fno-range-check`.
- `smartmet-library-dynlib.spec` lists `BuildRequires: gcc-gfortran,
  lapack-devel` and `Requires: blas, lapack, libgfortran` (updated after
  the user installed lapack-devel and noticed blas came with it).

### WMS plugin (`~/hub/brainstorm/plugins/wms/`)

- **New**: `GridFrontSource` — pulls scalar field + u/v from the
  querydata engine, calls `Fmi::Dynlib::detectFrontsMaxGrad`, returns
  `FrontCurve` curves in WGS84. Plugs into existing `WeatherFrontsLayer`
  via the new `"front_source": "grid"` option.
- **New**: `WeatherObjectsLayer` (`"layer_type": "weather_objects"`) —
  the single-layer design the user requested. Accepts an `objects`
  array where each entry selects one detector type. Supported:
  `fronts`, `fronts_maxcurv`, `jet_axes`, `trough_axes`,
  `convergence_lines`, `deformation_lines`, `vorticity_lines`.
- LayerFactory wires `weather_objects` as a new layer type.
- Makefile has a new `DYNLIB_STAGE` variable defaulting to
  `$HOME/smartmet-dynlib-staging`, with `-I`, `-L`, and `-Wl,-rpath`
  pointing there so the local build links without installing the RPM.
- My new files (`WeatherObjectsLayer.o`, `GridFrontSource.o`,
  `LayerFactory.o`, `WeatherFrontsLayer.o`) all compile cleanly.

### Not fixed

- A full `make wms.so` build fails in `Plugin.o` due to a *pre-existing*
  macgyver↔spine template version skew on the installed headers
  (`Fmi::Cache::Cache` expects different template parameters than
  `SmartMetCache.h` provides). This is unrelated to dynlib integration
  and was broken before I touched the repo. The user's local
  environment needs either an updated `smartmet-library-spine`
  install or a rebuild-and-install of the sibling repos in `~/hub`.
- Full WMS integration tests (in `test/input/*.get`) against live
  grid data were not added — they require running backends and a
  well-known producer. Unit tests for dynlib itself (10 passing)
  exercise the C → Fortran boundary, which is where the integration
  risk lives.

## Design choices I made without asking

1. **Vendor the Fortran, don't link upstream as RPM.** dynlib isn't
   packaged anywhere. Vendored as a subtree-style copy under
   `third_party/dynlib/` with upstream LICENSE preserved and a
   `UPSTREAM` file documenting the commit and refresh procedure.

2. **One C ABI per detector.** Alternatives considered: one generic
   `dynlib_detect(kind_enum, ...)` dispatch. I went with per-detector
   entry points (`dynlib_detect_fronts_maxgrad`, `_maxcurv`,
   `_jet_axes`, `_lines`, `_block_indicator`) because each has
   different inputs (scalar field vs wind-only) and different output
   shapes (three-type typed vs single-list vs scalar field).

3. **C++ API uses `Fmi::Matrix<double>`, not `NFmiDataMatrix<float>`.**
   The user flagged `NFmiDataMatrix` as a `vector<vector<>>` with bad
   cache locality and no contiguous pointer. `Fmi::Matrix<double>` is
   contiguous, `i` fastest, and its `&m(0,0)` pointer hands directly
   to the Fortran shim with no copy. The shim transposes internally
   because upstream dynlib declares fields as `(nz, ny, nx)` in
   Fortran column-major, which is `(x-slowest, y-fastest)` in C
   memory — the transpose of `Fmi::Matrix`.

4. **Shim always initialises outputs to NaN and resets module-level
   config on each call.** dynlib's Fortran subroutines declare
   outputs as `intent(out)`, and `line_locate` only writes the valid
   prefix — unused slots are uninitialised stack memory. The shim
   NaN-fills outputs and the C decoder stops at the first NaN /
   non-monotone offset. Config variables (`nsmooth`, `frint_thres`,
   `frspd_thres`) live at module scope in upstream; the shim resets
   them to the documented defaults on every call so a previous
   caller's override cannot leak into the next detection.

5. **Trough detection = vorticity lines.** dynlib does not ship a
   dedicated trough detector. Spensberger's group treats trough
   axes as maxima of cyclonic relative vorticity, so
   `detectTroughAxes` is an inline alias for `detectVorticityLines`.
   Flag for morning: if you want signed trough detection (sign the
   vorticity based on hemisphere), add it at the C++ layer — the
   Fortran `vorline` subroutine is already sign-agnostic.

6. **`WeatherObjectsLayer` is a new layer, not an extension of
   `WeatherFrontsLayer`.** The user's suggestion was "one layer
   listing all the objects". WeatherFrontsLayer has a mature SVG
   path + glyph rendering pipeline with typed styling for the six
   classical front/ridge/trough categories; mixing in jet axes and
   vorticity lines would have complicated that interface. The split
   is: fronts keep their existing rendering (now reachable via
   `front_source: "grid"`); jets/troughs/lines render as simple
   `<path>` elements with CSS classes in the new layer. A user can
   use both layers in the same view for a combined display.

7. **Edge-row grid spacing mirrors the nearest interior.** The
   upstream `dx(j,i)`/`dy(j,i)` convention is the distance between
   `(j, i-1)` and `(j, i+1)`. At `i=0` and `i=nx-1` that would
   require out-of-grid samples. I use the same value as the nearest
   interior row. This keeps `dx`/`dy` well-defined everywhere and is
   accurate to first order at the boundary. An alternative is to
   throw at the boundary, but that would make detection fail for
   any grid that is not artificially padded.

8. **Smoothing / threshold overrides via JSON, not layer-level
   metaparameters.** `ComputedFields` uses a `"tfp"` block on the
   isoband layer to scope TFP options; I put detection options
   directly in each `objects[]` entry of `weather_objects`. This
   keeps the layer's schema flat and gives each object independent
   tuning.

## Follow-up additions (2026-04-18 afternoon)

Added C ABI + C++ wrappers for the four detectors that were deferred
at hand-off, plus a jet-axis variant:

- **`detectJetAxesFFThres`** — `jetaxis_ff_thres` variant with the
  fixed 30 m/s wind-speed mask. Exposed via a new `variant` argument
  on `dynlib_detect_jet_axes`; C++ has both `detectJetAxes` (default)
  and `detectJetAxesFFThres`.
- **`detectRossbyWaveBreakingGradRev`** — wraps `rwb_by_grad_rev`.
  Returns a `RwbResult` struct with all seven per-cell fields
  (anticyc/cyc flag + |grad| + |dPV/dy| + tested mask) as
  `Fmi::Matrix<double>`. Takes a user mask; passes `nullptr` to test
  everywhere.
- **`detectCyclonesByContour`** — wraps `cyclone_by_contour_fortran`.
  The C++ side does the required MSL sort + (i, j) unpacking; caller
  just supplies MSL, orography, lon/lat. Cyclone meta saved/restored
  around the call so a previous call's overrides cannot leak state.
- **`detectPrecipitationBlobs`** — wraps `blobs_fortran`. Same sort
  protocol with a sign flip (upstream triggers on `val < 0`). Returns
  a `BlobResult` with mask + metadata, precip values restored to
  positive sign on output.

Six new unit tests added. All 17 tests pass.

**One gotcha, still unresolved:**
`frontline_at_maxgrad` + the upstream `line_locate` subroutine declare
their output buffers as `intent(out)` but only write the valid prefix.
Under gfortran -O2 the shim's pre-call NaN fill can be elided as a
dead store, so unused slots carry stack leftovers from neighbouring
frames. With the small original test binary this never surfaced; with
more test cases the stack layout shifts and the front detector started
returning 0 lines for inputs that should return 78. The Makefile now
compiles the vendored Fortran with `-O0 -g` to keep output
deterministic. A proper fix would be (a) refactoring the shim to not
rely on `intent(out)` pre-fill semantics, (b) patching upstream
`line_locate` to declare its buffers `intent(inout)`, or (c) porting
the numerics to C++ as originally considered.

## Deferred
- **Caching.** Each `WeatherObjectsLayer::generate` re-fetches grid
  data and re-runs detection every render. For animations this is
  probably prohibitively slow. Consider a `Fmi::Cache::Cache` keyed
  on `(producer, time, object_spec_hash)`.
- **Coordinate system edge cases.** The current grid-index → lon/lat
  conversion is bilinear on the coordinate matrix. Fine for
  continuous grids, unsafe across the antimeridian or across polar
  singularities. If you expose this for global products,
  `gis::Interrupt` handling is needed.
- **Grid-engine backend.** I only wired the querydata engine. The
  grid engine (`Engine::Grid`) uses a different API. An analogous
  data-fetch path is needed if dynlib should work against
  GRIB/NetCDF-backed producers. The `GridFrontSource` + layer
  factoring makes this easy to add without touching the C++/Fortran
  boundary.
- **SVG styling.** The `weather_objects` layer emits
  `<path class="jet-axis"/>` etc. with no arrow heads or direction
  decoration. Jets typically get an arrowhead at the downstream end.
  This is SVG work; no detection changes needed.
- **ECMWF / ERA5 tuning.** Upstream thresholds (`frint_thres =
  -2e-11`, `frspd_thres = 1.5 m/s`, `nsmooth = 2`) were tuned for
  ERA-Interim. Harmonie/MEPS will probably need different values.
  The JSON schema already exposes the overrides.

## Alternatives I did not take (listed for discussion)

- **Reimplement the numerics in C++ instead of linking Fortran.** I
  chose the link route because the audit showed kernels are clean
  (no netCDF / no Python coupling at the Fortran level) and MIT-
  licensed, so link is cheaper than a port. A pure-C++ port is still
  attractive if you want SmartMet to control the algorithm code path
  and avoid the gfortran dependency. Cost: weeks.
- **Embed CPython and call dynlib as a Python library.** Fastest to
  prototype but unsuitable for a long-running multi-threaded plugin
  (GIL, interpreter lifetime, NumPy in the closed tree). Rejected.
- **One unified layer that also handles blocking / cyclone masks.**
  Would require merging line-feature rendering with isoband-style
  contouring of a scalar field in the same layer. Probably better as
  a dedicated `blocking` object that emits computed fields consumable
  by the existing `IsobandLayer` — more work but cleaner separation.
  Deferred.
- **Migrate `WeatherFrontsLayer` and `ComputedFields::computeTFP` to
  `Fmi::Matrix<double>` too.** Consistent, avoids the
  `NFmiDataMatrix<float>` per-row allocations inside the plugin.
  Non-trivial touch across IsobandLayer/IsolineLayer. Deferred.

## How to validate in the morning

```bash
# 1. Library
cd ~/hub/smartmet-library-dynlib
make && make test
# Expect: 10 test cases, all pass.

# 2. WMS individual files compile
cd ~/hub/brainstorm/plugins/wms
make obj/GridFrontSource.o obj/WeatherObjectsLayer.o \
     obj/WeatherFrontsLayer.o obj/LayerFactory.o
# Expect: clean builds.

# 3. Full WMS build
make
# Will currently fail in Plugin.o with a macgyver↔spine template
# version skew unrelated to this work. Fix by rebuilding sibling
# libraries in ~/hub and reinstalling, then rerun.

# 4. Example product JSONs
ls ~/hub/smartmet-library-dynlib/examples/wms/
```
