#pragma once

#include <macgyver/Matrix.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace Fmi
{
namespace Dynlib
{

// ---------------------------------------------------------------------------
// Common types
// ---------------------------------------------------------------------------

enum class FrontType : uint8_t
{
  Cold = 0,
  Warm = 1,
  Stationary = 2
};

// A single point on a detected feature. Coordinates are in grid-index
// space (j = row, i = column), both fractional because detection uses
// interpolated zero crossings. `value` is the native scalar sampled
// along the feature (e.g. TFP for fronts, wind speed for jets).
struct FeaturePoint
{
  double j;
  double i;
  double value;
};

// A single-type line feature (e.g. one jet axis, one convergence line,
// one trough axis).
struct LineFeature
{
  std::vector<FeaturePoint> points;
};

// A typed front line. Reuses FeaturePoint.
struct FrontLine
{
  FrontType type = FrontType::Stationary;
  std::vector<FeaturePoint> points;
};

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------

struct FrontOptions
{
  // Pre-allocated buffer capacities. Must be large enough to hold the
  // full detected result; upstream dynlib aborts with a Fortran stop
  // if these are exceeded.
  std::size_t max_points_per_type = 10000;
  std::size_t max_lines_per_type = 500;

  // Frontal intensity threshold in the input field's native units.
  // Negative preserves upstream default.
  double intensity_threshold = -1.0;

  // Speed threshold in m/s separating warm / stationary / cold.
  // Negative preserves upstream default.
  double speed_threshold = -1.0;

  // 1-2-1 smoothing passes on field + wind before detection. Negative
  // preserves upstream default.
  int smoothing_passes = -1;
};

struct LineOptions
{
  std::size_t max_points = 10000;
  std::size_t max_lines = 500;
  int smoothing_passes = -1;
};

// ---------------------------------------------------------------------------
// Front detection
// ---------------------------------------------------------------------------

// Jenkner et al. (2010) max-gradient front detection. Recommended
// default for modern NWP grids.
std::vector<FrontLine> detectFrontsMaxGrad(const Fmi::Matrix<double>& field,
                                           const Fmi::Matrix<double>& u,
                                           const Fmi::Matrix<double>& v,
                                           const Fmi::Matrix<double>& dx,
                                           const Fmi::Matrix<double>& dy,
                                           const FrontOptions& options = {});

// Berry et al. (2007) max-curvature front detection. Classic Hewson-
// style placement at the leading edge of a gradient zone.
std::vector<FrontLine> detectFrontsMaxCurv(const Fmi::Matrix<double>& field,
                                           const Fmi::Matrix<double>& u,
                                           const Fmi::Matrix<double>& v,
                                           const Fmi::Matrix<double>& dx,
                                           const Fmi::Matrix<double>& dy,
                                           const FrontOptions& options = {});

// ---------------------------------------------------------------------------
// Single-list line detectors (jet axes, conv/def/vor lines)
// ---------------------------------------------------------------------------

std::vector<LineFeature> detectJetAxes(const Fmi::Matrix<double>& u,
                                       const Fmi::Matrix<double>& v,
                                       const Fmi::Matrix<double>& dx,
                                       const Fmi::Matrix<double>& dy,
                                       const LineOptions& options = {});

// Jet axes with a fixed 30 m/s wind-speed threshold. Upstream
// `jetaxis_ff_thres`: same zero-shear condition as detectJetAxes
// but only accepts points where the smoothed wind speed exceeds
// 30 m/s. Useful on noisy grids where the curvature-only mask
// returns spurious lines in weak flow.
std::vector<LineFeature> detectJetAxesFFThres(const Fmi::Matrix<double>& u,
                                              const Fmi::Matrix<double>& v,
                                              const Fmi::Matrix<double>& dx,
                                              const Fmi::Matrix<double>& dy,
                                              const LineOptions& options = {});

std::vector<LineFeature> detectConvergenceLines(const Fmi::Matrix<double>& u,
                                                const Fmi::Matrix<double>& v,
                                                const Fmi::Matrix<double>& dx,
                                                const Fmi::Matrix<double>& dy,
                                                const LineOptions& options = {});

std::vector<LineFeature> detectDeformationLines(const Fmi::Matrix<double>& u,
                                                const Fmi::Matrix<double>& v,
                                                const Fmi::Matrix<double>& dx,
                                                const Fmi::Matrix<double>& dy,
                                                const LineOptions& options = {});

// Vorticity lines. Often used as a trough-axis proxy: pass (-u, -v) to
// reverse the sign convention if needed.
std::vector<LineFeature> detectVorticityLines(const Fmi::Matrix<double>& u,
                                              const Fmi::Matrix<double>& v,
                                              const Fmi::Matrix<double>& dx,
                                              const Fmi::Matrix<double>& dy,
                                              const LineOptions& options = {});

// Trough-axis alias — vorticity lines. dynlib does not ship a
// dedicated trough detector; Spensberger's group treats trough axes
// as maxima of cyclonic relative vorticity.
inline std::vector<LineFeature> detectTroughAxes(const Fmi::Matrix<double>& u,
                                                 const Fmi::Matrix<double>& v,
                                                 const Fmi::Matrix<double>& dx,
                                                 const Fmi::Matrix<double>& dy,
                                                 const LineOptions& options = {})
{
  return detectVorticityLines(u, v, dx, dy, options);
}

// ---------------------------------------------------------------------------
// Rossby wave breaking (local y-gradient reversal)
// ---------------------------------------------------------------------------

struct RwbOptions
{
  // Cutoff for the magnitude of d(field)/dy: a reversal is counted as
  // "significant" only when -d(field)/dy > ddy_thres. Must be >= 0. A
  // zero threshold flags every negative y-gradient as significant.
  double ddy_thres = 0.0;
};

// Per-cell outputs of the gradient-reversal RWB detector. All seven
// fields have the same shape as the input scalar field. Flag fields
// are 0/1 doubles (Fortran int8 internally, cast on return).
struct RwbResult
{
  Fmi::Matrix<double> anticyclonic_flag;
  Fmi::Matrix<double> cyclonic_flag;
  Fmi::Matrix<double> anticyclonic_gradmag;
  Fmi::Matrix<double> cyclonic_gradmag;
  Fmi::Matrix<double> anticyclonic_dfield_dy;
  Fmi::Matrix<double> cyclonic_dfield_dy;
  Fmi::Matrix<double> tested;
};

// Rossby wave breaking by local y-gradient reversal. `field` is
// typically isentropic PV; `latitudes` lists the latitude (degrees)
// of each grid row (must have size equal to field.height()). `mask`
// selects which cells are tested — pass nullptr to test everywhere.
RwbResult detectRossbyWaveBreakingGradRev(const Fmi::Matrix<double>& field,
                                          const std::vector<double>& latitudes,
                                          const Fmi::Matrix<double>& dx,
                                          const Fmi::Matrix<double>& dy,
                                          const RwbOptions& options = {},
                                          const Fmi::Matrix<double>* mask = nullptr);

// ---------------------------------------------------------------------------
// Cyclone detection (Wernli & Schwierz contour algorithm)
// ---------------------------------------------------------------------------

struct CycloneOptions
{
  // Maximum number of cyclones to find per call.
  std::size_t max_cyclones = 500;

  // Native-unit thresholds forwarded to dynlib's config module on
  // each call and restored afterwards. Negative means "keep the
  // upstream default".
  double min_size_km2 = -1.0;       // default 800
  double max_size_km2 = -1.0;       // default 4.5e6
  double max_orography_m = -1.0;    // default 1500
  double min_distance_km = -1.0;    // default 750
  double min_prominence = -1.0;     // default 200 (in units of MSL)
};

struct Cyclone
{
  double latitude;
  double longitude;
  double min_value;     // minimum MSL inside the cyclone
  double outer_value;   // MSL at the outermost contour
  double size_km2;
};

struct CycloneResult
{
  // Same shape as the input MSL. 0 = outside any cyclone; positive
  // integer-valued entries index into `cyclones` (1-based on input,
  // but indexed 0..size()-1 in the vector — cell value k marks the
  // (k-1)-th cyclone).
  Fmi::Matrix<double> mask;
  std::vector<Cyclone> cyclones;
};

// Cyclone detection from an MSL (or similar) field. `orography` has
// the same grid and is used to suppress minima over high terrain.
// `longitudes` (size nx) and `latitudes` (size ny) are degrees.
CycloneResult detectCyclonesByContour(const Fmi::Matrix<double>& msl,
                                      const Fmi::Matrix<double>& orography,
                                      const std::vector<double>& longitudes,
                                      const std::vector<double>& latitudes,
                                      const Fmi::Matrix<double>& dx,
                                      const Fmi::Matrix<double>& dy,
                                      const CycloneOptions& options = {});

// ---------------------------------------------------------------------------
// Precipitation blob detection
// ---------------------------------------------------------------------------

struct BlobOptions
{
  std::size_t max_blobs = 500;

  // Minimum distance between two blob centres (km). The algorithm
  // merges neighbours closer than this into the first-found blob.
  double min_distance_km = 100.0;
};

struct Blob
{
  double latitude;
  double longitude;
  double peak_value;     // maximum precip at the seed cell
  double outer_value;    // precip at the outermost accreted cell
  double size_km2;
};

struct BlobResult
{
  Fmi::Matrix<double> mask;
  std::vector<Blob> blobs;
};

BlobResult detectPrecipitationBlobs(const Fmi::Matrix<double>& precip,
                                    const std::vector<double>& longitudes,
                                    const std::vector<double>& latitudes,
                                    const Fmi::Matrix<double>& dx,
                                    const Fmi::Matrix<double>& dy,
                                    const BlobOptions& options = {});

// ---------------------------------------------------------------------------
// Blocking indicator (Masato et al. 2012)
// ---------------------------------------------------------------------------

// Returns a field of the same shape as the input. Positive entries
// mark candidate blocks. Typical input: 500 hPa geopotential or PV on
// 330 K. Downstream consumers threshold / contour as desired.
Fmi::Matrix<double> blockingIndicator(const Fmi::Matrix<double>& field,
                                      const Fmi::Matrix<double>& dx,
                                      const Fmi::Matrix<double>& dy);

// ---------------------------------------------------------------------------
// Grid helpers
// ---------------------------------------------------------------------------

// Build a double-grid-spacing matrix for a regular WGS84 lon/lat grid.
// lon0, lat0 are the grid origin, dlon and dlat the spacing in degrees.
// Outputs dx(i,j) = 2 * R * cos(lat(j)) * dlon_rad and dy(i,j) = 2 * R * dlat_rad.
void latLonDoubleGridSpacing(std::size_t nx,
                             std::size_t ny,
                             double lat0,
                             double dlat,
                             Fmi::Matrix<double>& dx,
                             Fmi::Matrix<double>& dy);

}  // namespace Dynlib
}  // namespace Fmi
