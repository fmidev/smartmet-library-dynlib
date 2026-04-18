#include "Dynlib.h"
#include "DynlibC.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace Fmi
{
namespace Dynlib
{

namespace
{
constexpr double kEarthRadius = 6371000.0;
constexpr double kPi = 3.14159265358979323846;

void checkSameShape(const Fmi::Matrix<double>& a,
                    const Fmi::Matrix<double>& b,
                    const char* name)
{
  if (a.width() != b.width() || a.height() != b.height())
    throw std::invalid_argument(std::string("Dynlib: shape mismatch for ") + name);
}

// Decode a single list of (j, i, value) points into a vector of
// LineFeature. The offsets in off_t are scanned as a monotonically
// non-decreasing sequence in [0, no]; the loop stops at the first
// NaN, out-of-range value, or non-monotone step.
std::vector<LineFeature> decodeLines(const double* off_t,
                                     const double* pts_t,
                                     int32_t nf,
                                     int32_t no)
{
  // Upstream line_locate declares the offset buffer as intent(out),
  // which allows the compiler to clobber any pre-call NaN fill — in
  // practice uninitialised slots contain stack garbage (often
  // subnormal doubles). We therefore accept only values that are
  // actual non-negative integers in range. The sentinel at the end
  // of the valid prefix equals the total point count, so a value
  // that exceeds `no` is also treated as out-of-range garbage.
  std::vector<int32_t> offsets;
  offsets.reserve(16);
  double prev = -1.0;
  for (int32_t f = 0; f < nf; ++f)
  {
    const double v = off_t[f];
    if (!std::isfinite(v))
      break;
    if (v < 0.0 || v > static_cast<double>(no))
      break;
    if (v != std::floor(v))   // garbage subnormals fail this check
      break;
    if (v < prev)             // offsets must be non-decreasing
      break;
    offsets.push_back(static_cast<int32_t>(v));
    prev = v;
  }

  std::vector<LineFeature> out;
  for (std::size_t k = 0; k + 1 < offsets.size(); ++k)
  {
    const int32_t start = offsets[k];
    const int32_t end = offsets[k + 1];
    if (end <= start)
      continue;

    LineFeature line;
    line.points.reserve(static_cast<std::size_t>(end - start));
    for (int32_t p = start; p < end; ++p)
    {
      const double j = pts_t[p * 3 + 0];
      const double i = pts_t[p * 3 + 1];
      const double val = pts_t[p * 3 + 2];
      if (std::isnan(j) || std::isnan(i))
        break;
      line.points.push_back({j, i, val});
    }
    if (!line.points.empty())
      out.push_back(std::move(line));
  }
  return out;
}

// Three-type variant used for fronts. t = 0 cold, 1 warm, 2 stationary.
std::vector<FrontLine> decodeTypedFronts(const double* off,
                                         const double* pts,
                                         int32_t nf,
                                         int32_t no)
{
  std::vector<FrontLine> out;
  for (int32_t t = 0; t < 3; ++t)
  {
    const double* off_t = &off[static_cast<std::size_t>(t) * nf];
    const double* pts_t = &pts[static_cast<std::size_t>(t) * no * 3];
    auto lines = decodeLines(off_t, pts_t, nf, no);
    for (auto& l : lines)
    {
      FrontLine fl;
      fl.type = static_cast<FrontType>(t);
      fl.points = std::move(l.points);
      out.push_back(std::move(fl));
    }
  }
  return out;
}

// Validate input shapes against the reference field and return (nx, ny)
// as int32_t. Throws if grid is too small.
std::pair<int32_t, int32_t> checkAndGetDims(const Fmi::Matrix<double>& ref,
                                            const Fmi::Matrix<double>* const* peers,
                                            const char* const* names,
                                            std::size_t n_peers)
{
  for (std::size_t k = 0; k < n_peers; ++k)
    checkSameShape(ref, *peers[k], names[k]);
  const int32_t nx = static_cast<int32_t>(ref.width());
  const int32_t ny = static_cast<int32_t>(ref.height());
  if (nx < 3 || ny < 3)
    throw std::invalid_argument("Dynlib: grid must be at least 3x3");
  return {nx, ny};
}

}  // namespace

// ---------------------------------------------------------------------------
// Front detection
// ---------------------------------------------------------------------------

std::vector<FrontLine> detectFrontsMaxGrad(const Fmi::Matrix<double>& field,
                                           const Fmi::Matrix<double>& u,
                                           const Fmi::Matrix<double>& v,
                                           const Fmi::Matrix<double>& dx,
                                           const Fmi::Matrix<double>& dy,
                                           const FrontOptions& options)
{
  const Fmi::Matrix<double>* peers[] = {&u, &v, &dx, &dy};
  const char* names[] = {"u", "v", "dx", "dy"};
  auto [nx, ny] = checkAndGetDims(field, peers, names, 4);

  const int32_t no = static_cast<int32_t>(options.max_points_per_type);
  const int32_t nf = static_cast<int32_t>(options.max_lines_per_type);

  std::vector<double> pts(static_cast<std::size_t>(3) * no * 3);
  std::vector<double> off(static_cast<std::size_t>(3) * nf);

  dynlib_detect_fronts_maxgrad(nx, ny, no, nf,
                               &field(0, 0), &u(0, 0), &v(0, 0),
                               &dx(0, 0), &dy(0, 0),
                               pts.data(), off.data(),
                               options.intensity_threshold,
                               options.speed_threshold,
                               options.smoothing_passes);

  return decodeTypedFronts(off.data(), pts.data(), nf, no);
}

std::vector<FrontLine> detectFrontsMaxCurv(const Fmi::Matrix<double>& field,
                                           const Fmi::Matrix<double>& u,
                                           const Fmi::Matrix<double>& v,
                                           const Fmi::Matrix<double>& dx,
                                           const Fmi::Matrix<double>& dy,
                                           const FrontOptions& options)
{
  const Fmi::Matrix<double>* peers[] = {&u, &v, &dx, &dy};
  const char* names[] = {"u", "v", "dx", "dy"};
  auto [nx, ny] = checkAndGetDims(field, peers, names, 4);

  const int32_t no = static_cast<int32_t>(options.max_points_per_type);
  const int32_t nf = static_cast<int32_t>(options.max_lines_per_type);

  std::vector<double> pts(static_cast<std::size_t>(3) * no * 3);
  std::vector<double> off(static_cast<std::size_t>(3) * nf);

  dynlib_detect_fronts_maxcurv(nx, ny, no, nf,
                               &field(0, 0), &u(0, 0), &v(0, 0),
                               &dx(0, 0), &dy(0, 0),
                               pts.data(), off.data(),
                               options.intensity_threshold,
                               options.speed_threshold,
                               options.smoothing_passes);

  return decodeTypedFronts(off.data(), pts.data(), nf, no);
}

// ---------------------------------------------------------------------------
// Jet axes and line detectors
// ---------------------------------------------------------------------------

namespace
{
// kind_code 0 = jet axes, 1 = convergence, 2 = deformation, 3 = vorticity,
// 4 = jet axes with ff>=30 m/s mask.
std::vector<LineFeature> runLineDetector(int32_t kind_code,
                                         const Fmi::Matrix<double>& u,
                                         const Fmi::Matrix<double>& v,
                                         const Fmi::Matrix<double>& dx,
                                         const Fmi::Matrix<double>& dy,
                                         const LineOptions& options)
{
  const Fmi::Matrix<double>* peers[] = {&v, &dx, &dy};
  const char* names[] = {"v", "dx", "dy"};
  auto [nx, ny] = checkAndGetDims(u, peers, names, 3);

  const int32_t no = static_cast<int32_t>(options.max_points);
  const int32_t nf = static_cast<int32_t>(options.max_lines);

  std::vector<double> pts(static_cast<std::size_t>(no) * 3);
  std::vector<double> off(static_cast<std::size_t>(nf));

  if (kind_code == 0 || kind_code == 4)
  {
    const int32_t variant = (kind_code == 4 ? 1 : 0);
    dynlib_detect_jet_axes(variant, nx, ny, no, nf,
                           &u(0, 0), &v(0, 0), &dx(0, 0), &dy(0, 0),
                           pts.data(), off.data(), options.smoothing_passes);
  }
  else
    dynlib_detect_lines(kind_code, nx, ny, no, nf,
                        &u(0, 0), &v(0, 0), &dx(0, 0), &dy(0, 0),
                        pts.data(), off.data(), options.smoothing_passes);

  return decodeLines(off.data(), pts.data(), nf, no);
}
}  // namespace

std::vector<LineFeature> detectJetAxes(const Fmi::Matrix<double>& u,
                                       const Fmi::Matrix<double>& v,
                                       const Fmi::Matrix<double>& dx,
                                       const Fmi::Matrix<double>& dy,
                                       const LineOptions& options)
{
  return runLineDetector(0, u, v, dx, dy, options);
}

std::vector<LineFeature> detectJetAxesFFThres(const Fmi::Matrix<double>& u,
                                              const Fmi::Matrix<double>& v,
                                              const Fmi::Matrix<double>& dx,
                                              const Fmi::Matrix<double>& dy,
                                              const LineOptions& options)
{
  return runLineDetector(4, u, v, dx, dy, options);
}

std::vector<LineFeature> detectConvergenceLines(const Fmi::Matrix<double>& u,
                                                const Fmi::Matrix<double>& v,
                                                const Fmi::Matrix<double>& dx,
                                                const Fmi::Matrix<double>& dy,
                                                const LineOptions& options)
{
  return runLineDetector(1, u, v, dx, dy, options);
}

std::vector<LineFeature> detectDeformationLines(const Fmi::Matrix<double>& u,
                                                const Fmi::Matrix<double>& v,
                                                const Fmi::Matrix<double>& dx,
                                                const Fmi::Matrix<double>& dy,
                                                const LineOptions& options)
{
  return runLineDetector(2, u, v, dx, dy, options);
}

std::vector<LineFeature> detectVorticityLines(const Fmi::Matrix<double>& u,
                                              const Fmi::Matrix<double>& v,
                                              const Fmi::Matrix<double>& dx,
                                              const Fmi::Matrix<double>& dy,
                                              const LineOptions& options)
{
  return runLineDetector(3, u, v, dx, dy, options);
}

// ---------------------------------------------------------------------------
// Rossby wave breaking (gradient-reversal)
// ---------------------------------------------------------------------------

RwbResult detectRossbyWaveBreakingGradRev(const Fmi::Matrix<double>& field,
                                          const std::vector<double>& latitudes,
                                          const Fmi::Matrix<double>& dx,
                                          const Fmi::Matrix<double>& dy,
                                          const RwbOptions& options,
                                          const Fmi::Matrix<double>* mask)
{
  const Fmi::Matrix<double>* peers[] = {&dx, &dy};
  const char* names[] = {"dx", "dy"};
  auto [nx, ny] = checkAndGetDims(field, peers, names, 2);

  if (latitudes.size() != static_cast<std::size_t>(ny))
    throw std::invalid_argument("Dynlib: latitudes length must equal field.height()");

  if (options.ddy_thres < 0.0)
    throw std::invalid_argument("Dynlib: RwbOptions.ddy_thres must be non-negative");

  Fmi::Matrix<double> mask_all(nx, ny);
  if (mask != nullptr)
  {
    checkSameShape(field, *mask, "mask");
  }
  else
  {
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        mask_all(i, j) = 1.0;
  }

  RwbResult out;
  out.anticyclonic_flag = Fmi::Matrix<double>(nx, ny);
  out.cyclonic_flag = Fmi::Matrix<double>(nx, ny);
  out.anticyclonic_gradmag = Fmi::Matrix<double>(nx, ny);
  out.cyclonic_gradmag = Fmi::Matrix<double>(nx, ny);
  out.anticyclonic_dfield_dy = Fmi::Matrix<double>(nx, ny);
  out.cyclonic_dfield_dy = Fmi::Matrix<double>(nx, ny);
  out.tested = Fmi::Matrix<double>(nx, ny);

  const double* mask_ptr = (mask != nullptr ? &(*mask)(0, 0) : &mask_all(0, 0));

  dynlib_detect_rwb_grad_rev(nx, ny,
                             &field(0, 0), mask_ptr, latitudes.data(),
                             options.ddy_thres,
                             &dx(0, 0), &dy(0, 0),
                             &out.anticyclonic_flag(0, 0),
                             &out.cyclonic_flag(0, 0),
                             &out.anticyclonic_gradmag(0, 0),
                             &out.cyclonic_gradmag(0, 0),
                             &out.anticyclonic_dfield_dy(0, 0),
                             &out.cyclonic_dfield_dy(0, 0),
                             &out.tested(0, 0));
  return out;
}

// ---------------------------------------------------------------------------
// Cyclone / precipitation-blob helpers (shared sort step)
// ---------------------------------------------------------------------------

namespace
{
// Build the (sorted values, iis, jjs) triple from a field. `ascending`
// true sorts low-to-high (cyclones want this); false sorts high-to-low.
// The triple is what dynlib's cyclone / blob Fortran subroutines
// consume — they walk sorted_values in order, re-reading the full
// field for neighbour lookups.
struct SortedField
{
  std::vector<double> sorted_values;
  std::vector<int32_t> iis;
  std::vector<int32_t> jjs;
};

SortedField buildSortedField(const Fmi::Matrix<double>& field, bool ascending)
{
  const int32_t nx = static_cast<int32_t>(field.width());
  const int32_t ny = static_cast<int32_t>(field.height());
  const std::size_t n = static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny);

  std::vector<std::size_t> order(n);
  for (std::size_t k = 0; k < n; ++k)
    order[k] = k;

  const double* data = &field(0, 0);
  if (ascending)
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return data[a] < data[b]; });
  else
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return data[a] > data[b]; });

  SortedField out;
  out.sorted_values.resize(n);
  out.iis.resize(n);
  out.jjs.resize(n);
  for (std::size_t k = 0; k < n; ++k)
  {
    const std::size_t lin = order[k];
    out.sorted_values[k] = data[lin];
    out.iis[k] = static_cast<int32_t>(lin % static_cast<std::size_t>(nx));
    out.jjs[k] = static_cast<int32_t>(lin / static_cast<std::size_t>(nx));
  }
  return out;
}
}  // namespace

// ---------------------------------------------------------------------------
// Cyclone detection
// ---------------------------------------------------------------------------

CycloneResult detectCyclonesByContour(const Fmi::Matrix<double>& msl,
                                      const Fmi::Matrix<double>& orography,
                                      const std::vector<double>& longitudes,
                                      const std::vector<double>& latitudes,
                                      const Fmi::Matrix<double>& dx,
                                      const Fmi::Matrix<double>& dy,
                                      const CycloneOptions& options)
{
  const Fmi::Matrix<double>* peers[] = {&orography, &dx, &dy};
  const char* names[] = {"orography", "dx", "dy"};
  auto [nx, ny] = checkAndGetDims(msl, peers, names, 3);

  if (longitudes.size() != static_cast<std::size_t>(nx))
    throw std::invalid_argument("Dynlib: longitudes length must equal msl.width()");
  if (latitudes.size() != static_cast<std::size_t>(ny))
    throw std::invalid_argument("Dynlib: latitudes length must equal msl.height()");

  const int32_t nn = static_cast<int32_t>(options.max_cyclones);

  const auto sorted = buildSortedField(msl, /*ascending=*/true);

  CycloneResult out;
  out.mask = Fmi::Matrix<double>(nx, ny);
  std::vector<double> meta(static_cast<std::size_t>(5) * static_cast<std::size_t>(nn), 0.0);

  dynlib_detect_cyclones(nx, ny, nn,
                         &msl(0, 0),
                         sorted.sorted_values.data(),
                         sorted.iis.data(),
                         sorted.jjs.data(),
                         &orography(0, 0),
                         longitudes.data(),
                         latitudes.data(),
                         &dx(0, 0), &dy(0, 0),
                         options.min_size_km2,
                         options.max_size_km2,
                         options.max_orography_m,
                         options.min_distance_km,
                         options.min_prominence,
                         &out.mask(0, 0),
                         meta.data());

  for (int32_t m = 0; m < nn; ++m)
  {
    const double lat = meta[m * 5 + 0];
    const double lon = meta[m * 5 + 1];
    const double minv = meta[m * 5 + 2];
    const double outer = meta[m * 5 + 3];
    const double size = meta[m * 5 + 4];
    // Upstream zero-initialises unused slots.
    if (lat == 0.0 && lon == 0.0 && minv == 0.0 && outer == 0.0 && size == 0.0)
      break;
    out.cyclones.push_back({lat, lon, minv, outer, size});
  }
  return out;
}

// ---------------------------------------------------------------------------
// Precipitation blob detection
// ---------------------------------------------------------------------------

BlobResult detectPrecipitationBlobs(const Fmi::Matrix<double>& precip,
                                    const std::vector<double>& longitudes,
                                    const std::vector<double>& latitudes,
                                    const Fmi::Matrix<double>& dx,
                                    const Fmi::Matrix<double>& dy,
                                    const BlobOptions& options)
{
  const Fmi::Matrix<double>* peers[] = {&dx, &dy};
  const char* names[] = {"dx", "dy"};
  auto [nx, ny] = checkAndGetDims(precip, peers, names, 2);

  if (longitudes.size() != static_cast<std::size_t>(nx))
    throw std::invalid_argument("Dynlib: longitudes length must equal precip.width()");
  if (latitudes.size() != static_cast<std::size_t>(ny))
    throw std::invalid_argument("Dynlib: latitudes length must equal precip.height()");

  const int32_t nn = static_cast<int32_t>(options.max_blobs);

  // Upstream triggers on `val < 0` in the sorted array. Build a
  // sign-flipped ascending sort so that the strongest precipitation
  // cells end up first as the most negative entries.
  Fmi::Matrix<double> neg(nx, ny);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      neg(i, j) = -precip(i, j);
  const auto sorted = buildSortedField(neg, /*ascending=*/true);

  BlobResult out;
  out.mask = Fmi::Matrix<double>(nx, ny);
  std::vector<double> meta(static_cast<std::size_t>(5) * static_cast<std::size_t>(nn), 0.0);

  dynlib_detect_blobs(nx, ny, nn,
                      &precip(0, 0),
                      sorted.sorted_values.data(),
                      sorted.iis.data(),
                      sorted.jjs.data(),
                      longitudes.data(),
                      latitudes.data(),
                      &dx(0, 0), &dy(0, 0),
                      options.min_distance_km,
                      &out.mask(0, 0),
                      meta.data());

  for (int32_t m = 0; m < nn; ++m)
  {
    const double lat = meta[m * 5 + 0];
    const double lon = meta[m * 5 + 1];
    const double peak = -meta[m * 5 + 2];   // flip sign back
    const double outer = -meta[m * 5 + 3];
    const double size = meta[m * 5 + 4];
    if (lat == 0.0 && lon == 0.0 && size == 0.0)
      break;
    out.blobs.push_back({lat, lon, peak, outer, size});
  }
  return out;
}

// ---------------------------------------------------------------------------
// Blocking indicator
// ---------------------------------------------------------------------------

Fmi::Matrix<double> blockingIndicator(const Fmi::Matrix<double>& field,
                                      const Fmi::Matrix<double>& dx,
                                      const Fmi::Matrix<double>& dy)
{
  const Fmi::Matrix<double>* peers[] = {&dx, &dy};
  const char* names[] = {"dx", "dy"};
  auto [nx, ny] = checkAndGetDims(field, peers, names, 2);

  Fmi::Matrix<double> out(nx, ny);
  dynlib_block_indicator(nx, ny,
                         &field(0, 0), &dx(0, 0), &dy(0, 0),
                         &out(0, 0));
  return out;
}

// ---------------------------------------------------------------------------
// Grid helpers
// ---------------------------------------------------------------------------

void latLonDoubleGridSpacing(std::size_t nx,
                             std::size_t ny,
                             double lat0,
                             double dlat,
                             Fmi::Matrix<double>& dx,
                             Fmi::Matrix<double>& dy)
{
  // dx(i,j) = 2 * R * cos(lat_j) * dlon_rad
  // dy(i,j) = 2 * R * dlat_rad (independent of i,j on a regular grid)
  const double dlon_rad = (2.0 * kPi / 360.0) * dlat;  // caller usually sets dlon = dlat
  const double dlat_rad = (2.0 * kPi / 360.0) * dlat;
  const double dy_const = 2.0 * kEarthRadius * dlat_rad;

  dx = Fmi::Matrix<double>(static_cast<int>(nx), static_cast<int>(ny));
  dy = Fmi::Matrix<double>(static_cast<int>(nx), static_cast<int>(ny));

  for (int j = 0; j < static_cast<int>(ny); ++j)
  {
    const double lat_j = lat0 + j * dlat;
    const double lat_r = (2.0 * kPi / 360.0) * lat_j;
    const double dx_j = 2.0 * kEarthRadius * std::cos(lat_r) * dlon_rad;
    for (int i = 0; i < static_cast<int>(nx); ++i)
    {
      dx(i, j) = dx_j;
      dy(i, j) = dy_const;
    }
  }
}

}  // namespace Dynlib
}  // namespace Fmi
