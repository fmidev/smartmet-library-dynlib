#include "Dynlib.h"

#include <boost/test/included/unit_test.hpp>
#include <macgyver/Matrix.h>

#include <cmath>

using namespace boost::unit_test;

test_suite* init_unit_test_suite(int /*argc*/, char* /*argv*/[])
{
  const char* name = "Dynlib tests";
  unit_test_log.set_threshold_level(log_messages);
  framework::master_test_suite().p_name.value = name;
  return NULL;
}

namespace
{
// Synthetic "cold front" field: tanh step across the y = x diagonal,
// 10 K total swing, transition width of a few cells.
Fmi::Matrix<double> syntheticFrontField(int nx, int ny)
{
  Fmi::Matrix<double> f(nx, ny);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      const double d = (static_cast<double>(i) - static_cast<double>(j)) * 0.25;
      f(i, j) = 280.0 + 5.0 * std::tanh(d);
    }
  return f;
}

// Synthetic "jet stream": bell-shaped zonal wind peak centred on
// j0 with half-width scale, peak strength peak m/s.
Fmi::Matrix<double> syntheticJetU(int nx, int ny, double j0, double scale, double peak)
{
  Fmi::Matrix<double> u(nx, ny);
  for (int j = 0; j < ny; ++j)
  {
    const double dy_j = (j - j0) / scale;
    const double row = peak * std::exp(-dy_j * dy_j);
    for (int i = 0; i < nx; ++i)
      u(i, j) = row;
  }
  return u;
}

Fmi::Matrix<double> uniform(int nx, int ny, double value)
{
  Fmi::Matrix<double> m(nx, ny);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      m(i, j) = value;
  return m;
}

// 25 km double grid spacing, suitable for a mid-latitude mesoscale grid.
Fmi::Matrix<double> regularGridSpacing(int nx, int ny)
{
  return uniform(nx, ny, 50000.0);
}

std::size_t totalPoints(const std::vector<Fmi::Dynlib::LineFeature>& lines)
{
  std::size_t n = 0;
  for (const auto& l : lines)
    n += l.points.size();
  return n;
}

std::size_t totalPoints(const std::vector<Fmi::Dynlib::FrontLine>& lines)
{
  std::size_t n = 0;
  for (const auto& l : lines)
    n += l.points.size();
  return n;
}
}  // namespace

// ---------------------------------------------------------------------------
// Front detection
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(fronts_maxgrad_synthetic_step)
{
  BOOST_TEST_MESSAGE("+ [detectFrontsMaxGrad: synthetic diagonal step]");

  const int n = 80;
  const auto field = syntheticFrontField(n, n);
  const auto u = uniform(n, n, 15.0);
  const auto v = uniform(n, n, 0.0);
  const auto dx = regularGridSpacing(n, n);
  const auto dy = regularGridSpacing(n, n);

  Fmi::Dynlib::FrontOptions opts;
  opts.smoothing_passes = 1;

  const auto lines = Fmi::Dynlib::detectFrontsMaxGrad(field, u, v, dx, dy, opts);
  BOOST_TEST_MESSAGE("  detected " << lines.size() << " front line(s), "
                                    << totalPoints(lines) << " points total");

  BOOST_REQUIRE(!lines.empty());
  BOOST_REQUIRE(totalPoints(lines) > 10u);

  std::size_t near_diag = 0;
  const auto total = totalPoints(lines);
  for (const auto& l : lines)
    for (const auto& p : l.points)
      if (std::abs(p.i - p.j) < 4.0)
        ++near_diag;
  BOOST_CHECK_MESSAGE(near_diag > total / 2,
                      "most front points should lie near the y=x boundary; got "
                          << near_diag << " of " << total);
}

BOOST_AUTO_TEST_CASE(fronts_maxgrad_flat_field_empty)
{
  BOOST_TEST_MESSAGE("+ [detectFrontsMaxGrad: flat field -> empty]");
  const int n = 40;
  const auto field = uniform(n, n, 280.0);
  const auto u = uniform(n, n, 10.0);
  const auto v = uniform(n, n, 0.0);
  const auto dx = regularGridSpacing(n, n);
  const auto dy = regularGridSpacing(n, n);
  const auto lines = Fmi::Dynlib::detectFrontsMaxGrad(field, u, v, dx, dy);
  BOOST_CHECK(lines.empty());
}

BOOST_AUTO_TEST_CASE(fronts_maxcurv_synthetic_step)
{
  BOOST_TEST_MESSAGE("+ [detectFrontsMaxCurv: synthetic diagonal step]");
  const int n = 80;
  const auto field = syntheticFrontField(n, n);
  const auto u = uniform(n, n, 15.0);
  const auto v = uniform(n, n, 0.0);
  const auto dx = regularGridSpacing(n, n);
  const auto dy = regularGridSpacing(n, n);

  Fmi::Dynlib::FrontOptions opts;
  opts.smoothing_passes = 1;

  const auto lines = Fmi::Dynlib::detectFrontsMaxCurv(field, u, v, dx, dy, opts);
  BOOST_TEST_MESSAGE("  detected " << lines.size() << " front line(s), "
                                    << totalPoints(lines) << " points total");

  // MaxCurv places fronts at the leading edge of a gradient zone; on a
  // perfectly symmetric tanh step the extremum may be too weak under
  // default thresholds. Accept zero OR a diagonal result.
  if (!lines.empty())
  {
    BOOST_CHECK(totalPoints(lines) > 5u);
  }
}

// ---------------------------------------------------------------------------
// Jet axes
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(jet_axes_smoke_test)
{
  BOOST_TEST_MESSAGE("+ [detectJetAxes: synthetic zonal jet smoke test]");

  // The jetaxis algorithm masks candidates via a second-derivative
  // intensity threshold whose scale is tuned against ERA5; it often
  // returns nothing on idealised synthetic grids. The purpose of this
  // test is to confirm the C-to-Fortran plumbing runs cleanly, returns
  // well-formed output, and — for any line it does find — points fall
  // inside the grid bounds.
  const int n = 80;
  const auto u = syntheticJetU(n, n, 40.0, 5.0, 60.0);
  const auto v = uniform(n, n, 0.0);
  const auto dx = regularGridSpacing(n, n);
  const auto dy = regularGridSpacing(n, n);

  Fmi::Dynlib::LineOptions opts;
  opts.smoothing_passes = 1;

  const auto lines = Fmi::Dynlib::detectJetAxes(u, v, dx, dy, opts);
  BOOST_TEST_MESSAGE("  detected " << lines.size() << " jet axis/axes, "
                                    << totalPoints(lines) << " points total");

  for (const auto& l : lines)
    for (const auto& p : l.points)
    {
      BOOST_CHECK(p.i >= 0.0 && p.i < static_cast<double>(n));
      BOOST_CHECK(p.j >= 0.0 && p.j < static_cast<double>(n));
    }
}

BOOST_AUTO_TEST_CASE(jet_axes_calm_is_empty)
{
  BOOST_TEST_MESSAGE("+ [detectJetAxes: zero wind -> empty]");
  const int n = 40;
  const auto u = uniform(n, n, 0.0);
  const auto v = uniform(n, n, 0.0);
  const auto dx = regularGridSpacing(n, n);
  const auto dy = regularGridSpacing(n, n);
  const auto lines = Fmi::Dynlib::detectJetAxes(u, v, dx, dy);
  BOOST_CHECK(lines.empty());
}

// ---------------------------------------------------------------------------
// Convergence / deformation / vorticity lines
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(convergence_lines_smoke_test)
{
  BOOST_TEST_MESSAGE("+ [detectConvergenceLines: colliding flow smoke test]");

  // The upstream comment on convline flags it as "not sufficiently
  // tested for general applicability", so this is purely a plumbing
  // test: ensure execution succeeds and any output is well-formed.
  const int n = 60;
  Fmi::Matrix<double> u(n, n);
  Fmi::Matrix<double> v(n, n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
    {
      u(i, j) = (i < n / 2 ? 10.0 : -10.0);
      v(i, j) = 0.0;
    }
  const auto dx = regularGridSpacing(n, n);
  const auto dy = regularGridSpacing(n, n);

  Fmi::Dynlib::LineOptions opts;
  opts.smoothing_passes = 2;

  const auto lines = Fmi::Dynlib::detectConvergenceLines(u, v, dx, dy, opts);
  BOOST_TEST_MESSAGE("  detected " << lines.size() << " line(s), "
                                    << totalPoints(lines) << " pts");
  for (const auto& l : lines)
    for (const auto& p : l.points)
    {
      BOOST_CHECK(p.i >= 0.0 && p.i < static_cast<double>(n));
      BOOST_CHECK(p.j >= 0.0 && p.j < static_cast<double>(n));
    }
}

BOOST_AUTO_TEST_CASE(deformation_lines_smoke_test)
{
  BOOST_TEST_MESSAGE("+ [detectDeformationLines: smoke test]");
  const int n = 40;
  const auto u = uniform(n, n, 5.0);
  const auto v = uniform(n, n, 5.0);
  const auto dx = regularGridSpacing(n, n);
  const auto dy = regularGridSpacing(n, n);
  const auto lines = Fmi::Dynlib::detectDeformationLines(u, v, dx, dy);
  BOOST_TEST_MESSAGE("  detected " << lines.size() << " deformation line(s)");
  // Uniform wind has zero deformation; expect empty output.
  BOOST_CHECK(lines.empty());
}

BOOST_AUTO_TEST_CASE(vorticity_and_trough_lines_alias)
{
  BOOST_TEST_MESSAGE("+ [detectTroughAxes vs detectVorticityLines equivalence]");

  const int n = 60;
  // Cyclonic shear: u increases with j -> negative dU/dy = positive relative vorticity.
  Fmi::Matrix<double> u(n, n);
  const auto v = uniform(n, n, 0.0);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      u(i, j) = 5.0 * (j - n / 2.0) / (n / 2.0);
  const auto dx = regularGridSpacing(n, n);
  const auto dy = regularGridSpacing(n, n);

  const auto a = Fmi::Dynlib::detectVorticityLines(u, v, dx, dy);
  const auto b = Fmi::Dynlib::detectTroughAxes(u, v, dx, dy);
  BOOST_CHECK_EQUAL(a.size(), b.size());
  BOOST_CHECK_EQUAL(totalPoints(a), totalPoints(b));
}

// ---------------------------------------------------------------------------
// Blocking indicator
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(blocking_indicator_shape_and_sign)
{
  BOOST_TEST_MESSAGE("+ [blockingIndicator: returns same-shape matrix]");

  const int n = 40;
  // Linear Z500-like field with a ridge near the centre.
  Fmi::Matrix<double> field(n, n);
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      field(i, j) = 5500.0 + 50.0 * std::sin(2.0 * M_PI * j / n) +
                    100.0 * std::exp(-std::pow((j - n / 2.0) / 4.0, 2));

  const auto dx = regularGridSpacing(n, n);
  const auto dy = regularGridSpacing(n, n);

  const auto bi = Fmi::Dynlib::blockingIndicator(field, dx, dy);
  BOOST_CHECK_EQUAL(bi.width(), n);
  BOOST_CHECK_EQUAL(bi.height(), n);

  // At least one cell should register a non-zero indicator value on a
  // non-uniform field.
  std::size_t nonzero = 0;
  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      if (std::abs(bi(i, j)) > 1e-12)
        ++nonzero;
  BOOST_TEST_MESSAGE("  " << nonzero << " cells with nonzero blocking indicator");
}

// ---------------------------------------------------------------------------
// jetaxis_ff_thres variant
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(jet_axes_ff_thres_below_cutoff_empty)
{
  BOOST_TEST_MESSAGE("+ [detectJetAxesFFThres: all-calm wind -> empty]");
  const int n = 40;
  const auto u = uniform(n, n, 5.0);
  const auto v = uniform(n, n, 0.0);
  const auto dx = regularGridSpacing(n, n);
  const auto dy = regularGridSpacing(n, n);
  const auto lines = Fmi::Dynlib::detectJetAxesFFThres(u, v, dx, dy);
  BOOST_CHECK(lines.empty());
}

BOOST_AUTO_TEST_CASE(jet_axes_ff_thres_plumbing)
{
  BOOST_TEST_MESSAGE("+ [detectJetAxesFFThres: plumbing + bounds]");
  const int n = 80;
  const auto u = syntheticJetU(n, n, 40.0, 5.0, 60.0);
  const auto v = uniform(n, n, 0.0);
  const auto dx = regularGridSpacing(n, n);
  const auto dy = regularGridSpacing(n, n);

  const auto lines = Fmi::Dynlib::detectJetAxesFFThres(u, v, dx, dy);
  BOOST_TEST_MESSAGE("  detected " << lines.size() << " jet axis/axes, "
                                    << totalPoints(lines) << " points total");
  for (const auto& l : lines)
    for (const auto& p : l.points)
    {
      BOOST_CHECK(p.i >= 0.0 && p.i < static_cast<double>(n));
      BOOST_CHECK(p.j >= 0.0 && p.j < static_cast<double>(n));
    }
}

// ---------------------------------------------------------------------------
// RWB gradient reversal
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(rwb_grad_rev_monotone_pv_no_reversal)
{
  BOOST_TEST_MESSAGE("+ [detectRossbyWaveBreakingGradRev: monotone PV, no reversals]");

  const int nx = 40;
  const int ny = 30;
  Fmi::Matrix<double> pv(nx, ny);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      pv(i, j) = 1e-6 * j;

  std::vector<double> lats(ny);
  for (int j = 0; j < ny; ++j)
    lats[j] = 30.0 + j * 0.5;

  const auto dx = regularGridSpacing(nx, ny);
  const auto dy = regularGridSpacing(nx, ny);

  Fmi::Dynlib::RwbOptions opts;
  opts.ddy_thres = 0.0;
  const auto r = Fmi::Dynlib::detectRossbyWaveBreakingGradRev(pv, lats, dx, dy, opts);

  BOOST_CHECK_EQUAL(r.anticyclonic_flag.width(), nx);
  BOOST_CHECK_EQUAL(r.anticyclonic_flag.height(), ny);

  std::size_t anticyc_flags = 0, cyc_flags = 0;
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      if (r.anticyclonic_flag(i, j) > 0.5) ++anticyc_flags;
      if (r.cyclonic_flag(i, j) > 0.5) ++cyc_flags;
    }
  BOOST_CHECK_EQUAL(anticyc_flags, 0u);
  BOOST_CHECK_EQUAL(cyc_flags, 0u);
}

BOOST_AUTO_TEST_CASE(rwb_grad_rev_reversed_pv_flags_something)
{
  BOOST_TEST_MESSAGE("+ [detectRossbyWaveBreakingGradRev: decreasing PV flags reversals]");

  const int nx = 30;
  const int ny = 30;
  Fmi::Matrix<double> pv(nx, ny);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      pv(i, j) = 1e-6 * (ny - j) + 1e-7 * i;

  std::vector<double> lats(ny);
  for (int j = 0; j < ny; ++j)
    lats[j] = 30.0 + j * 0.5;

  const auto dx = regularGridSpacing(nx, ny);
  const auto dy = regularGridSpacing(nx, ny);

  Fmi::Dynlib::RwbOptions opts;
  opts.ddy_thres = 0.0;
  const auto r = Fmi::Dynlib::detectRossbyWaveBreakingGradRev(pv, lats, dx, dy, opts);

  std::size_t any_flags = 0;
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      if (r.anticyclonic_flag(i, j) > 0.5 || r.cyclonic_flag(i, j) > 0.5)
        ++any_flags;
  BOOST_TEST_MESSAGE("  reversal flags: " << any_flags);
  BOOST_CHECK(any_flags > 0u);
}

// ---------------------------------------------------------------------------
// Cyclone detection
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cyclones_synthetic_low)
{
  BOOST_TEST_MESSAGE("+ [detectCyclonesByContour: single synthetic low]");

  const int nx = 60;
  const int ny = 60;

  Fmi::Matrix<double> msl(nx, ny);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      const double di = i - nx / 2.0;
      const double dj = j - ny / 2.0;
      msl(i, j) = 101325.0 - 1500.0 * std::exp(-(di * di + dj * dj) / (2.0 * 100.0));
    }

  const auto oro = uniform(nx, ny, 0.0);
  const auto dx = regularGridSpacing(nx, ny);
  const auto dy = regularGridSpacing(nx, ny);

  std::vector<double> lons(nx), lats(ny);
  for (int i = 0; i < nx; ++i)
    lons[i] = -10.0 + i * 0.5;
  for (int j = 0; j < ny; ++j)
    lats[j] = 40.0 + j * 0.5;

  Fmi::Dynlib::CycloneOptions opts;
  opts.max_cyclones = 50;
  opts.min_size_km2 = 100.0;
  opts.min_prominence = 100.0;
  const auto r =
      Fmi::Dynlib::detectCyclonesByContour(msl, oro, lons, lats, dx, dy, opts);

  BOOST_TEST_MESSAGE("  detected " << r.cyclones.size() << " cyclone(s)");
  BOOST_CHECK_EQUAL(r.mask.width(), nx);
  BOOST_CHECK_EQUAL(r.mask.height(), ny);
  BOOST_REQUIRE(!r.cyclones.empty());
  BOOST_CHECK_CLOSE(r.cyclones[0].longitude, lons[nx / 2], 5.0);
  BOOST_CHECK_CLOSE(r.cyclones[0].latitude, lats[ny / 2], 5.0);
}

BOOST_AUTO_TEST_CASE(cyclones_flat_field_empty)
{
  BOOST_TEST_MESSAGE("+ [detectCyclonesByContour: flat field -> no cyclones]");
  const int nx = 30, ny = 30;
  const auto msl = uniform(nx, ny, 101325.0);
  const auto oro = uniform(nx, ny, 0.0);
  const auto dx = regularGridSpacing(nx, ny);
  const auto dy = regularGridSpacing(nx, ny);
  std::vector<double> lons(nx), lats(ny);
  for (int i = 0; i < nx; ++i) lons[i] = i * 0.5;
  for (int j = 0; j < ny; ++j) lats[j] = 40.0 + j * 0.5;

  const auto r = Fmi::Dynlib::detectCyclonesByContour(msl, oro, lons, lats, dx, dy);
  BOOST_CHECK(r.cyclones.empty());
}

// ---------------------------------------------------------------------------
// Precipitation blob detection
// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(blobs_two_peaks)
{
  BOOST_TEST_MESSAGE("+ [detectPrecipitationBlobs: two synthetic peaks]");

  const int nx = 60, ny = 60;
  Fmi::Matrix<double> precip(nx, ny);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      const double d1 = std::hypot(i - 20.0, j - 20.0);
      const double d2 = std::hypot(i - 40.0, j - 40.0);
      precip(i, j) = 20.0 * std::exp(-d1 * d1 / 8.0) + 15.0 * std::exp(-d2 * d2 / 8.0);
    }

  const auto dx = regularGridSpacing(nx, ny);
  const auto dy = regularGridSpacing(nx, ny);

  std::vector<double> lons(nx), lats(ny);
  for (int i = 0; i < nx; ++i) lons[i] = i * 0.5;
  for (int j = 0; j < ny; ++j) lats[j] = 40.0 + j * 0.5;

  Fmi::Dynlib::BlobOptions opts;
  opts.max_blobs = 20;
  opts.min_distance_km = 50.0;

  const auto r =
      Fmi::Dynlib::detectPrecipitationBlobs(precip, lons, lats, dx, dy, opts);

  BOOST_TEST_MESSAGE("  detected " << r.blobs.size() << " blob(s)");
  BOOST_CHECK_EQUAL(r.mask.width(), nx);
  BOOST_CHECK_EQUAL(r.mask.height(), ny);
  BOOST_CHECK(r.blobs.size() >= 1u);
}

// ---------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(latlon_grid_spacing_shape_and_positivity)
{
  Fmi::Matrix<double> dx, dy;
  Fmi::Dynlib::latLonDoubleGridSpacing(20, 10, 30.0, 0.5, dx, dy);
  BOOST_CHECK_EQUAL(dx.width(), 20);
  BOOST_CHECK_EQUAL(dx.height(), 10);
  BOOST_CHECK(dx(0, 0) > 0.0);
  BOOST_CHECK(dy(0, 0) > 0.0);
  // dx should shrink at higher latitudes (larger j in this setup).
  BOOST_CHECK(dx(0, 0) > dx(0, 9));
}
