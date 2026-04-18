// End-to-end smoke test against a real ECMWF querydata file.
//
// Reads ~/hub/europe_pressurelevels_12utc.sqd (three 12 UTC Europe
// snapshots on pressure levels), runs every detector that is covered
// by the available fields, and writes a GeoJSON FeatureCollection
// per timestep containing the detected features in lon/lat. Blocking
// is reported as summary statistics since it is a scalar field and
// belongs on an isoband renderer.
//
// Caveats on the rotated lat/lon sample file:
//   * Line LOCATIONS are geographically correct. Grid-relative
//     finite differences find the same zero crossings regardless
//     of projection.
//   * Cold/warm classification depends on the sign of u · ∇θ. The
//     detector computes ∇θ in grid-relative coordinates, but the
//     querydata wind components u, v are in geographic orientation
//     (ECMWF GRIB convention). For rotated grids they are NOT
//     aligned, so the sign of frspd may flip near the map edges.
//     A proper fix would rotate u, v into grid-relative coordinates
//     before the dynlib call — deferred until we integrate this with
//     an engine that exposes the projection matrix.
//
// Build:  make examples
// Run:    ./examples/demo [path-to-sqd]   (default: ~/hub/europe_pressurelevels_12utc.sqd)

#include "Dynlib.h"

#include <newbase/NFmiFastQueryInfo.h>
#include <newbase/NFmiGrid.h>
#include <newbase/NFmiLevel.h>
#include <newbase/NFmiLevelType.h>
#include <newbase/NFmiParameterName.h>
#include <newbase/NFmiQueryData.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
constexpr double kEarthRadius = 6371000.0;
constexpr double kPi = 3.14159265358979323846;

double toRad(double deg) { return deg * kPi / 180.0; }

// Great-circle distance between two lon/lat points in metres.
double haversine(double lon1, double lat1, double lon2, double lat2)
{
  const double dlat = toRad(lat2 - lat1);
  const double dlon = toRad(lon2 - lon1);
  const double a = std::sin(dlat / 2) * std::sin(dlat / 2) +
                   std::cos(toRad(lat1)) * std::cos(toRad(lat2)) *
                       std::sin(dlon / 2) * std::sin(dlon / 2);
  return 2.0 * kEarthRadius * std::asin(std::min(1.0, std::sqrt(a)));
}

struct GridView
{
  int nx = 0;
  int ny = 0;
  std::vector<double> lons;  // size nx*ny, row-major (i fastest)
  std::vector<double> lats;
  Fmi::Matrix<double> dx;    // double grid spacing, metres
  Fmi::Matrix<double> dy;

  double lon(int i, int j) const { return lons[i + nx * j]; }
  double lat(int i, int j) const { return lats[i + nx * j]; }
};

// Build the lon/lat arrays and double-grid-spacing matrices for the
// currently-selected time/param of `info`.
GridView buildGridView(NFmiFastQueryInfo& info)
{
  GridView gv;
  gv.nx = static_cast<int>(info.GridXNumber());
  gv.ny = static_cast<int>(info.GridYNumber());
  gv.lons.resize(static_cast<std::size_t>(gv.nx) * gv.ny);
  gv.lats.resize(gv.lons.size());

  unsigned long idx = 0;
  for (int j = 0; j < gv.ny; ++j)
    for (int i = 0; i < gv.nx; ++i)
    {
      const NFmiPoint ll = info.LatLon(idx++);
      gv.lons[i + gv.nx * j] = ll.X();
      gv.lats[i + gv.nx * j] = ll.Y();
    }

  gv.dx = Fmi::Matrix<double>(gv.nx, gv.ny);
  gv.dy = Fmi::Matrix<double>(gv.nx, gv.ny);
  for (int j = 0; j < gv.ny; ++j)
    for (int i = 0; i < gv.nx; ++i)
    {
      const int im = (i == 0 ? 0 : i - 1);
      const int ip = (i == gv.nx - 1 ? gv.nx - 1 : i + 1);
      const int jm = (j == 0 ? 0 : j - 1);
      const int jp = (j == gv.ny - 1 ? gv.ny - 1 : j + 1);
      const double dxm = haversine(gv.lon(im, j), gv.lat(im, j),
                                    gv.lon(ip, j), gv.lat(ip, j));
      const double dym = haversine(gv.lon(i, jm), gv.lat(i, jm),
                                    gv.lon(i, jp), gv.lat(i, jp));
      // Mirror-boundary: scale up to "double" spacing when at the edge
      // (neighbour on one side was clamped to self).
      gv.dx(i, j) = (i == 0 || i == gv.nx - 1) ? dxm * 2.0 : dxm;
      gv.dy(i, j) = (j == 0 || j == gv.ny - 1) ? dym * 2.0 : dym;
    }
  return gv;
}

// Bolton (1980) equivalent potential temperature. T_K in Kelvin,
// p_hPa in hectopascals, RH_pct in percent (0..100). Returns θe in
// Kelvin. Uses Bolton equations 10 (e_s), 22 (T_L), 43 (θe). θe has
// sharper gradients at frontal zones than dry θ because it collapses
// the latent-heat contrast between moist warm and dry cold air masses
// into a single scalar.
double boltonThetaE(double T_K, double p_hPa, double RH_pct)
{
  const double RH = std::max(1e-3, RH_pct) / 100.0;
  const double T_L = 1.0 / (1.0 / (T_K - 55.0) - std::log(RH) / 2840.0) + 55.0;
  const double T_C = T_K - 273.15;
  const double e_s = 6.112 * std::exp(17.67 * T_C / (T_C + 243.5));   // hPa
  const double e = RH * e_s;                                           // hPa
  const double r = 622.0 * e / std::max(1e-6, p_hPa - e);              // g/kg
  return T_K * std::pow(1000.0 / p_hPa, 0.2854 * (1.0 - 0.28e-3 * r)) *
         std::exp(r * (1.0 + 0.81e-3 * r) * (3.376 / T_L - 0.00254));
}

// Extract the current time/param/level's grid values into an
// Fmi::Matrix<double>. Missing cells are filled with the nearest
// valid neighbour, since dynlib's numerics do not handle NaN inputs
// gracefully. Returns the fraction of missing cells (caller can warn).
double extractField(NFmiFastQueryInfo& info, Fmi::Matrix<double>& out)
{
  const int nx = static_cast<int>(info.GridXNumber());
  const int ny = static_cast<int>(info.GridYNumber());
  out = Fmi::Matrix<double>(nx, ny);
  std::size_t miss = 0;
  double fallback = 0.0;
  unsigned long idx = 0;
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      info.LocationIndex(idx++);
      const float v = info.FloatValue();
      if (v == kFloatMissing)
      {
        ++miss;
        out(i, j) = fallback;
      }
      else
      {
        out(i, j) = static_cast<double>(v);
        fallback = out(i, j);
      }
    }
  return static_cast<double>(miss) / (nx * ny);
}

// Select (param, level). Returns true on success. Prints a warning to
// stderr and returns false if not found.
bool select(NFmiFastQueryInfo& info, FmiParameterName param, double pressure_hpa)
{
  if (!info.Param(param))
  {
    std::cerr << "  warning: parameter " << param << " not in file\n";
    return false;
  }
  if (!info.Level(NFmiLevel(kFmiPressureLevel, static_cast<float>(pressure_hpa))))
  {
    std::cerr << "  warning: pressure level " << pressure_hpa << " not in file\n";
    return false;
  }
  return true;
}

// Convert a feature point (fractional j, i) into geographic lon/lat by
// bilinear interpolation of the grid's own lon/lat matrices.
void toLonLat(const GridView& gv, double j, double i, double& lon, double& lat)
{
  const int i0 = std::max(0, std::min(gv.nx - 1, static_cast<int>(std::floor(i))));
  const int j0 = std::max(0, std::min(gv.ny - 1, static_cast<int>(std::floor(j))));
  const int i1 = std::min(gv.nx - 1, i0 + 1);
  const int j1 = std::min(gv.ny - 1, j0 + 1);
  const double ai = i - i0;
  const double aj = j - j0;
  const double lon00 = gv.lon(i0, j0), lon10 = gv.lon(i1, j0);
  const double lon01 = gv.lon(i0, j1), lon11 = gv.lon(i1, j1);
  const double lat00 = gv.lat(i0, j0), lat10 = gv.lat(i1, j0);
  const double lat01 = gv.lat(i0, j1), lat11 = gv.lat(i1, j1);
  lon = (1 - ai) * (1 - aj) * lon00 + ai * (1 - aj) * lon10 +
        (1 - ai) * aj * lon01 + ai * aj * lon11;
  lat = (1 - ai) * (1 - aj) * lat00 + ai * (1 - aj) * lat10 +
        (1 - ai) * aj * lat01 + ai * aj * lat11;
}

// Append a GeoJSON LineString feature from a detector point list.
void writeLineString(std::ostream& os,
                     bool& first,
                     const GridView& gv,
                     const std::vector<Fmi::Dynlib::FeaturePoint>& pts,
                     const std::string& detector,
                     const std::string& subtype = "")
{
  if (pts.size() < 2) return;
  if (!first) os << ",\n";
  first = false;
  os << "{\"type\":\"Feature\",\"properties\":{\"detector\":\"" << detector << "\"";
  if (!subtype.empty())
    os << ",\"subtype\":\"" << subtype << "\"";
  os << "},\"geometry\":{\"type\":\"LineString\",\"coordinates\":[";
  for (std::size_t k = 0; k < pts.size(); ++k)
  {
    double lon, lat;
    toLonLat(gv, pts[k].j, pts[k].i, lon, lat);
    if (k > 0) os << ",";
    os << "[" << lon << "," << lat << "]";
  }
  os << "]}}";
}

const char* frontTypeName(Fmi::Dynlib::FrontType t)
{
  switch (t)
  {
    case Fmi::Dynlib::FrontType::Cold: return "cold";
    case Fmi::Dynlib::FrontType::Warm: return "warm";
    default: return "stationary";
  }
}

// Per-timestep bundle of detected features for later SVG rendering.
struct FeatureBundle
{
  std::string stamp;
  std::vector<Fmi::Dynlib::FrontLine> fronts;
  std::vector<Fmi::Dynlib::LineFeature> jets_default;
  std::vector<Fmi::Dynlib::LineFeature> jets_ff30;
  std::vector<Fmi::Dynlib::LineFeature> convergence;
  std::vector<Fmi::Dynlib::LineFeature> deformation;
  std::vector<Fmi::Dynlib::LineFeature> troughs;
};

// SVG panel dimensions and geographic window. Europe-centred, wide
// enough to match the FMI chart crop.
constexpr double kLonMin = -30.0;
constexpr double kLonMax = 45.0;
constexpr double kLatMin = 30.0;
constexpr double kLatMax = 72.0;
constexpr int kPanelW = 900;
constexpr int kPanelH = 600;

double xOf(double lon)
{
  return (lon - kLonMin) / (kLonMax - kLonMin) * kPanelW;
}
double yOf(double lat)
{
  return (kLatMax - lat) / (kLatMax - kLatMin) * kPanelH;
}

void svgPath(std::ostream& os,
             const GridView& gv,
             const std::vector<Fmi::Dynlib::FeaturePoint>& pts,
             const char* cssClass)
{
  if (pts.size() < 2) return;
  os << "<path class=\"" << cssClass << "\" d=\"";
  bool moved = false;
  double prev_lon = 0;
  for (std::size_t k = 0; k < pts.size(); ++k)
  {
    double lon, lat;
    toLonLat(gv, pts[k].j, pts[k].i, lon, lat);
    // Break segments that wrap across the antimeridian so the SVG
    // doesn't draw a long diagonal across the panel.
    const bool jump = moved && std::fabs(lon - prev_lon) > 30.0;
    if (!moved || jump)
    {
      os << "M" << xOf(lon) << "," << yOf(lat) << " ";
      moved = true;
    }
    else
      os << "L" << xOf(lon) << "," << yOf(lat) << " ";
    prev_lon = lon;
  }
  os << "\"/>\n";
}

void writeSvgPanel(std::ostream& os,
                   int panel_index,
                   const GridView& gv,
                   const FeatureBundle& fb)
{
  const int ox = panel_index * kPanelW;
  os << "<g transform=\"translate(" << ox << ",0)\">\n";
  os << "<rect width=\"" << kPanelW << "\" height=\"" << kPanelH
     << "\" class=\"panel-bg\"/>\n";

  // Graticule: 10-degree meridians and parallels.
  os << "<g class=\"graticule\">\n";
  for (int lon = -40; lon <= 50; lon += 10)
    os << "<line x1=\"" << xOf(lon) << "\" y1=\"0\" x2=\"" << xOf(lon)
       << "\" y2=\"" << kPanelH << "\"/>\n";
  for (int lat = 30; lat <= 80; lat += 10)
    os << "<line x1=\"0\" y1=\"" << yOf(lat) << "\" x2=\"" << kPanelW
       << "\" y2=\"" << yOf(lat) << "\"/>\n";
  os << "</g>\n";

  // City landmarks for geographic orientation (no coastline dataset).
  struct City { const char* name; double lon; double lat; };
  const City cities[] = {
      {"London", -0.1, 51.5},     {"Paris", 2.3, 48.9},
      {"Berlin", 13.4, 52.5},     {"Madrid", -3.7, 40.4},
      {"Rome", 12.5, 41.9},       {"Moscow", 37.6, 55.8},
      {"Helsinki", 24.9, 60.2},   {"Stockholm", 18.1, 59.3},
      {"Reykjavik", -21.9, 64.1}, {"Istanbul", 29.0, 41.0},
      {"Kyiv", 30.5, 50.4},       {"Warsaw", 21.0, 52.2}};
  os << "<g class=\"city\">\n";
  for (const auto& c : cities)
  {
    const double cx = xOf(c.lon), cy = yOf(c.lat);
    if (cx >= 0 && cx <= kPanelW && cy >= 0 && cy <= kPanelH)
    {
      os << "<circle cx=\"" << cx << "\" cy=\"" << cy << "\" r=\"2\"/>\n";
      os << "<text x=\"" << (cx + 4) << "\" y=\"" << (cy - 4) << "\">" << c.name
         << "</text>\n";
    }
  }
  os << "</g>\n";

  // Detected features. Fronts + jet axes + trough axes are shown;
  // convergence/deformation lines are in the GeoJSON but omitted here
  // to keep the composite comparable to the meteorologist's chart,
  // which doesn't draw them.
  for (const auto& l : fb.fronts)
  {
    const char* cls =
        l.type == Fmi::Dynlib::FrontType::Cold      ? "front-cold"
        : l.type == Fmi::Dynlib::FrontType::Warm    ? "front-warm"
                                                    : "front-stat";
    svgPath(os, gv, l.points, cls);
  }
  for (const auto& l : fb.jets_ff30) svgPath(os, gv, l.points, "jet-ff30");
  for (const auto& l : fb.troughs)   svgPath(os, gv, l.points, "trough");

  // Title.
  os << "<text x=\"10\" y=\"24\" class=\"title\">" << fb.stamp << "</text>\n";
  os << "</g>\n";
}

void writeSvgComposite(const std::string& path,
                       const GridView& gv,
                       const std::vector<FeatureBundle>& bundles)
{
  std::ofstream os(path);
  const int total_w = kPanelW * static_cast<int>(bundles.size());
  const int total_h = kPanelH + 60;  // room for legend row
  os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
     << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
     << "width=\"" << total_w << "\" height=\"" << total_h << "\" "
     << "viewBox=\"0 0 " << total_w << " " << total_h << "\">\n"
     << "<defs><style>\n"
     << " .panel-bg { fill:#f8f8f8; stroke:#888; stroke-width:1 }\n"
     << " .graticule line { stroke:#ccc; stroke-width:0.5 }\n"
     << " .city circle { fill:#666 }\n"
     << " .city text { font:10px sans-serif; fill:#444 }\n"
     << " .title { font:14px sans-serif; fill:#222 }\n"
     << " path { fill:none; stroke-linecap:round; stroke-linejoin:round }\n"
     << " .front-cold { stroke:#00b8e6; stroke-width:2.5 }\n"   // cyan (FMI convention)
     << " .front-warm { stroke:#e61c1c; stroke-width:2.5 }\n"   // red  (FMI convention)
     << " .front-stat { stroke:#905020; stroke-width:2; stroke-dasharray:6,3 }\n"
     << " .jet-ff30    { stroke:#008040; stroke-width:2 }\n"
     << " .jet-default { stroke:#00a060; stroke-width:1; stroke-dasharray:4,3 }\n"
     << " .conv        { stroke:#a000a0; stroke-width:1 }\n"
     << " .defline     { stroke:#505050; stroke-width:1; stroke-dasharray:2,2 }\n"
     << " .trough      { stroke:#000000; stroke-width:1.5; stroke-dasharray:8,3 }\n"
     << " .legend-swatch { stroke-width:3 }\n"
     << " .legend-text   { font:12px sans-serif; fill:#222 }\n"
     << "</style></defs>\n";

  for (std::size_t k = 0; k < bundles.size(); ++k)
    writeSvgPanel(os, static_cast<int>(k), gv, bundles[k]);

  // Legend row below the panels.
  struct LegendItem
  {
    const char* label;
    const char* cls;
  };
  const LegendItem items[] = {
      {"cold front", "front-cold"}, {"warm front", "front-warm"},
      {"stationary front", "front-stat"},
      {"jet axis (|V|>=30 m/s)", "jet-ff30"},
      {"trough axis", "trough"}};
  int lx = 10;
  const int ly = kPanelH + 35;
  os << "<g>\n";
  for (const auto& it : items)
  {
    os << "<line class=\"legend-swatch " << it.cls << "\" x1=\"" << lx
       << "\" y1=\"" << ly << "\" x2=\"" << (lx + 30) << "\" y2=\"" << ly << "\"/>\n";
    os << "<text class=\"legend-text\" x=\"" << (lx + 35) << "\" y=\""
       << (ly + 4) << "\">" << it.label << "</text>\n";
    lx += 35 + 6 * static_cast<int>(std::strlen(it.label)) + 20;
  }
  os << "</g>\n</svg>\n";
}

}  // namespace

int main(int argc, char** argv)
{
  const std::string path =
      (argc > 1 ? argv[1] : std::string(std::getenv("HOME"))
                                + "/hub/europe_pressurelevels_12utc.sqd");

  NFmiQueryData qdata(path);
  NFmiFastQueryInfo info(&qdata);

  std::cerr << "Opened " << path << "\n";
  std::cerr << "Grid: " << info.GridXNumber() << " x " << info.GridYNumber()
            << ", " << info.SizeTimes() << " time(s)\n";

  // Build grid view once (geometry does not change across timesteps).
  info.First();
  const GridView gv = buildGridView(info);

  std::vector<FeatureBundle> bundles;

  for (info.ResetTime(); info.NextTime();)
  {
    const NFmiMetTime t = info.ValidTime();
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d_%02dZ",
                  t.GetYear(), t.GetMonth(), t.GetDay(), t.GetHour());
    const std::string stamp = buf;
    const std::string outpath = "/tmp/dynlib_demo_" + stamp + ".geojson";
    std::cerr << "\n=== " << stamp << " -> " << outpath << " ===\n";

    FeatureBundle bundle;
    bundle.stamp = stamp;

    std::ofstream fout(outpath);
    fout << "{\"type\":\"FeatureCollection\",\"features\":[\n";
    bool first = true;

    Fmi::Matrix<double> field, u, v;

    // Fronts at 850 hPa. Prefer a stored θe field; fall back to a
    // Bolton-derived θe from T + RH; last-resort fall back to dry θ.
    const double p_hPa = 850.0;
    std::string theta_src;
    double miss_f = 1.0;
    if (select(info, kFmiPseudoAdiabaticPotentialTemperature, 850))
    {
      miss_f = extractField(info, field);
      theta_src = "θe (stored)";
    }
    if (miss_f > 0.5)
    {
      // Derive θe from T + RH at this level.
      Fmi::Matrix<double> T_K, RH;
      double miss_T = 1.0, miss_R = 1.0;
      if (select(info, kFmiTemperature, 850)) miss_T = extractField(info, T_K);
      if (select(info, kFmiHumidity, 850)) miss_R = extractField(info, RH);
      if (miss_T < 0.5 && miss_R < 0.5)
      {
        const int nx = T_K.width(), ny = T_K.height();
        field = Fmi::Matrix<double>(nx, ny);
        for (int j = 0; j < ny; ++j)
          for (int i = 0; i < nx; ++i)
          {
            // Querydata stores temperature in Celsius.
            const double Tk = T_K(i, j) + 273.15;
            field(i, j) = boltonThetaE(Tk, p_hPa, RH(i, j));
          }
        miss_f = 0.0;
        theta_src = "θe (Bolton from T+RH)";
      }
    }
    if (miss_f > 0.5 && select(info, kFmiPotentialTemperature, 850))
    {
      miss_f = extractField(info, field);
      theta_src = "θ (dry, fallback)";
    }
    if (miss_f < 0.5)
    {
      if (select(info, kFmiWindUMS, 850)) extractField(info, u);
      if (select(info, kFmiWindVMS, 850)) extractField(info, v);
      // Upstream defaults were tuned for ERA-Interim (~1.5°). At this
      // file's ~0.25° resolution we need more smoothing and a positive
      // |∇field| threshold to suppress noisy short segments. |∇θe| at
      // frontal zones is typically several 1e-5 K/m, so 1.5e-5 lets
      // the main baroclinic zones through while dropping weak noise.
      Fmi::Dynlib::FrontOptions opts;
      opts.smoothing_passes = 8;
      opts.intensity_threshold = 1.5e-5;
      const auto lines = Fmi::Dynlib::detectFrontsMaxGrad(field, u, v, gv.dx, gv.dy, opts);
      std::size_t pts = 0;
      for (const auto& l : lines) pts += l.points.size();
      std::cerr << "  fronts_maxgrad(850 " << theta_src << "): "
                << lines.size() << " lines / " << pts << " pts\n";
      for (const auto& l : lines)
        writeLineString(fout, first, gv, l.points, "fronts_maxgrad",
                        frontTypeName(l.type));
      bundle.fronts = lines;
    }
    else
    {
      std::cerr << "  fronts: no usable scalar field at 850 hPa (all NaN)\n";
    }

    // Jet axes at 300 hPa.
    if (select(info, kFmiWindUMS, 300))
    {
      extractField(info, u);
      if (select(info, kFmiWindVMS, 300)) extractField(info, v);
      Fmi::Dynlib::LineOptions opts;
      opts.smoothing_passes = 2;
      const auto jets = Fmi::Dynlib::detectJetAxes(u, v, gv.dx, gv.dy, opts);
      const auto jets30 = Fmi::Dynlib::detectJetAxesFFThres(u, v, gv.dx, gv.dy, opts);
      std::cerr << "  jet_axes(300): " << jets.size()
                << "  jet_axes_ff>=30(300): " << jets30.size() << "\n";
      for (const auto& l : jets)
        writeLineString(fout, first, gv, l.points, "jet_axes");
      for (const auto& l : jets30)
        writeLineString(fout, first, gv, l.points, "jet_axes_ff30");
      bundle.jets_default = jets;
      bundle.jets_ff30 = jets30;
    }

    // Low-level convergence / deformation / vorticity (trough) lines at 1000 hPa.
    if (select(info, kFmiWindUMS, 1000))
    {
      extractField(info, u);
      if (select(info, kFmiWindVMS, 1000)) extractField(info, v);
      const auto conv = Fmi::Dynlib::detectConvergenceLines(u, v, gv.dx, gv.dy);
      const auto defl = Fmi::Dynlib::detectDeformationLines(u, v, gv.dx, gv.dy);
      const auto vort = Fmi::Dynlib::detectVorticityLines(u, v, gv.dx, gv.dy);
      std::cerr << "  conv(1000): " << conv.size() << "  def(1000): " << defl.size()
                << "  vor/trough(1000): " << vort.size() << "\n";
      for (const auto& l : conv)
        writeLineString(fout, first, gv, l.points, "convergence_lines");
      for (const auto& l : defl)
        writeLineString(fout, first, gv, l.points, "deformation_lines");
      for (const auto& l : vort)
        writeLineString(fout, first, gv, l.points, "trough_axes");
      bundle.convergence = conv;
      bundle.deformation = defl;
      bundle.troughs = vort;
    }

    // Blocking indicator on 500 hPa geopotential height. Scalar field;
    // just report statistics, not a GeoJSON feature.
    if (select(info, kFmiGeopHeight, 500))
    {
      extractField(info, field);
      const auto bi = Fmi::Dynlib::blockingIndicator(field, gv.dx, gv.dy);
      double bi_max = -1e30;
      std::size_t pos = 0;
      for (int j = 0; j < bi.height(); ++j)
        for (int i = 0; i < bi.width(); ++i)
        {
          if (bi(i, j) > bi_max) bi_max = bi(i, j);
          if (bi(i, j) > 0) ++pos;
        }
      std::cerr << "  blocking(500): max=" << bi_max
                << "  positive_cells=" << pos << " / "
                << (bi.width() * bi.height()) << "\n";
    }

    fout << "\n]}\n";
    fout.close();
    bundles.push_back(std::move(bundle));
  }

  const std::string svg_path = "/tmp/dynlib_demo_composite.svg";
  writeSvgComposite(svg_path, gv, bundles);
  std::cerr << "\nComposite SVG: " << svg_path << "\n";
  std::cerr << "Done.\n";
  return 0;
}
