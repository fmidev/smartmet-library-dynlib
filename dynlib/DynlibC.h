#pragma once

#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

// ----------------------------------------------------------------------
// Array conventions
// ----------------------------------------------------------------------
// All input arrays are row-major (i.e. C order): element (i, j) of an
// (nx, ny) field is at index i + nx*j. dx/dy are the "double grid
// spacing" in metres — dx(i,j) is the distance between (i-1, j) and
// (i+1, j), analogous for dy. For a regular WGS84 lon/lat grid that
// is typically 2*R*cos(lat)*dlon and 2*R*dlat.
//
// Line-detector outputs
// ---------------------
// Point buffers hold triples (j, i, value). Unused slots are NaN.
// Offset buffers hold the start index of each line as a 0-based offset
// into the point buffer. Unused slots are NaN. The final written slot
// acts as a sentinel equal to the total number of points (so offsets
// form a monotonically non-decreasing sequence; the decoder stops at
// the first NaN or the first non-monotone value).
//
// Three-type front output (fronts_maxgrad / fronts_maxcurv)
//   pts: 3 * no * 3 doubles indexed pts[t*no*3 + p*3 + c]
//   off: 3 * nf     doubles indexed off[t*nf + f]
//   t=0 cold, t=1 warm, t=2 stationary.
//
// Single-list output (jet_axes / lines)
//   pts: no * 3 doubles indexed pts[p*3 + c]
//   off: nf     doubles indexed off[f]
// ----------------------------------------------------------------------

// Jenkner et al. (2010) max-gradient front detection.
void dynlib_detect_fronts_maxgrad(int32_t nx, int32_t ny, int32_t no, int32_t nf,
                                  const double* field, const double* u, const double* v,
                                  const double* dx, const double* dy,
                                  double* pts_out, double* off_out,
                                  double frint_thres, double frspd_thres,
                                  int32_t nsmooth);

// Berry et al. (2007) max-curvature front detection. Same API shape.
void dynlib_detect_fronts_maxcurv(int32_t nx, int32_t ny, int32_t no, int32_t nf,
                                  const double* field, const double* u, const double* v,
                                  const double* dx, const double* dy,
                                  double* pts_out, double* off_out,
                                  double frint_thres, double frspd_thres,
                                  int32_t nsmooth);

// Jet axis detection (single list, no type split).
// variant: 0 = jetaxis (shear second-derivative mask; default),
//          1 = jetaxis_ff_thres (same mask + fixed 30 m/s wind-speed cutoff).
void dynlib_detect_jet_axes(int32_t variant,
                            int32_t nx, int32_t ny, int32_t no, int32_t nf,
                            const double* u, const double* v,
                            const double* dx, const double* dy,
                            double* pts_out, double* off_out,
                            int32_t nsmooth);

// Convergence / deformation / vorticity line detection.
// kind_code: 1 = convergence, 2 = deformation, 3 = vorticity (trough proxy).
void dynlib_detect_lines(int32_t kind_code,
                         int32_t nx, int32_t ny, int32_t no, int32_t nf,
                         const double* u, const double* v,
                         const double* dx, const double* dy,
                         double* pts_out, double* off_out,
                         int32_t nsmooth);

// Masato et al. (2012) blocking indicator (gradient-reversal variant).
// Returns a field of the same shape as the input; positive values
// mark candidate blocks.
void dynlib_block_indicator(int32_t nx, int32_t ny,
                            const double* field,
                            const double* dx, const double* dy,
                            double* out_field);

// Rossby wave breaking by local y-gradient reversal (Spensberger).
// Inputs: scalar field `pv` (PV on an isentrope, or any similar
// field), a per-cell `mask` where >=0.5 enables the test, per-row
// `latitudes` (length ny), and a y-gradient threshold. Outputs are
// seven same-shape double fields (written by caller; not allocated):
//   resa    : anticyclonic flag (0 or 1) — threshold-screened
//   resc    : cyclonic flag (0 or 1)
//   resai   : |grad(field)| at significant anticyclonic reversals
//   resci   : |grad(field)| at significant cyclonic reversals
//   resaiy  : |dfield/dy| of all anticyclonic reversals (unthresholded)
//   resciy  : |dfield/dy| of all cyclonic reversals (unthresholded)
//   tested  : flag (0 or 1) for each point actually tested.
void dynlib_detect_rwb_grad_rev(int32_t nx, int32_t ny,
                                const double* pv,
                                const double* mask,
                                const double* latitudes,
                                double ddy_thres,
                                const double* dx, const double* dy,
                                double* resa_out, double* resc_out,
                                double* resai_out, double* resci_out,
                                double* resaiy_out, double* resciy_out,
                                double* tested_out);

// Cyclone detection (Wernli & Schwierz contour method). The caller
// must pre-sort MSL values ascending and supply the sorted values
// along with the original (i, j) indices as iis/jjs (0-based, length
// nx*ny). Config overrides are in native units — pass negative to
// keep upstream defaults (km^2, m, km, Pa).
//   mask_out: same shape as msl; 0 = no cyclone, 1..N = cyclone id.
//   meta_out: 5 * nn doubles indexed meta[m*5 + c] where
//             c=0 lat, 1 lon, 2 min SLP, 3 outer SLP, 4 size (km^2).
//             Unused cyclone slots are all-zero.
void dynlib_detect_cyclones(int32_t nx, int32_t ny, int32_t nn,
                            const double* msl,
                            const double* msl_sorted,
                            const int32_t* iis,
                            const int32_t* jjs,
                            const double* oro,
                            const double* lon,
                            const double* lat,
                            const double* dx, const double* dy,
                            double cyc_minsize_km2,
                            double cyc_maxsize_km2,
                            double cyc_maxoro_m,
                            double cyc_mindist_km,
                            double cyc_minprominence,
                            double* mask_out,
                            double* meta_out);

// Precipitation blob detection. Same calling protocol as
// dynlib_detect_cyclones, but without orography and with a single
// user-provided minimum inter-blob distance (km). The sorted array
// carries sign-flipped precipitation values (upstream triggers on
// val < 0), so pass -precip and its ascending sort order. precip
// itself is passed unmodified.
void dynlib_detect_blobs(int32_t nx, int32_t ny, int32_t nn,
                         const double* precip,
                         const double* precip_sorted_neg,
                         const int32_t* iis,
                         const int32_t* jjs,
                         const double* lon,
                         const double* lat,
                         const double* dx, const double* dy,
                         double blob_mindist_km,
                         double* mask_out,
                         double* meta_out);

#ifdef __cplusplus
}  // extern "C"
#endif
