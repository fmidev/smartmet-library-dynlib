! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ISO_C_BINDING shim over the dynlib feature-detection subroutines.
!
! Array-order convention
! ----------------------
! dynlib's Fortran subroutines declare 2D fields as (ny, nx) and scalar
! 3D fields as (nz, ny, nx). Under Fortran column-major storage this
! places the j (y) index in the fastest-varying position for a 2D slice.
! The C callers in this library supply arrays laid out in row-major
! order with the i (x) index fastest (matching Fmi::Matrix<double>).
! The shim therefore transposes inputs on entry.
!
! Output packing
! --------------
! dynlib returns line features as two arrays: a flat point array of
! (j, i, value) triples and a separate offset array whose entries mark
! the start index of each line within the point array. The upstream
! only writes the valid prefix; unused slots are left uninitialised.
! To make the output unambiguous for the C decoder, the shim
! initialises both arrays to NaN before each call AND re-copies the
! upstream result into the C output buffers, so NaN in an offset slot
! means "no more lines".
!
! The repacking also transposes the output array axes from Fortran
! column-major layout (type fastest, component slowest) to C
! row-major layout (component fastest, type slowest). The C decoder
! reads pts[t*no*3 + p*3 + c].
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module dynlib_wrapper
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_int32_t
  implicit none
  private
  public :: dynlib_detect_fronts_maxgrad_c
  public :: dynlib_detect_fronts_maxcurv_c
  public :: dynlib_detect_jet_axes_c
  public :: dynlib_detect_lines_c
  public :: dynlib_block_indicator_c
  public :: dynlib_detect_rwb_grad_rev_c
  public :: dynlib_detect_cyclones_c
  public :: dynlib_detect_blobs_c

contains

  ! -------------------------------------------------------------------
  ! Internal helper: set thresholds / smoothing from caller. A negative
  ! value means "keep the upstream default".
  ! -------------------------------------------------------------------
  subroutine apply_config(frint_thres, frspd_thres, nsmooth_in)
    use kind,   only: nr, ni
    use config, only: frint_thres_cfg => frint_thres, &
                      frspd_thres_cfg => frspd_thres, &
                      nsmooth
    real(c_double), intent(in) :: frint_thres, frspd_thres
    integer(c_int), intent(in) :: nsmooth_in

    ! Upstream defaults from dynlib's config.f90. Reset on every call
    ! so that a previous caller's override does not leak state into
    ! the current detection. Callers pass -1 to keep the default;
    ! non-negative values override.
    real(kind=nr), parameter :: default_frint = -2.0e-11_nr
    real(kind=nr), parameter :: default_frspd = 1.5_nr
    integer(ni),   parameter :: default_nsmooth = 2_ni

    frint_thres_cfg = default_frint
    frspd_thres_cfg = default_frspd
    nsmooth         = default_nsmooth

    if (frint_thres >= 0.0_c_double) frint_thres_cfg = real(frint_thres, nr)
    if (frspd_thres >= 0.0_c_double) frspd_thres_cfg = real(frspd_thres, nr)
    if (nsmooth_in  >= 0_c_int)      nsmooth         = int(nsmooth_in,  ni)
  end subroutine apply_config

  ! -------------------------------------------------------------------
  ! Repack a 3-type line result (used by front detection) from the
  ! Fortran buffers fr(1, 3, no, 3), froff(1, 3, nf) into the C output
  ! buffers pts_out(3, no, 3) [indexed (c, p, t)] and off_out(nf, 3).
  ! -------------------------------------------------------------------
  subroutine repack_typed_lines(fr, froff, nx_unused, no, nf, pts_out, off_out)
    use kind,   only: nr, ni
    use consts, only: nan
    real(kind=nr), intent(in)  :: fr(1, 3_ni, no, 3_ni)
    real(kind=nr), intent(in)  :: froff(1, 3_ni, nf)
    integer(c_int), intent(in) :: nx_unused, no, nf
    real(c_double), intent(out) :: pts_out(3, no, 3)
    real(c_double), intent(out) :: off_out(nf, 3)
    integer(ni) :: t, p
    pts_out = nan
    off_out = nan
    do t = 1, 3
      do p = 1, no
        pts_out(1, p, t) = fr(1, t, p, 1)   ! j
        pts_out(2, p, t) = fr(1, t, p, 2)   ! i
        pts_out(3, p, t) = fr(1, t, p, 3)   ! value
      end do
      off_out(:, t) = froff(1, t, :)
    end do
  end subroutine repack_typed_lines

  ! -------------------------------------------------------------------
  ! Repack a single-type line result (jet axes / conv / def / vor /
  ! maxcurv) from fr(1, no, 3), froff(1, nf) into pts_out(3, no) and
  ! off_out(nf). The C decoder reads pts[p*3 + c] and off[f].
  ! -------------------------------------------------------------------
  subroutine repack_lines(fr, froff, no, nf, pts_out, off_out)
    use kind,   only: nr, ni
    use consts, only: nan
    real(kind=nr), intent(in)  :: fr(1, no, 3_ni)
    real(kind=nr), intent(in)  :: froff(1, nf)
    integer(c_int), intent(in) :: no, nf
    real(c_double), intent(out) :: pts_out(3, no)
    real(c_double), intent(out) :: off_out(nf)
    integer(ni) :: p
    pts_out = nan
    off_out = nan
    do p = 1, no
      pts_out(1, p) = fr(1, p, 1)   ! j
      pts_out(2, p) = fr(1, p, 2)   ! i
      pts_out(3, p) = fr(1, p, 3)   ! value
    end do
    off_out(:) = froff(1, :)
  end subroutine repack_lines

  ! ===================================================================
  ! Public entry points
  ! ===================================================================

  subroutine dynlib_detect_fronts_maxgrad_c( &
      nx, ny, no, nf, field, u, v, dx, dy, &
      pts_out, off_out, frint_thres, frspd_thres, nsmooth_in) &
      bind(c, name="dynlib_detect_fronts_maxgrad")
    use kind,   only: nr, ni
    use consts, only: nan
    use detect, only: frontline_at_maxgrad

    integer(c_int), value, intent(in)  :: nx, ny, no, nf, nsmooth_in
    real(c_double),        intent(in)  :: field(nx, ny), u(nx, ny), v(nx, ny)
    real(c_double),        intent(in)  :: dx(nx, ny), dy(nx, ny)
    real(c_double),        intent(out) :: pts_out(3, no, 3)
    real(c_double),        intent(out) :: off_out(nf, 3)
    real(c_double), value, intent(in)  :: frint_thres, frspd_thres

    real(kind=nr) :: field_f(1, ny, nx), u_f(1, ny, nx), v_f(1, ny, nx)
    real(kind=nr) :: dx_f(ny, nx), dy_f(ny, nx)
    real(kind=nr) :: fr(1, 3_ni, no, 3_ni), froff(1, 3_ni, nf)

    call apply_config(frint_thres, frspd_thres, nsmooth_in)

    field_f(1,:,:) = transpose(field)
    u_f    (1,:,:) = transpose(u)
    v_f    (1,:,:) = transpose(v)
    dx_f          = transpose(dx)
    dy_f          = transpose(dy)

    fr    = nan
    froff = nan

    call frontline_at_maxgrad(fr, froff, int(nx, ni), int(ny, ni), 1_ni, &
                              int(no, ni), int(nf, ni), field_f, u_f, v_f, dx_f, dy_f)

    call repack_typed_lines(fr, froff, nx, no, nf, pts_out, off_out)
  end subroutine dynlib_detect_fronts_maxgrad_c

  subroutine dynlib_detect_fronts_maxcurv_c( &
      nx, ny, no, nf, field, u, v, dx, dy, &
      pts_out, off_out, frint_thres, frspd_thres, nsmooth_in) &
      bind(c, name="dynlib_detect_fronts_maxcurv")
    use kind,   only: nr, ni
    use consts, only: nan
    use detect, only: frontline_at_maxcurv

    integer(c_int), value, intent(in)  :: nx, ny, no, nf, nsmooth_in
    real(c_double),        intent(in)  :: field(nx, ny), u(nx, ny), v(nx, ny)
    real(c_double),        intent(in)  :: dx(nx, ny), dy(nx, ny)
    real(c_double),        intent(out) :: pts_out(3, no, 3)
    real(c_double),        intent(out) :: off_out(nf, 3)
    real(c_double), value, intent(in)  :: frint_thres, frspd_thres

    real(kind=nr) :: field_f(1, ny, nx), u_f(1, ny, nx), v_f(1, ny, nx)
    real(kind=nr) :: dx_f(ny, nx), dy_f(ny, nx)
    real(kind=nr) :: fr(1, 3_ni, no, 3_ni), froff(1, 3_ni, nf)

    call apply_config(frint_thres, frspd_thres, nsmooth_in)

    field_f(1,:,:) = transpose(field)
    u_f    (1,:,:) = transpose(u)
    v_f    (1,:,:) = transpose(v)
    dx_f          = transpose(dx)
    dy_f          = transpose(dy)

    fr    = nan
    froff = nan

    call frontline_at_maxcurv(fr, froff, int(nx, ni), int(ny, ni), 1_ni, &
                              int(no, ni), int(nf, ni), field_f, u_f, v_f, dx_f, dy_f)

    call repack_typed_lines(fr, froff, nx, no, nf, pts_out, off_out)
  end subroutine dynlib_detect_fronts_maxcurv_c

  ! Jet axes. Inputs: wind components + grid spacing. No scalar field,
  ! no per-type output. Uses the same point / offset packing as the
  ! other single-list line detectors. variant: 0 = jetaxis (shear
  ! second-derivative mask), 1 = jetaxis_ff_thres (shear mask + fixed
  ! 30 m/s wind speed threshold).
  subroutine dynlib_detect_jet_axes_c( &
      variant, nx, ny, no, nf, u, v, dx, dy, pts_out, off_out, nsmooth_in) &
      bind(c, name="dynlib_detect_jet_axes")
    use kind,   only: nr, ni
    use consts, only: nan
    use detect, only: jetaxis, jetaxis_ff_thres

    integer(c_int), value, intent(in)  :: variant, nx, ny, no, nf, nsmooth_in
    real(c_double),        intent(in)  :: u(nx, ny), v(nx, ny)
    real(c_double),        intent(in)  :: dx(nx, ny), dy(nx, ny)
    real(c_double),        intent(out) :: pts_out(3, no)
    real(c_double),        intent(out) :: off_out(nf)

    real(kind=nr) :: u_f(1, ny, nx), v_f(1, ny, nx)
    real(kind=nr) :: dx_f(ny, nx), dy_f(ny, nx)
    real(kind=nr) :: ja(1, no, 3_ni), jaoff(1, nf)

    call apply_config(-1.0_c_double, -1.0_c_double, nsmooth_in)

    u_f (1,:,:) = transpose(u)
    v_f (1,:,:) = transpose(v)
    dx_f       = transpose(dx)
    dy_f       = transpose(dy)

    ja    = nan
    jaoff = nan

    if (variant == 1_c_int) then
      call jetaxis_ff_thres(ja, jaoff, int(nx, ni), int(ny, ni), 1_ni, &
                            int(no, ni), int(nf, ni), u_f, v_f, dx_f, dy_f)
    else
      call jetaxis(ja, jaoff, int(nx, ni), int(ny, ni), 1_ni, &
                   int(no, ni), int(nf, ni), u_f, v_f, dx_f, dy_f)
    end if

    call repack_lines(ja, jaoff, no, nf, pts_out, off_out)
  end subroutine dynlib_detect_jet_axes_c

  ! Convergence / deformation / vorticity line detection. Same inputs
  ! and output shape as jet axes. The kind_code selects which detector:
  !   1 = convline, 2 = defline, 3 = vorline.
  subroutine dynlib_detect_lines_c( &
      kind_code, nx, ny, no, nf, u, v, dx, dy, pts_out, off_out, nsmooth_in) &
      bind(c, name="dynlib_detect_lines")
    use kind,   only: nr, ni
    use consts, only: nan
    use detect, only: convline, defline, vorline

    integer(c_int), value, intent(in)  :: kind_code, nx, ny, no, nf, nsmooth_in
    real(c_double),        intent(in)  :: u(nx, ny), v(nx, ny)
    real(c_double),        intent(in)  :: dx(nx, ny), dy(nx, ny)
    real(c_double),        intent(out) :: pts_out(3, no)
    real(c_double),        intent(out) :: off_out(nf)

    real(kind=nr) :: u_f(1, ny, nx), v_f(1, ny, nx)
    real(kind=nr) :: dx_f(ny, nx), dy_f(ny, nx)
    real(kind=nr) :: fr(1, no, 3_ni), froff(1, nf)

    call apply_config(-1.0_c_double, -1.0_c_double, nsmooth_in)

    u_f (1,:,:) = transpose(u)
    v_f (1,:,:) = transpose(v)
    dx_f       = transpose(dx)
    dy_f       = transpose(dy)

    fr    = nan
    froff = nan

    select case (kind_code)
    case (1)
      call convline(fr, froff, int(nx, ni), int(ny, ni), 1_ni, &
                    int(no, ni), int(nf, ni), u_f, v_f, dx_f, dy_f)
    case (2)
      call defline (fr, froff, int(nx, ni), int(ny, ni), 1_ni, &
                    int(no, ni), int(nf, ni), u_f, v_f, dx_f, dy_f)
    case default
      call vorline (fr, froff, int(nx, ni), int(ny, ni), 1_ni, &
                    int(no, ni), int(nf, ni), u_f, v_f, dx_f, dy_f)
    end select

    call repack_lines(fr, froff, no, nf, pts_out, off_out)
  end subroutine dynlib_detect_lines_c

  ! Blocking indicator (Masato et al. 2012, gradient-reversal variant).
  ! Returns a field of the same shape as the input (the indicator
  ! value; positive means a candidate block). Caller contours / masks
  ! it downstream.
  subroutine dynlib_block_indicator_c( &
      nx, ny, field, dx, dy, out_field) &
      bind(c, name="dynlib_block_indicator")
    use kind,   only: nr, ni
    use detect, only: block_indicator_grad_rev

    integer(c_int), value, intent(in)  :: nx, ny
    real(c_double),        intent(in)  :: field(nx, ny), dx(nx, ny), dy(nx, ny)
    real(c_double),        intent(out) :: out_field(nx, ny)

    real(kind=nr) :: field_f(1, ny, nx), dx_f(ny, nx), dy_f(ny, nx)
    real(kind=nr) :: res_f(1, ny, nx)

    field_f(1,:,:) = transpose(field)
    dx_f          = transpose(dx)
    dy_f          = transpose(dy)

    call block_indicator_grad_rev(res_f, int(nx, ni), int(ny, ni), 1_ni, &
                                  field_f, dx_f, dy_f)

    out_field = transpose(res_f(1, :, :))
  end subroutine dynlib_block_indicator_c

  ! Rossby wave breaking detection by local y-gradient reversal
  ! (Spensberger gradient-reversal method). Inputs are the scalar
  ! field (PV or similar), a user mask (1 = test, 0 = skip), per-row
  ! latitudes, the y-gradient threshold, and grid spacing. Outputs
  ! are seven same-shape double fields — Fortran uses int8 internally
  ! but the shim casts to double so the C side can reuse the regular
  ! Matrix<double> plumbing.
  subroutine dynlib_detect_rwb_grad_rev_c( &
      nx, ny, pv, mask, latitudes, ddy_thres, dx, dy, &
      resa_out, resc_out, resai_out, resci_out, resaiy_out, resciy_out, &
      tested_out) &
      bind(c, name="dynlib_detect_rwb_grad_rev")
    use kind,   only: nr, ni
    use detect, only: rwb_by_grad_rev

    integer(c_int), value, intent(in)  :: nx, ny
    real(c_double),        intent(in)  :: pv(nx, ny), mask(nx, ny)
    real(c_double),        intent(in)  :: latitudes(ny)
    real(c_double), value, intent(in)  :: ddy_thres
    real(c_double),        intent(in)  :: dx(nx, ny), dy(nx, ny)
    real(c_double),        intent(out) :: resa_out(nx, ny), resc_out(nx, ny)
    real(c_double),        intent(out) :: resai_out(nx, ny), resci_out(nx, ny)
    real(c_double),        intent(out) :: resaiy_out(nx, ny), resciy_out(nx, ny)
    real(c_double),        intent(out) :: tested_out(nx, ny)

    real(kind=nr) :: pv_f(1, ny, nx), dx_f(ny, nx), dy_f(ny, nx)
    real(kind=nr) :: lats_f(ny)
    integer(kind=1) :: mask_f(1, ny, nx)
    integer(kind=1) :: resa_f(1, ny, nx), resc_f(1, ny, nx), tested_f(1, ny, nx)
    real(kind=nr) :: resai_f(1, ny, nx), resci_f(1, ny, nx)
    real(kind=nr) :: resaiy_f(1, ny, nx), resciy_f(1, ny, nx)
    integer(ni) :: j, i

    pv_f(1,:,:) = transpose(pv)
    dx_f       = transpose(dx)
    dy_f       = transpose(dy)
    lats_f     = latitudes

    ! mask is a double 0/1 on the C side. Anything >= 0.5 enables the test.
    do i = 1_ni, nx
      do j = 1_ni, ny
        if (mask(i, j) >= 0.5_c_double) then
          mask_f(1, j, i) = 1_1
        else
          mask_f(1, j, i) = 0_1
        end if
      end do
    end do

    call rwb_by_grad_rev(resa_f, resc_f, resai_f, resci_f, resaiy_f, resciy_f, &
                         tested_f, int(nx, ni), int(ny, ni), 1_ni, &
                         pv_f, mask_f, lats_f, real(ddy_thres, nr), dx_f, dy_f)

    resa_out   = transpose(real(resa_f  (1, :, :), c_double))
    resc_out   = transpose(real(resc_f  (1, :, :), c_double))
    resai_out  = transpose(real(resai_f (1, :, :), c_double))
    resci_out  = transpose(real(resci_f (1, :, :), c_double))
    resaiy_out = transpose(real(resaiy_f(1, :, :), c_double))
    resciy_out = transpose(real(resciy_f(1, :, :), c_double))
    tested_out = transpose(real(tested_f(1, :, :), c_double))
  end subroutine dynlib_detect_rwb_grad_rev_c

  ! Cyclone detection (Wernli & Schwierz contour algorithm).
  ! The caller must pre-sort MSL values ascending and pass the sorted
  ! values (msls) together with the original (i, j) grid indices as
  ! (iis, jjs) — both 0-based, length nx*ny. Config overrides for
  ! cyc_minsize/maxsize/maxoro/mindist/minprominence are saved and
  ! restored so a previous caller's values cannot leak state.
  subroutine dynlib_detect_cyclones_c( &
      nx, ny, nn, msl, msls, iis, jjs, oro, lon, lat, dx, dy, &
      minsize_in, maxsize_in, maxoro_in, mindist_in, minprominence_in, &
      mask_out, meta_out) &
      bind(c, name="dynlib_detect_cyclones")
    use kind,   only: nr, ni
    use config, only: cyc_minsize, cyc_maxsize, cyc_maxoro, cyc_mindist, cyc_minprominence
    use detect, only: cyclone_by_contour_fortran

    integer(c_int), value, intent(in)  :: nx, ny, nn
    real(c_double),        intent(in)  :: msl(nx, ny)
    real(c_double),        intent(in)  :: msls(nx * ny)
    integer(c_int32_t),    intent(in)  :: iis(nx * ny), jjs(nx * ny)
    real(c_double),        intent(in)  :: oro(nx, ny)
    real(c_double),        intent(in)  :: lon(nx), lat(ny)
    real(c_double),        intent(in)  :: dx(nx, ny), dy(nx, ny)
    real(c_double), value, intent(in)  :: minsize_in, maxsize_in
    real(c_double), value, intent(in)  :: maxoro_in, mindist_in, minprominence_in
    real(c_double),        intent(out) :: mask_out(nx, ny)
    real(c_double),        intent(out) :: meta_out(5, nn)

    real(kind=nr) :: msl_f(1, ny, nx), oro_f(ny, nx), dx_f(ny, nx), dy_f(ny, nx)
    real(kind=nr) :: msls_f(1, ny * nx)
    integer(kind=ni) :: iis_f(1, ny * nx), jjs_f(1, ny * nx)
    real(kind=nr) :: lon_f(nx), lat_f(ny)
    real(kind=nr) :: meta_f(1, nn, 5_ni)
    integer(kind=ni) :: mask_f(1, ny, nx)
    integer(ni) :: m
    real(kind=nr) :: save_minsize, save_maxsize, save_maxoro, save_mindist, save_minprom

    msl_f(1,:,:) = transpose(msl)
    oro_f       = transpose(oro)
    dx_f        = transpose(dx)
    dy_f        = transpose(dy)
    lon_f       = lon
    lat_f       = lat

    msls_f(1, :) = msls
    iis_f (1, :) = iis
    jjs_f (1, :) = jjs

    ! Save + override cyclone config. Negative input keeps default.
    save_minsize = cyc_minsize
    save_maxsize = cyc_maxsize
    save_maxoro  = cyc_maxoro
    save_mindist = cyc_mindist
    save_minprom = cyc_minprominence
    if (minsize_in       >= 0.0_c_double) cyc_minsize       = real(minsize_in,       nr)
    if (maxsize_in       >= 0.0_c_double) cyc_maxsize       = real(maxsize_in,       nr)
    if (maxoro_in        >= 0.0_c_double) cyc_maxoro        = real(maxoro_in,        nr)
    if (mindist_in       >= 0.0_c_double) cyc_mindist       = real(mindist_in,       nr)
    if (minprominence_in >= 0.0_c_double) cyc_minprominence = real(minprominence_in, nr)

    mask_f = 0_ni

    call cyclone_by_contour_fortran(mask_f, meta_f, &
                                    int(nx, ni), int(ny, ni), 1_ni, int(nn, ni), &
                                    msl_f, msls_f, iis_f, jjs_f, &
                                    oro_f, lon_f, lat_f, dx_f, dy_f)

    ! Restore caller-invisible global config.
    cyc_minsize       = save_minsize
    cyc_maxsize       = save_maxsize
    cyc_maxoro        = save_maxoro
    cyc_mindist       = save_mindist
    cyc_minprominence = save_minprom

    mask_out = transpose(real(mask_f(1, :, :), c_double))
    do m = 1_ni, int(nn, ni)
      meta_out(1, m) = meta_f(1, m, 1_ni)  ! lat
      meta_out(2, m) = meta_f(1, m, 2_ni)  ! lon
      meta_out(3, m) = meta_f(1, m, 3_ni)  ! min SLP
      meta_out(4, m) = meta_f(1, m, 4_ni)  ! outer SLP
      meta_out(5, m) = meta_f(1, m, 5_ni)  ! size
    end do
  end subroutine dynlib_detect_cyclones_c

  ! Precipitation blob detection. Same sort-precomputation protocol as
  ! cyclones; the caller pre-sorts precipitation ascending with the
  ! original (i, j) indices. `precip_sorted` must contain the values
  ! with sign flipped (upstream tests `val < 0` as its trigger for a
  ! "large positive" precipitation maximum). The sign flip is the
  ! caller's responsibility because the upstream code reads the
  ! original precip field unchanged for neighbour lookups.
  subroutine dynlib_detect_blobs_c( &
      nx, ny, nn, precip, precip_sorted, iis, jjs, lon, lat, dx, dy, &
      blob_mindist_km, mask_out, meta_out) &
      bind(c, name="dynlib_detect_blobs")
    use kind,   only: nr, ni
    use detect, only: blobs_fortran

    integer(c_int), value, intent(in)  :: nx, ny, nn
    real(c_double),        intent(in)  :: precip(nx, ny)
    real(c_double),        intent(in)  :: precip_sorted(nx * ny)
    integer(c_int32_t),    intent(in)  :: iis(nx * ny), jjs(nx * ny)
    real(c_double),        intent(in)  :: lon(nx), lat(ny)
    real(c_double),        intent(in)  :: dx(nx, ny), dy(nx, ny)
    real(c_double), value, intent(in)  :: blob_mindist_km
    real(c_double),        intent(out) :: mask_out(nx, ny)
    real(c_double),        intent(out) :: meta_out(5, nn)

    real(kind=nr) :: precip_f(1, ny, nx), dx_f(ny, nx), dy_f(ny, nx)
    real(kind=nr) :: precip_sorted_f(1, ny * nx)
    integer(kind=ni) :: iis_f(1, ny * nx), jjs_f(1, ny * nx)
    real(kind=nr) :: lon_f(nx), lat_f(ny)
    real(kind=nr) :: meta_f(1, nn, 5_ni)
    integer(kind=ni) :: mask_f(1, ny, nx)
    integer(ni) :: m

    precip_f(1,:,:) = transpose(precip)
    dx_f           = transpose(dx)
    dy_f           = transpose(dy)
    lon_f          = lon
    lat_f          = lat

    precip_sorted_f(1, :) = precip_sorted
    iis_f         (1, :) = iis
    jjs_f         (1, :) = jjs

    mask_f = 0_ni

    call blobs_fortran(mask_f, meta_f, &
                       int(nx, ni), int(ny, ni), 1_ni, int(nn, ni), &
                       precip_f, precip_sorted_f, iis_f, jjs_f, &
                       lon_f, lat_f, dx_f, dy_f, real(blob_mindist_km, nr))

    mask_out = transpose(real(mask_f(1, :, :), c_double))
    do m = 1_ni, int(nn, ni)
      meta_out(1, m) = meta_f(1, m, 1_ni)  ! lat
      meta_out(2, m) = meta_f(1, m, 2_ni)  ! lon
      meta_out(3, m) = meta_f(1, m, 3_ni)  ! peak value (flipped sign inside)
      meta_out(4, m) = meta_f(1, m, 4_ni)  ! outermost value
      meta_out(5, m) = meta_f(1, m, 5_ni)  ! size
    end do
  end subroutine dynlib_detect_blobs_c

end module dynlib_wrapper
