MODULE DYNBAS_PIECES_MODULE

  INTERFACE APPEND
		MODULE PROCEDURE APPEND_MATRIX_MATRIX_CHARACTER, APPEND_MATRIX_VECTOR_CHARACTER, &
		                 APPEND_VECTOR_VECTOR_CHARACTER, APPEND_MATRIX_MATRIX_INTEGER, &
		                 APPEND_MATRIX_VECTOR_INTEGER, APPEND_VECTOR_VECTOR_INTEGER, &
		                 APPEND_MATRIX_MATRIX_REAL, APPEND_MATRIX_VECTOR_REAL, APPEND_VECTOR_VECTOR_REAL
	END INTERFACE

  INTERFACE NEAREST_ELEMENT
		MODULE PROCEDURE NEAREST_ELEMENT_REAL, NEAREST_ELEMENT_INTEGER
	END INTERFACE

	INTERFACE SLICE
		MODULE PROCEDURE SLICE_MATRIX_REAL, SLICE_VECTOR_REAL, SLICE_MATRIX_INTEGER, &
		                 SLICE_VECTOR_INTEGER
	END INTERFACE

	INTERFACE SUBMATRIX
		MODULE PROCEDURE SUBMATRIX_SCALAR_REAL, SUBMATRIX_VECTOR_REAL, SUBMATRIX_SCALAR_INTEGER, &
		                 SUBMATRIX_VECTOR_INTEGER, SUBMATRIX_INPLACE_SCALAR_REAL, &
										 SUBMATRIX_INPLACE_VECTOR_REAL, SUBMATRIX_INPLACE_SCALAR_INTEGER, &
								     SUBMATRIX_INPLACE_VECTOR_INTEGER, SUBVECTOR_INPLACE_VECTOR_REAL
	END INTERFACE

  INTERFACE VECTOR_LOC
		MODULE PROCEDURE VECTOR_LOC_REAL, VECTOR_LOC_INTEGER
	END INTERFACE

  INTERFACE READ_SED
    MODULE PROCEDURE READ_SED_SCALAR_, READ_SED_VECTOR_
  END INTERFACE

  INTERFACE APPLY_VD
    MODULE PROCEDURE APPLY_VD_VECTOR_, APPLY_VD_MATRIX_
  END INTERFACE

  INTERFACE FILTER_EXIST
    MODULE PROCEDURE FILTER_EXIST_SCALAR_, FILTER_EXIST_VECTOR_
  END INTERFACE

  INTERFACE GET_FILTER
    MODULE PROCEDURE GET_FILTER_SCALAR_, GET_FILTER_VECTOR_
  END INTERFACE

  INTERFACE EFFECTIVE_FLUX
    MODULE PROCEDURE EFLUX_SCALAR_VECTOR_, EFLUX_VECTOR_VECTOR_, EFLUX_SCALAR_MATRIX_,&
                     EFLUX_VECTOR_MATRIX_
  END INTERFACE

  INTERFACE MEAN
		MODULE PROCEDURE MEAN_VECTOR_, MEAN_MATRIX_, MEAN_MATRIX_REDUX_, WEIGHTED_MEAN_VECTOR_,&
		                 WEIGHTED_MEAN_MATRIX_, WEIGHTED_MEAN_MATRIX_REDUX_
	END INTERFACE

  INTERFACE LINEAR_INTER
    MODULE PROCEDURE LINEAR_INTER_VEC, LINEAR_INTER_MAT
  END INTERFACE

  INTERFACE TRAPZO_INTEG
    MODULE PROCEDURE TRAPZO_INTEG_VEC, TRAPZO_INTEG_MAT
  END INTERFACE

!	BC FILTERS PARAMETERS ============================================================================
	INTEGER, ALLOCATABLE, PUBLIC, PROTECTED             :: ni(:), nf(:), np(:)
	REAL, ALLOCATABLE, PUBLIC, PROTECTED                :: r(:, :)
	CHARACTER(LEN = 64), ALLOCATABLE, PUBLIC, PROTECTED :: fid(:)
	LOGICAL, PUBLIC, PROTECTED                          :: filters_ready = .false.
	LOGICAL, SAVE, PUBLIC                               :: connected(0:99) = .false.
  REAL, PARAMETER, PUBLIC                             :: pi = 4.0 * atan(1.0)
  REAL, PARAMETER, PUBLIC                             :: c_cgs = 2.99792458E10
  REAL, PARAMETER, PUBLIC                             :: h = 6.6261e-27
  REAL, PARAMETER, PUBLIC                             :: Lo = 3.826E33

  CONTAINS

  SUBROUTINE HELP()
    write(6, *)
    write(6, *) "Usage"
    write(6, *) "-----"
    write(6, *) "  dynbas [options] <input_file>"
    write(6, *)
    write(6, *) "  <input_file>         : a plain text file with columns: wlength, flux, sigma."
    write(6, *) "  [options]            : any of the following ones."
    write(6, *)
    write(6, *) "Options"
    write(6, *) "-------"
    write(6, *)
    write(6, *) "  --output-dir=<path>  : path to the directory where log (output files) shall be"
    write(6, *) "                         placed."
    write(6, *) "  --ID=<value>         : ID for the output file as in: dynbasfit_ID_inputfile.log"
    write(6, *) "  --passbands=<values> : a comma-separated list of filter IDs as defined in the"
    write(6, *) "                         BC/CB srcs."
    write(6, *) "  --L-passband=<value> : filter ID for the luminosity-weighted quantities."
    write(6, *) "  --photon-count       : if given and passbands given, the models are computed in"
    write(6, *) "                         photon counts."
    write(6, *) "  --LOSVD-wl=<values>  : the wavelength range where the LOSVD will be computed."
    write(6, *) "                         Defaults to [3900, 4000]."
    write(6, *) "  --LOSVD-grid=<values>: the valid range of LOSVD. If four values given, the last"
    write(6, *) "                         two are taken as the wavelength range where the LOSVD will"
    write(6, *) "                         be computed."
    write(6, *) "  --age-grid=<values>  : age range where the models will be constraint. If a third"
    write(6, *) "                         value is given, it will be set as the step in age indexes."
    write(6, *) "                         e.g., if given --age-grid=0.0,14e9,3, the age grid will be"
    write(6, *) "                         resampled as 1,4,7,...,196."
    write(6, *) "  --ext-curve=<value>  : extinction curve model to use. Valid values are:"
    write(6, *) "                           * CF  : Charlot & Fall 2000."
    write(6, *) "                           * CCM : Cardelli, Clayton & Mathis 1989."
    write(6, *) "                         Defaults to CCM."
    write(6, *) "  --Av-grid=<values>   : dust extinction range (Av) in where the models will be"
    write(6, *) "                         constraint. If a third value is given, it will be taken"
    write(6, *) "                         as the number of values in the grid."
    write(6, *) "  --Zsun=<value>       : assumed value for solar metallicity."
    write(6, *) "  --models=<pattern>   : a set of wilcards (*, ?, [...]) matching .ised files in"
    write(6, *) "                         $DB_MODELS path."
    write(6, *) "  --models!=<pattern>  : likewise, but the pattern will be evaluated for exclusion."
    write(6, *) "  -v                   : verbose mode. Use this option to print on screen"
    write(6, *) "                         meaningful information of the current run and SED fitting"
    write(6, *) "                         setups."
    write(6, *) "  -h, --help           : display this message."
    write(6, *)
    stop
  END SUBROUTINE

  PURE FUNCTION GAUSSIAN(x, mu, sigma)

		implicit none

		real, intent(in) :: x, mu, sigma
		real             :: gaussian

		gaussian = (1. / (sqrt(2 * pi) * sigma)) * exp(- (x - mu) ** 2 / sigma ** 2 / 2)

	END FUNCTION

  SUBROUTINE APPLY_VD_VECTOR_(wlength, flux, vel_disp)

    implicit none

    integer             :: i, j, k, m1, m2, nwl
    integer, parameter  :: m = 6, nx = 100
    real, intent(in)    :: wlength(:)
    real, intent(inout) :: flux(:)
    real, intent(in)    :: vel_disp
    real                :: wlength_(size(wlength)+2*nx), flux_(size(wlength)+2*nx)
    real                :: u(size(wlength)+2*nx), g(size(wlength)+2*nx)
    real                :: xmax
    real                :: c_mks = 3.00e5

    nwl = size(wlength) + 2 * nx

    wlength_(:nx)         = [(wlength(1) - i, i=nx, 1, -1)]
    wlength_(nx+1:nwl-nx) = wlength
    wlength_(nwl-nx+1:)   = [(wlength(nwl-2*nx) + i, i=1, nx)]
    flux_(:nx)            = flux(1)
    flux_(nx+1:nwl-nx)    = flux
    flux_(nwl-nx+1:)      = flux(nwl-2*nx)

    where(flux_ < 0.0) flux_ = 0.0

    do i=1, nwl
      if(flux_(i)==0.0) cycle
      xmax = c_mks*wlength_(i) / (c_mks-m*vel_disp)
      CALL LOCATE(wlength_, nwl, xmax, j)
      m2 = j + 1
      m1 = 2 * i - m2

      if(m1<1)   m1 = 1
      if(m2>nwl) m2 = nwl

      k = 0
      do j=m2, m1, -1
        k = k + 1
        if(flux_(j)==0.0) cycle
        u(k) = (wlength_(i) / wlength_(j)-1.0) * c_mks
        g(k) = flux_(j) * GAUSSIAN(u(k), 0.0, vel_disp)
      end do

      if(i>=nx+1 .and. i<=nwl-nx) flux(i-nx) = TRAPZO_INTEG(u(:k), g(:k))
    enddo

  END SUBROUTINE

  SUBROUTINE APPLY_VD_MATRIX_(wlength, flux, vel_disp)

    implicit none

    integer             :: i, j, k, nwl, nsed, m1, m2
    integer, parameter  :: m = 6, nx = 70
    real, intent(in)    :: wlength(:)
    real, intent(inout) :: flux(:, :)
    real, intent(in)    :: vel_disp
    real                :: wlength_(size(wlength)+2*nx), flux_(size(wlength)+2*nx, size(flux, 2))
    real                :: u(size(wlength)+2*nx), g(size(wlength)+2*nx, size(flux, 2))
    real                :: xmax
    real                :: c_mks = 3.00e5

    nwl = size(wlength) + 2 * nx
    nsed = size(flux, 2)

    wlength_(:nx)         = [(wlength(1)-i, i=nx, 1, -1)]
    wlength_(nx+1:nwl-nx) = wlength
    wlength_(nwl-nx+1:)   = [(wlength(nwl-2*nx)+i, i=1, nx)]
    flux_(:nx, :)         = spread(flux(1, :), 1, nx)
    flux_(nx+1:nwl-nx, :) = flux
    flux_(nwl-nx+1:, :)   = spread(flux(nwl-2*nx, :), 1, nx)

    where(flux_<0.0) flux_ = 0.0

    do i=1, nwl
      if(any(flux_(i, :)==0.0)) cycle
      xmax = c_mks*wlength_(i) / (c_mks-m*vel_disp)
      CALL LOCATE(wlength_, nwl, xmax, j)
      m2 = j + 1
      m1 = 2 * i - m2

      if (m1<1)   m1 = 1
      if (m2>nwl) m2 = nwl

      k = 0
      do j=m2, m1, -1
        k = k + 1
        if(any(flux_(j, :)==0.0)) cycle
        u(k)    = (wlength_(i) / wlength_(j)-1.0) * c_mks
        g(k, :) = flux_(j, :) * GAUSSIAN(u(k), 0.0, vel_disp)
      end do

      if(i>=nx+1 .and. i<=nwl-nx) flux(i-nx, :) = TRAPZO_INTEG(u(:k), g(:k, :))
    end do

  END SUBROUTINE

  ELEMENTAL FUNCTION EXT(id, wlength, Av, Rv)

    implicit none

    character(*), intent(in)   :: id
    real, intent(in)           :: wlength, Av
    real, intent(in), optional :: Rv
    real                       :: Rv_, EXT

    Rv_ = 3.1
    if(present(Rv)) Rv_ = Rv

    select case(id)
    case("CF", "cf")
      EXT = CF(wlength, Av)
    case default
      EXT = CCM(wlength, Av, Rv_)
    end select

    CONTAINS
    ELEMENTAL FUNCTION CF(wlength, Av)

      implicit none

      real, intent(in) :: wlength, Av
      real             :: CF

      CF = exp(-Av/1.086 * (wlength/5500.0)**(-0.7))

    END FUNCTION

    ELEMENTAL FUNCTION CCM(wlength, Av, Rv)

      implicit none

      real, intent(in)           :: wlength, Av
      real, intent(in), optional :: Rv
      real                       :: x, y, a, b, fa, fb, CCM

      x = 1.e4 / wlength

  !	Infrared:  Eq. 2
      if (x >= 0.3 .and. x < 1.1) then
        a= 0.574*x**1.61
        b= -0.527*x**1.61
  !	Optical/NIR:  Eq. 3
      else if (x >= 1.1 .and. x < 3.3) then
        y=x-1.82
        a=1.+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7
        b=1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
  !	Ultraviolet:  Eq.4
      else if (x >= 3.3 .and. x < 8.) then
        fa=0.
        fb=0.
          if(x.ge.5.9) then
            fa=-0.04473*(x-5.9)**2-0.009779*(x-5.9)**3
            fb=0.2130*(x-5.9)**2+0.1207*(x-5.9)**3
          end if
        a=1.752-0.316*x-0.104/((x-4.67)**2+0.341)+fa
        b=-3.090+1.825*x+1.206/((x-4.62)**2+0.263)+fb
  !	Far-UV:  Eq. 5
      else if (x.ge.8.and.x.le.10.) then
        a=-1.073-0.628*(x-8.)+0.137*(x-8.)**2-0.070*(x-8)**3
        b=13.670+4.257*(x-8)-0.420*(x-8)**2+0.374*(x-8)**3
      else
        CCM = 0.
        return
      end if

      if(present(Rv)) then
        CCM = 10 ** (-0.4 * Av * (a + b / Rv))
      else
        CCM = 10 ** (-0.4 * Av * (a + b / 3.1))
      end if

    END FUNCTION

  END FUNCTION

  FUNCTION EFLUX_SCALAR_VECTOR_(filter_ids, wlength, flux, z, photon_count) RESULT(EFLUX)

		implicit none

		integer                       :: n, pos(2)
		integer, intent(in)           :: filter_ids
		real, intent(in)              :: wlength(:), flux(:), z
		real, allocatable             :: wl_res(:), res(:)
		real, allocatable             :: flux_(:)
		real                          :: resolution, wlr_min, wlr_max, EFLUX
    logical, intent(in), optional :: photon_count

    EFLUX = 0.0

		CALL GET_FILTER(filter_ids, wl_res, res, echo=.false.)

    n = size(wl_res)
    wlr_min = wl_res(1)/(z+1)
    wlr_max = wl_res(n)/(z+1)
    pos(1) = maxloc(wlength, 1, wlength<wlr_min)
    pos(2) = minloc(wlength, 1, wlength>wlr_max)

    if(any(pos<=0)) return

    resolution = MEAN(wlength(pos(1)+1:pos(2))-wlength(pos(1):pos(2)-1))
    if(pos(2)-pos(1)+1<0.95/resolution*(wlr_max-wlr_min)) return
    if(any(flux(pos(1):pos(2))<=0.0)) return

    flux_ = LINEAR_INTER(wlength(pos(1):pos(2)), flux(pos(1):pos(2)), wl_res/(z+1))

    if(present(photon_count) .and. photon_count) then
      EFLUX = MEAN(wl_res*flux_, weights=res)/(h*c_cgs)*Lo
    else
      EFLUX = MEAN(flux_, weights=res)
    end if

	END FUNCTION

  FUNCTION EFLUX_SCALAR_MATRIX_(filter_ids, wlength, flux, z, photon_count) RESULT(EFLUX)

		implicit none

		integer                       :: n, pos(2)
		integer, intent(in)           :: filter_ids
		real, intent(in)              :: wlength(:), flux(:, :), z
		real, allocatable             :: wl_res(:), res(:)
		real, allocatable             :: flux_(:, :)
		real                          :: resolution, wlr_min, wlr_max, EFLUX(size(flux, 2))
    logical, allocatable          :: masked(:)
    logical, intent(in), optional :: photon_count

    EFLUX = 0.0

		CALL GET_FILTER(filter_ids, wl_res, res, echo=.false.)

    n = size(wl_res)
    wlr_min = wl_res(1)/(z+1)
    wlr_max = wl_res(n)/(z+1)
    pos(1) = maxloc(wlength, 1, wlength<wlr_min)
    pos(2) = minloc(wlength, 1, wlength>wlr_max)

    if(any(pos<=0)) return

    resolution = MEAN(wlength(pos(1)+1:pos(2))-wlength(pos(1):pos(2)-1))
    if(pos(2)-pos(1)+1<0.95/resolution*(wlr_max-wlr_min)) return

    masked = any(flux(pos(1):pos(2), :)<=0.0, dim=1)
    flux_  = LINEAR_INTER(wlength(pos(1):pos(2)), flux(pos(1):pos(2), :), wl_res/(z+1))

    if(present(photon_count) .and. photon_count) then
      where(.not.masked)
        EFLUX = MEAN(spread(wl_res, 2, size(flux, 2))*flux_, weights=spread(res, 2, size(flux, 2)), dim=1)/(h*c_cgs)*Lo
      end where
    else
      where(.not.masked)
        EFLUX = MEAN(flux_, weights=spread(res, 2, size(flux, 2)), dim=1)
      end where
    end if

	END FUNCTION

	FUNCTION EFLUX_VECTOR_VECTOR_(filter_ids, wlength, flux, z, photon_count) RESULT(EFLUX)

		implicit none

		integer, intent(in)           :: filter_ids(:)
		integer                       :: i(size(filter_ids), 2), j, k, l, pos(2)
		real, intent(in)              :: wlength(:), flux(:), z
		real, allocatable             :: wl_res(:), res(:)
		real, allocatable             :: flux_(:)
		real                          :: resolution, wlr_min, wlr_max
    real                          :: EFLUX(size(filter_ids))
    logical, intent(in), optional :: photon_count

    EFLUX = 0.0

		CALL GET_FILTER(filter_ids, wl_res, res, i, echo=.false.)

    do j=1, size(filter_ids)
      k = i(j, 1)
      l = i(j, 2)
      wlr_min = minval(wl_res(k:l))/(z+1)
      wlr_max = maxval(wl_res(k:l))/(z+1)

      pos = [maxloc(wlength, 1, wlength<wlr_min), minloc(wlength, 1, wlength>wlr_max)]

      if(any(pos<=0)) cycle

      resolution = MEAN(wlength(pos(1)+1:pos(2))-wlength(pos(1):pos(2)-1))
      if(pos(2)-pos(1)+1<0.95/resolution*(wlr_max-wlr_min)) cycle
      if(any(flux(pos(1):pos(2))<=0.0)) cycle

      flux_ = LINEAR_INTER(wlength(pos(1):pos(2)), flux(pos(1):pos(2)), wl_res(k:l)/(z+1))

      if(present(photon_count) .and. photon_count) then
        EFLUX(j) = MEAN(wl_res(k:l)*flux_, weights=res(k:l))/(h*c_cgs)*Lo
      else
        EFLUX(j) = MEAN(flux_, weights=res(k:l))
      end if
    end do

	END FUNCTION

	FUNCTION EFLUX_VECTOR_MATRIX_(filter_ids, wlength, flux, z, photon_count) RESULT(EFLUX)

		implicit none

		integer, intent(in)           :: filter_ids(:)
		integer                       :: i(size(filter_ids), 2), j, k, l, pos(2)
		real, intent(in)              :: wlength(:), flux(:, :), z
		real, allocatable             :: wl_res(:), res(:)
		real, allocatable             :: res_(:, :), flux_(:, :)
		real                          :: resolution, wlr_min, wlr_max
    real                          :: EFLUX(size(filter_ids), size(flux, 2))
    logical, allocatable          :: masked(:)
    logical, intent(in), optional :: photon_count

    EFLUX = 0.0

		CALL GET_FILTER(filter_ids, wl_res, res, i, echo=.false.)
    res_ = spread(res, 2, size(flux, 2))

    do j=1, size(filter_ids)
      k = i(j, 1)
      l = i(j, 2)
      wlr_min = minval(wl_res(k:l))/(z+1)
      wlr_max = maxval(wl_res(k:l))/(z+1)

      pos = [maxloc(wlength, 1, wlength<wlr_min), minloc(wlength, 1, wlength>wlr_max)]

      if(any(pos<=0)) cycle

      resolution = MEAN(wlength(pos(1)+1:pos(2))-wlength(pos(1):pos(2)-1))
      if(pos(2)-pos(1)+1<0.95/resolution*(wlr_max-wlr_min)) cycle

      masked = any(flux(pos(1):pos(2), :)<=0.0, dim=1)
      flux_  = LINEAR_INTER(wlength(pos(1):pos(2)), flux(pos(1):pos(2), :), wl_res(k:l)/(z+1))

      if(present(photon_count) .and. photon_count) then
        where(.not.masked)
          EFLUX(j, :) = MEAN(spread(wl_res, 2, size(flux, 2))*flux_, weights=res_(k:l, :), dim=1)/(h*c_cgs)*Lo
        end where
      else
        where(.not.masked)
          EFLUX(j, :) = MEAN(flux_, weights=res_(k:l, :), dim=1)
        end where
      end if
    end do

	END FUNCTION

  FUNCTION ERROR_PROPAGATION(filter_id, wlength, sigma_flux, photon_count) RESULT(ERROR)

		implicit none

		integer                       :: n, pos(2)
		integer, intent(in)           :: filter_id
		real, intent(in)              :: wlength(:), sigma_flux(:)
		real, allocatable             :: wl_res(:), res(:), weights(:), res_(:)
		real                          :: resolution, ERROR
    logical, intent(in), optional :: photon_count

    ERROR = 0.0

		CALL GET_FILTER(filter_id, wl_res, res, echo = .false.)

    n = size(wl_res)
    pos(1) = minloc(wlength, 1, wlength > wl_res(1))
    pos(2) = maxloc(wlength, 1, wlength < wl_res(n))

    if(any(pos == 0)) return

    resolution = MEAN(wlength(pos(1)+1:pos(2))-wlength(pos(1):pos(2)-1))

    if(pos(2)-pos(1)+1 < 0.95/resolution*(wl_res(n) - wl_res(1))) return
    if(any(sigma_flux(pos(1):pos(2)) <= 0.0)) return

		allocate(res_(pos(2)-pos(1)+1), weights(pos(2)-pos(1)+1))

		res_    = LINEAR_INTER(wl_res, res, wlength(pos(1):pos(2)))
		weights = res_ / TRAPZO_INTEG(wlength(pos(1):pos(2)), res_)

    if(present(photon_count) .and. photon_count) then
      ERROR = sqrt(TRAPZO_INTEG(wlength(pos(1):pos(2)), (weights*wlength(pos(1):pos(2))*sigma_flux(pos(1):pos(2)))**2))/(h*c_cgs)*Lo
    else
      ERROR = sqrt(TRAPZO_INTEG(wlength(pos(1):pos(2)), (weights*sigma_flux(pos(1):pos(2)))**2))
    end if

  END FUNCTION

  SUBROUTINE FILTER_INFO(filter_id, area, eff_wl, mean_wl, fwhm, width, filter_name)

		implicit none

		integer                                    :: i, loc1(1), loc2(1), pos(2)
		integer                                    :: unit, nhr, nwl
		integer, intent(in)                        :: filter_id
		real, intent(out), optional                :: area, eff_wl, mean_wl, fwhm, width
		real, dimension(:), allocatable            :: wl_res, res, vega_wl, vega_flux
		real                                       :: res_max, wl_i, wl_f
		real                                       :: area_
		character(len = 64), intent(out), optional :: filter_name
		character(len = :), allocatable            :: dir_inp
		save                                       :: vega_wl, vega_flux

    allocate(dir_inp, source="../data/")

		CALL GET_FILTER(filter_id, wl_res, res, echo = .false.)

		res_max = 0.0
		unit    = AVAILABLE_UNIT()

		if(present(area).or.present(mean_wl).or.present(width)) then
			area_   = TRAPZO_INTEG(wl_res, res)
			res_max = maxval(res)

			if(present(area))    area    = area_
			if(present(mean_wl)) mean_wl = TRAPZO_INTEG(wl_res, wl_res * res) / area_
			if(present(width))   width   = area_ / maxval(res)
		end if
		if(present(eff_wl)) then
			if(.not.allocated(vega_flux)) then
				CALL INQUIRE_FILE(dir_inp//"A0V_KURUCZ_92.SED", nheader = nhr, nrow = nwl)

				open(unit, file = dir_inp//"A0V_KURUCZ_92.SED", action = "read")
				allocate(vega_wl(nwl), vega_flux(nwl))
				do i = 1, nhr
					read(unit, *)
				end do
				do i = 1, nwl
					read(unit, *) vega_wl(i), vega_flux(i)
				end do
				close(unit)
			end if

      pos(1) = minloc(vega_wl - wl_res(1), 1, vega_wl > wl_res(1))
      pos(2) = maxloc(vega_wl - wl_res(size(wl_res)), 1, vega_wl < wl_res(size(wl_res)))

      vega_wl   = vega_wl(pos(1):pos(2))
			vega_flux = vega_flux(pos(1):pos(2)) * LINEAR_INTER(wl_res, res, vega_wl)
			eff_wl    = TRAPZO_INTEG(vega_wl, vega_wl * vega_flux) / TRAPZO_INTEG(vega_wl, vega_flux)
		end if
		if(present(fwhm)) then
			if(res_max == 0.0) res_max = maxval(res)
			loc1 = maxloc(res)

			loc2 = minval(abs(res(:loc1(1)) - res_max / 2))
			wl_i = wl_res(loc2(1))
			loc2 = minval(abs(res(loc1(1):) - res_max / 2))
			wl_f = wl_res(loc2(1))

			fwhm = wl_f - wl_i
		end if
		if(present(filter_name)) filter_name = fid(filter_id)

	END SUBROUTINE

  FUNCTION FILTER_EXIST_SCALAR_(filter_ids) RESULT(EXIST)

    implicit none

    integer, intent(in) :: filter_ids
    logical             :: EXIST

    if(.not.filters_ready) CALL READ_BC_FILTERS_(echo=.false.)

    EXIST = filter_ids<=size(np)

  END FUNCTION

  FUNCTION FILTER_EXIST_VECTOR_(filter_ids) RESULT(EXIST)

    implicit none

    integer, intent(in) :: filter_ids(:)
    logical             :: EXIST(size(filter_ids))

    if(.not.filters_ready) CALL READ_BC_FILTERS_(echo=.false.)

    EXIST = filter_ids<=size(np)

  END FUNCTION

	SUBROUTINE GET_FILTER_SCALAR_(filter_ids, wlength, response, filter_name, echo)

		implicit none

		integer, intent(in)                  :: filter_ids
    real, allocatable, intent(out)       :: wlength(:), response(:)
    character(64), intent(out), optional :: filter_name
    logical, intent(in)                  :: echo

    if(.not. FILTER_EXIST(filter_ids)) then
      write(6, *) "STOP GET_FILTER: filter "//STRING(filter_ids)//" is not defined in $DB_FILTERS."
      stop
    end if

		allocate(wlength(np(filter_ids)), response(np(filter_ids)))
		wlength  = r(ni(filter_ids) : nf(filter_ids), 1)
		response = r(ni(filter_ids) : nf(filter_ids), 2)

    if(present(filter_name)) filter_name = fid(filter_ids)

    if(echo) write(6, *) np(filter_ids), "data points in ", trim(fid(filter_ids))

	END SUBROUTINE

  SUBROUTINE GET_FILTER_VECTOR_(filter_ids, wlength, response, positions, filter_name, echo)

		implicit none

    integer                                           :: i, n
		integer, intent(in)                               :: filter_ids(:)
    integer, intent(out)                              :: positions(size(filter_ids), 2)
    real, allocatable, intent(out)                    :: wlength(:), response(:)
    character(64), allocatable, intent(out), optional :: filter_name(:)
    logical, intent(in)                               :: echo

    if(.not. all(FILTER_EXIST(filter_ids))) then
      write(6, *) "STOP GET_FILTER: some filter(s) may not be defined in $DB_FILTERS."
      stop
    end if

    n = sum([(np(filter_ids(i)), i=1, size(filter_ids))])
		allocate(wlength(n), response(n))
    wlength(1:np(filter_ids(1)))  = r(ni(filter_ids(1)):nf(filter_ids(1)), 1)
  	response(1:np(filter_ids(1))) = r(ni(filter_ids(1)):nf(filter_ids(1)), 2)
    positions(1, :) = [1, np(filter_ids(1))]
    do i=2, size(filter_ids)
      wlength(sum(np(filter_ids(1:i-1)))+1:sum(np(filter_ids(1:i))))  = r(ni(filter_ids(i)):nf(filter_ids(i)), 1)
  		response(sum(np(filter_ids(1:i-1)))+1:sum(np(filter_ids(1:i)))) = r(ni(filter_ids(i)):nf(filter_ids(i)), 2)
      positions(i, :) = [sum(np(filter_ids(1:i-1)))+1, sum(np(filter_ids(1:i)))]
    end do

    if(present(filter_name)) filter_name = [(fid(filter_ids(i)), i=1, size(filter_ids))]

    if(echo) then
      do i=1, size(filter_ids)
        write(6, *) np(filter_ids(i)), "data points in ", trim(fid(filter_ids(i)))
      end do
    end if

	END SUBROUTINE

  SUBROUTINE READ_BC_FILTERS_(echo)

    implicit none

    integer              :: i, nfil, ntot, status
    integer              :: unit
    character(len = 256) :: filters_bin
    logical, intent(in)  :: echo

		CALL GETENV("DB_FILTERS", filters_bin)

    unit = AVAILABLE_UNIT()
		open(unit, file=trim(filters_bin), form="unformatted", action="read", iostat=status)
    CALL IOERR(status, trim(filters_bin))
		read(unit) nfil

		allocate(ni(nfil), nf(nfil), np(nfil), fid(0 : nfil))

		rewind(unit)
		read(unit) nfil, (ni(i), i = 1, nfil), (nf(i), i = 1, nfil), (np(i), i = 1, nfil)

		ntot = sum(np)
		allocate(r(ntot, 2))

		rewind(unit)
		read(unit) nfil, (ni(i), i = 1, nfil), (nf(i), i = 1, nfil), (np(i), i = 1, nfil),&
               (fid(i), i = 0, nfil), ntot, (r(i, 1), i = 1, ntot), (r(i, 2), i = 1, ntot)
		close(unit)

		filters_ready = .true.

		if(echo) then
			i = scan(filters_bin, "/", back = .true.)
			write(6, *) "BC filters data from '", trim(filters_bin(i:)), "':"
			write(6, *) nfil,           "filters found"
			write(6, *) ntot,           "total data points"
			write(6, *) sum(np) / nfil, "mean data points per filter"
			write(6, *)
		end if

  END SUBROUTINE

  SUBROUTINE READ_SED_SCALAR_(file_sed, wlength, flux, ages, metallicity, ised, echo)

    implicit none

    integer                                  :: i, j, nwl, nsed, nhdr, ider, unit
    real, allocatable, intent(out)           :: wlength(:), flux(:, :)
    real, allocatable, intent(out), optional :: ages(:)
    real, intent(out), optional              :: metallicity
    real, allocatable                        :: ages_(:)
    character(len = *), intent(in)           :: file_sed
    logical, intent(in)                      :: ised, echo

		integer                                  :: ml, mu, iseg, jo
    integer, parameter                       :: imf = 10
    real                                     :: xx(imf), lm(imf), um(imf), baux(imf), cc(imf)
    real                                     :: cn(imf), totm, totn, avs, tau
    character(len = 80)                      :: id

		unit = AVAILABLE_UNIT()

    if(ised) then

      open(unit, file = trim(file_sed), form = "unformatted", status = "old", iostat = ider)

      CALL IOERR(ider, file_sed)

      read(unit) nsed
      read(unit) nwl

      rewind(unit)
      allocate(ages_(nsed), wlength(nwl), flux(nwl, nsed))

      read(unit) nsed, (ages_(i), i = 1, nsed), ml, mu, iseg, &
                 (xx(i), lm(i), um(i), baux(i), cn(i), cc(i), i = 1, iseg), totm, totn, avs, jo, tau, id
      read(unit) nwl, (wlength(i), i = 1, nwl)

      do j = 1, nsed
        read(unit) nwl, (flux(i, j), i = 1, nwl)
      end do
      close(unit)

			if(present(ages)) then
				allocate(ages(nsed))
				ages = ages_
			end if
			if(present(metallicity)) read(id(index(id, "Z=", back = .true.) + 2:), *) metallicity

    else

      CALL INQUIRE_FILE(file_sed, nrow = nwl, ncolumn = nsed, nheader = nhdr)

			nsed = nsed - 1

      allocate(wlength(nwl), flux(nwl, nsed))

      open(unit, file = trim(file_sed), status = "old", action = "read")

      do i = 1, nhdr
        read(unit, *)
      end do

      do i = 1, nwl
        read(unit, *) wlength(i), (flux(i, j), j = 1, nsed)
      end do
      close(unit)
    end if

    if(echo) then
      i = index(file_sed, "/", back = .true.)
      write(6, "(/X,3A)") "S.E.D. file '", trim(file_sed(i + 1:)), "' contains:"
      write(6, "(X,I6,A)") nsed, " spectra."
      write(6, "(X,I6,A/)") nwl, " wavelength points."
    end if

  END SUBROUTINE

  SUBROUTINE READ_SED_VECTOR_(file_sed, wlength, flux, ages, metallicity, ised, echo)

    implicit none

    integer                                  :: i, j, k, nfiles, nwl, nsed, nhdr, ider, unit
    real, allocatable, intent(out)           :: wlength(:), flux(:, :, :)
    real, allocatable, intent(out), optional :: ages(:), metallicity(:)
    real, allocatable                        :: ages_(:), metallicity_(:)
    character(len = *), intent(in)           :: file_sed(:)
    logical, intent(in)                      :: ised, echo

		integer                                  :: ml, mu, iseg, jo
    integer, parameter                       :: imf = 10
    real                                     :: xx(imf), lm(imf), um(imf), baux(imf), cc(imf)
    real                                     :: cn(imf), totm, totn, avs, tau
    character(len = 80)                      :: id

    unit   = AVAILABLE_UNIT()
    nfiles = size(file_sed)

    if(ised) then

      do k = 1, nfiles

        open(unit, file = trim(file_sed(k)), form = "unformatted", status = "old", iostat = ider)

        CALL IOERR(ider, file_sed(k))

        if(k == 1) then
          read(unit) nsed
          read(unit) nwl

          rewind(unit)

          allocate(ages_(nsed), wlength(nwl), flux(nwl, nsed, nfiles), metallicity_(nfiles))

          read(unit) nsed, (ages_(i), i = 1, nsed), ml, mu, iseg, &
                    (xx(i), lm(i), um(i), baux(i), cn(i), cc(i), i = 1, iseg), totm, totn, avs, jo,&
                    tau, id
          read(unit) nwl, (wlength(i), i = 1, nwl)

        else

          read(unit) nsed, (ages_(i), i = 1, nsed), ml, mu, iseg, &
                    (xx(i), lm(i), um(i), baux(i), cn(i), cc(i), i = 1, iseg), totm, totn, avs, jo,&
                    tau, id
          read(unit)

        end if

        read(id(index(id, "Z=", back = .true.) + 2:), *) metallicity_(k)

        do j = 1, nsed
          read(unit) nwl, (flux(i, j, k), i = 1, nwl)
        end do

        close(unit)

      end do

      if(present(ages)) then
				allocate(ages(nsed))
				ages = ages_
			end if
			if(present(metallicity)) then
				allocate(metallicity(nfiles))
				metallicity = metallicity_
			end if

    else

      CALL INQUIRE_FILE(file_sed(1), nrow = nwl, ncolumn = nsed, nheader = nhdr)

      allocate(wlength(nwl), flux(nwl, nsed, nfiles))

      do k = 1, nfiles

        open(unit, file = trim(file_sed(k)), status = "old", action = "read")

        do i = 1, nhdr
          read(unit, *)
        end do

        do i = 1, nwl
          read(unit, *) wlength(i), (flux(i, j, k), j = 1, nsed)
        end do

        close(unit)

      end do

    end if

    if(echo) then
      write(6, "(/X,I2,A)") nfiles, " SED files containing:"
      write(6, "(X,I6,A)") nsed, " spectra."
      write(6, "(X,I6,A/)") nwl," wavelength points."
    end if

  END SUBROUTINE

  SUBROUTINE MONOCHROMATIC_SEDS(wlength, fluxes, filter_id, z, selection, mono_fluxes)

    implicit none

    integer                        :: nsed
    integer, intent(in)            :: filter_id
    integer, intent(in), optional  :: selection(:)
    real, intent(in)               :: z, wlength(:), fluxes(:, :)
    real, allocatable              :: sel_fssps(:, :)
    real, allocatable, intent(out) :: mono_fluxes(:)

    if(present(selection)) then
      CALL SUBMATRIX(fluxes, 2, selection, sel_fssps)
      nsed = size(selection)
    else
      nsed = size(fluxes, 2)
      allocate(sel_fssps(size(wlength), nsed))
      sel_fssps = fluxes
    end if

    allocate(mono_fluxes(nsed))
		mono_fluxes = EFFECTIVE_FLUX(filter_id, wlength, sel_fssps, z)

  END SUBROUTINE

  FUNCTION AVAILABLE_UNIT() RESULT(unit)
! This routine retrieve a disconnected unit for use.
!
! OUTPUT PARAMETERS:
! -----------------
!   * unit     : integer value with the resulting available unit.

    integer :: unit, last_unit = 0

		if(all(connected)) STOP "AVAILABLE UNIT: no available unit found."
    do while(connected(last_unit))
      last_unit = last_unit + 1
      if(last_unit == 100) last_unit = 0
    end do
    unit = last_unit

  END FUNCTION

  SUBROUTINE IOERR(status, filename)
! This routine checks for the status of some unit by verifying the value of iostat.
!
! INPUT PARAMETERS:
! ----------------
!   * status   : integer number with the iostat value. If different to zero an error message is
!                displayed on screen and the program is stopped.
!   * filename : name of the file to check if an error ocurred during an IO transaction.

    implicit none

    integer, intent(in)            :: status
    integer, parameter             :: nstat = 46
    integer                        :: i, iostat_values(nstat)
    character(len = *), intent(in) :: filename
    character(len = 256)           :: messages(nstat), msg(1)
    logical                        :: messages_ready = .false.

    save                           :: iostat_values, messages, messages_ready

    if(status == 0) return

    if(.not. messages_ready) then

      iostat_values(1) = - 1
      messages(1)      = "end of file reached in the file 'X'."

      iostat_values(2) = - 2
      messages(2)      = "end of record reached in the file 'X'."

      iostat_values(3) = 1
      messages(3)      = "the record specified in the file 'X' doesn't exist."

      iostat_values(5) = 2
      messages(5)      = "the file 'X' doesn't exist."

      iostat_values(4) = 6
      messages(4)      = "end of file reached in the file 'X'."

      iostat_values(6) = 10
      messages(6)      = "error reading on direct access the file 'X'."

      iostat_values(7) = 11
      messages(7)      = "error writing on direct access the file 'X'."

      iostat_values(8) = 12
      messages(8)      = "error reading on sequential the file 'X'."

      iostat_values(9) = 13
      messages(9)      = "error writing on sequential the file 'X'."

      iostat_values(10) = 14
      messages(10)      = "error opening the file 'X'."

      iostat_values(11) = 15
      messages(11)      = "permanent IO error found in the file 'X'."

      iostat_values(12) = 37
      messages(12)      = "out of memory when attempting IO operation in the file 'X'."

      iostat_values(13) = 38
      messages(13)       = "error rewinding file 'X'."

      iostat_values(14) = 39
      messages(14)      = "error on endfile operation in the file 'X'."

      iostat_values(15) = 40
      messages(15)      = "error backspacing the file 'X'."

      iostat_values(16) = 107
      messages(16)      = "opening an already existent file, 'X', with status='NEW'."

      iostat_values(17) = 119
      messages(17)      = "error backspacing the file 'X'."

      iostat_values(18) = 122
      messages(18)      = "incomplete record found on direct access the file 'X'."

      iostat_values(19) = 130
      messages(19)      = "action='READWRITE' specified on the OPEN statement to connect the pipe &
                           'X'."

      iostat_values(20) = 135
      messages(20)      = "this program is making calls to an unsupported version of the 'X'L &
                           Fortran run-time environment."

      iostat_values(21) = 139
      messages(21)      = "IO operation not allowed for the file 'X'. Check the specifier action on &
                           the OPEN statement."

      iostat_values(22) = 142
      messages(22)      = "close error in the file 'X'."

      iostat_values(23) = 144
      messages(23)      = "inquire error in the file 'X'."

      iostat_values(24) = 152
      messages(24)      = "access='DIRECT' specified for the file 'X', which can only be accessed &
                           sequentially."

      iostat_values(25) = 153
      messages(25)      = "invalid value of position in openning the pipe 'X'."

      iostat_values(26) = 156
      messages(26)      = "invalid value of recl in opening the file 'X'."

      iostat_values(27) = 159
      messages(27)      = "couldn't flush external the file 'X' input."

      iostat_values(28) = 165
      messages(28)      = "the record number of the next record in the file 'X' is out of range of &
                           the value specified for nextrec."

      iostat_values(29) = 169
      messages(29)      = "attempting asynchronous IO operation on the file 'X' for synchronous IO &
                           only."

      iostat_values(30) = 172
      messages(30)      = "couldn't open the file 'X' in asynchronous mode."

      iostat_values(31) = 173
      messages(31)      = "attempting asynchronous READ/WRITE in the file 'X' when another &
                           WRITE/READ is pending."

      iostat_values(32) = 174
      messages(32)      = "attempting a synchronous IO operation in the file 'X' when an &
                           asynchronous is pending."

      iostat_values(33) = 175
      messages(33)      = "invalid value of id in WAIT statement for the file 'X'."

      iostat_values(34) = 176
      messages(34)      = "WAIT statement couldn't be completed in the file 'X' because the &
                           corresponding asynchronous IO operation is in a different scoping unit."

      iostat_values(35) = 178
      messages(35)      = "attempting two asynchronous WRITE in the file 'X' in the same record."

      iostat_values(36) = 179
      messages(36)      = "attempting IO operation in the file 'X' when another IO operation is &
                           pending."

      iostat_values(37) = 181
      messages(37)      = "attempting to connect the file 'X' to multiple units in asynchronous &
                           mode."

      iostat_values(38) = 182
      messages(38)      = "invalid value of uwidth in the file 'X'. It must be set to either 32 or &
                           64."

      iostat_values(39) = 183
      messages(39)      = "the maximun record length for the file 'X' is out of the range of the &
                           value of recl in the INQUIRE statement."

      iostat_values(40) = 184
      messages(40)      = "attempting to transmit more bytes than the specified in the value of &
                           size/num in the IO statement on the file 'X'."

      iostat_values(41) = 185
      messages(41)      = "attempting to open the file 'X' to two units with different uwidth &
                           values."

      iostat_values(42) = 186
      messages(42)      = "error opening the file 'X'. The value of unit must be between 0 and &
                           2,147,483,647."

      iostat_values(43) = 192
      messages(43)      = "the value of the position in the file 'X' is out of the range of the &
                           value specified for pos in the INQUIRE statement."

      iostat_values(44) = 193
      messages(44)      = "the value of the size in the file 'X' is out of the range of the value &
                           specified for size in the INQUIRE statement."

      iostat_values(45) = 200
      messages(45)      = "error flushing the file 'X'."

      iostat_values(46) = 201
      messages(46)      = "attempting to flush the non-seekable file 'X'."
    end if

    if(.not. any(iostat_values == status)) then
      write(6, *) "STOP IOERR: unknown error on IO operation in file ", trim(filename), "."
      stop
    end if

    msg = pack(messages, iostat_values == status)

    i = index(msg(1), "X")
    write(6, "('STOP IOERR: ', A)") msg(1)(:i - 1)//trim(filename)//msg(1)(i + 1:)
    stop

  END SUBROUTINE

  SUBROUTINE INQUIRE_FILE(filename, nheader, nrow, ncolumn, status)
! This routine inquire by name the arrangement of a file under the following assumptions:
!    1) The file is an ASCII table, possibly containing one or more lines of header starting with #
!       and columns separated by spaces (" ").
!    2) There are no missing values footer neither.
!    3) The maximum number of columns is 1000. This can be changed via the parameter maxcol.
!
! INPUT PARAMETERS:
! ----------------
!   * filename : character string with the name of the file which is intented to be inquired.
!   * nheader  : integer number. If present, the number of lines corresponding to the header.
!   * nrow     : integer number. If present, the number of lines in the file without the header.
!   * ncolumn  : integer number. If present, the number of columns in the file.
!   * status   : integer number. If present, the value of the iostat parameter returned by open. If
!                different to zero, the file doesn't exist and the requested values are set to zero.

    implicit none

    integer                        :: in_nheader, unit, in_status, of_stat
    integer                        :: maxcol
    integer, intent(out), optional :: nrow, ncolumn, nheader, status
    real                           :: smallest = tiny(smallest)
    real, allocatable              :: row(:)
    character(len = *), intent(in) :: filename
    character(len = 70000)         :: reg
    logical                        :: present_status, present_nheader, present_ncolumn, present_nrow
    parameter(maxcol = 1000)

    reg             = ""
    present_status  = present(status)
    present_nheader = present(nheader)
    present_ncolumn = present(ncolumn)
    present_nrow    = present(nrow)

    unit = AVAILABLE_UNIT()
    open(unit, file = filename, action = "read", status = "old", iostat = in_status)
    if(present_status) status = in_status
    if(.not. present_nrow .and. .not. present_ncolumn .and. .not. present_nheader) then
      close(unit)
      return
    end if

    select case(in_status)
    case(0)
      allocate(row(maxcol))
      row = smallest

      in_nheader = 0
      do
        read(unit, "(A)", iostat = of_stat) reg

        if(index(adjustl(reg), "#") == 1) then
          in_nheader = in_nheader + 1
        else
          if(present_nheader) nheader = in_nheader
          if(.not. present_ncolumn .and. .not. present_nrow) then
	          close(unit)
	          return
          end if
          if(.not. present_ncolumn) exit

          ncolumn = size(SPLIT(reg, " "))
          if(.not. present_nrow) then
						close(unit)
	          return
	        end if
          exit
        end if

        if(of_stat == - 1) exit
      end do

      rewind(unit)
      of_stat = 0
      nrow    = - in_nheader - 1
      do while(of_stat /= - 1)
        read(unit, *, iostat = of_stat)
        nrow = nrow + 1
      end do
			close(unit)
    case default
      if(present_nheader) nheader = 0
      if(present_ncolumn) ncolumn = 0
      if(present_nrow)    nrow    = 0
      close(unit)
    end select

  END SUBROUTINE

  FUNCTION LINEAR_INTER_VEC(x, y, x_new) RESULT(y_new)

    implicit none

    integer          :: i, n
    real, intent(in) :: x(:), y(:), x_new(:)
    real             :: y_new(size(x_new)), m(size(x)-1), b(size(x)-1)

    n = size(x)
		m = (y(2:) - y(:n-1)) / (x(2:) - x(:n-1))
		b = y(:n-1) - m * x(:n-1)
    do i=1, n-1
      where(x(i) <= x_new .and. x_new <= x(i+1)) y_new = m(i) * x_new + b(i)
    end do

  END FUNCTION

  FUNCTION LINEAR_INTER_MAT(x, y, x_new) RESULT(y_new)

    implicit none

    integer          :: i, j, n, l
    real, intent(in) :: x(:), y(:, :), x_new(:)
    real             :: m(size(x)-1, size(y, 2)), b(size(x)-1, size(y, 2))
    real             :: y_new(size(x_new), size(y, 2))

    n = size(x)
    l = size(y, 2)

    m = (y(2:, :) - y(:n-1, :)) / spread(x(2:) - x(:n-1), 2, size(y, 2))
    b = y(:n-1, :) - m * spread(x(:n-1), 2, size(y, 2))
    do j=1, l
      do i=1, n-1
        where(x(i) <= x_new .and. x_new <= x(i+1)) y_new(:, j) = m(i, j) * x_new + b(i, j)
      end do
    end do

  END FUNCTION

	FUNCTION TRAPZO_INTEG_VEC(x, y)

		implicit none

		integer          :: n
		real, intent(in) :: x(:), y(:)
		real             :: trapzo_integ_vec

		n = size(x)
		trapzo_integ_vec = sum(abs(x(2:) - x(:n-1)) * (y(2:) + y(:n-1))) * 0.5

	END FUNCTION

  FUNCTION TRAPZO_INTEG_MAT(x, y)

		implicit none

		integer          :: n
		real, intent(in) :: x(:), y(:, :)
		real             :: trapzo_integ_mat(size(y, 2))

		n = size(x)
		trapzo_integ_mat = sum(spread(abs(x(2:) - x(:n-1)), 2, size(y, 2)) * (y(2:, :) + y(:n-1, :)), dim=1) * 0.5

	END FUNCTION

  SUBROUTINE GET_PATH(env_variable, join_with, path)

		implicit none

		integer                                    :: actual_length, status
		character(len=*), intent(in)               :: env_variable
		character(len=*), intent(in), optional     :: join_with
		character(len=256)                         :: long_string
		character(len=:), allocatable, intent(out) :: path

		CALL GET_ENVIRONMENT_VARIABLE(env_variable, long_string, actual_length, status)
		select case(status)
		case(-1)
			STOP "GET_PATH: the path is too long."
		case(+1)
			STOP "GET_PATH: environment variable does not exist."
		case(+2)
			STOP "GET_PATH: processor does not support environment variable."
		case default
			path = long_string(1:actual_length)
			if(path(actual_length:actual_length) /= "/") path = path//"/"
		end select

		if(present(join_with)) then
			if(join_with(1:1) == "/") then
				path = path//join_with(2:)
			else
				path = path//join_with
			end if
			actual_length = len(path)
			if(path(actual_length:actual_length) /= "/") path = path//"/"
		end if

	END SUBROUTINE

  SUBROUTINE PARSE_ARGS(short_opts, long_opts, opts, vals, args)
! This routine parse arguments given via command line. It works similar to getopt, if one argument
! starts with --, it is treated as a long-name option; if it starts with one -, is taken a short
! name option. Other wise, the argument goes to the args list.
!
! INPUT PARAMETERS:
! ----------------
!   * short_opts : character string value. If given, it should contain the possible options to be
!                  parsed. If one option requieres a value, it should be followed by a colon (:).
!   * long_opts  : character string vector. If given, it should contain the long versions of the
!                  possible options to be parsed. If one option requieres a value, it should be
!                  followed by a equal sign (=).
!
! OUTPUT PARAMETERS:
! -----------------
!   * opts       : character string vector, containing the name of the options given via command
!                  line.
!   * vals       : character string vector, containing the values of the options that requieres it.
!                  If one option doesn't requieres a value, the value in the corresponding position
!                  is set to blanks (" ").
!   * args       : character string vector containing the values of the arguments passed via command
!                  line.

		implicit none

		integer                                                :: i, j, k, narg, nopt, iopt, fopt
		character(len = *), intent(in), optional               :: short_opts, long_opts(:)
		character(len = *), allocatable, intent(out), optional :: opts(:), vals(:)
		character(len = *), allocatable, intent(out)           :: args(:)
		character(len = 256), allocatable                      :: args_(:), opts_(:), vals_(:)
		character(len = 256)                                   :: long_opt(1)
		character(len = :), allocatable                        :: short_opts_
		logical                                                :: has_val, must_have_val, invalid_opt
		logical                                                :: parsed, parse_short, parse_long
		logical, allocatable                                   :: mask(:), is_arg(:)

    parse_short = present(short_opts)
    parse_long  = present(long_opts)

    if((parse_short .or. parse_long) .and. .not. (present(opts) .and. present(vals))) then
      STOP "PARSE_ARGS: in order to parse options, 'opts' and 'vals' must be present."
    end if

		nopt = 0
		narg = COMMAND_ARGUMENT_COUNT()
		if(narg == 0) return

		if(parse_short) then
      nopt = len_trim(short_opts) - STRING_COUNT(short_opts, ":")

      allocate(short_opts_, source=short_opts//achar(0))
    end if
		if(parse_long) then
      nopt = nopt + size(long_opts)
      allocate(mask(size(long_opts)))
    end if

		allocate(is_arg(narg), args_(narg + 1))
    if(parse_short .or. parse_long) then
      allocate(opts_(nopt), vals_(nopt))

      opts_ = achar(0)
      vals_ = achar(0)
    end if

		args_  = achar(0)
		is_arg = .false.
		parsed = .false.

		do i = 1, narg
			CALL GET_COMMAND_ARGUMENT(i, args_(i))
		end do

		j = 0
		do i = 1, narg
			if(parsed) then
				parsed = .false.
				cycle
			end if

			if(parse_short .or. parse_long) j = j + 1
			if(j > nopt) STOP "PARSE_ARGS: size of options vector, 'nopt', is less than required."

			select case(index(args_(i)(:2), "-", back = .true.))
			case(2)
				if(any(is_arg)) then
					print*, "STOP PARSE_ARGS: invalid syntax. Expecting something like:"
					print*, "                 command [options] arguments"
					stop
				end if

				iopt = 3
				fopt = index(args_(i), "=", back = .true.) - 1

				has_val = fopt > 0
				if(.not. has_val) fopt = len_trim(args_(i))

        mask = long_opts == args_(i)(iopt:fopt+1)
        if(.not. any(mask)) then
          print*, "STOP PARSE_ARGS: invalid option '", args_(i)(iopt:fopt), "'."
          stop
        end if
				long_opt = pack(long_opts, mask)

				invalid_opt   = .not. any(mask)
				must_have_val = index(long_opt(1), "=") > 0

				if(invalid_opt) then
					print*, "STOP PARSE_ARGS: invalid option '", args_(i)(iopt:fopt), "'."
					stop
				else if(must_have_val .and. .not. has_val) then
					print*, "STOP PARSE_ARGS: option '", args_(i)(iopt:fopt), "' must have value."
					stop
				else if(.not. must_have_val .and. has_val) then
					print*, "STOP PARSE_ARGS: option '", args_(i)(iopt:fopt), "' must not have value."
					stop
				else
					opts_(j) = args_(i)(iopt:fopt)
					if(has_val) vals_(j) = args_(i)(fopt + 2:len_trim(args_(i)))
					if(.not. has_val) vals_(j) = ""
				end if

			case(1)
				if(any(is_arg)) then
					print*, "STOP PARSE_ARGS: invalid syntax. Expecting something like:"
					print*, "                 command [options] arguments"
					stop
				end if

				iopt = 2
				fopt = len_trim(args_(i))

				if(.not. parse_short) then
					print*, "STOP PARSE_ARGS: invalid option '", args_(i)(iopt:iopt), "'."
					stop
				end if

				do
					k = index(short_opts, args_(i)(iopt:iopt))
					invalid_opt   = k == 0
					must_have_val = short_opts_(k + 1:k + 1) == ":"

					if(invalid_opt) then
						print*, "STOP PARSE_ARGS: invalid option '", args_(i)(iopt:iopt), "'."
						stop
					else if(must_have_val .and. (iopt == fopt .and. args_(i + 1)(1:1) == "-")) then
						print*, "STOP PARSE_ARGS: option '", args_(i)(iopt:iopt), "' must have value."
						stop
					else
						opts_(j) = args_(i)(iopt:iopt)
						if(must_have_val .and. iopt < fopt) then
							vals_(j) = args_(i)(iopt + 1:fopt)
							iopt = fopt
						else if(must_have_val) then
							vals_(j) = trim(args_(i + 1))
							parsed   = .true.
						end if
						if(.not. must_have_val) vals_(j) = ""
					end if

					if(iopt < fopt) then
						iopt = iopt + 1
						j    = j + 1
					else
						exit
					end if
					if(j > nopt) STOP "PARSE_ARGS: size of options vector, 'nopt', is less than required."
				end do

			case(0)
				is_arg(i) = .true.
			end select

		end do

    if(parse_short .or. parse_long) then
      nopt = count(index(opts_, achar(0)) == 0)
      allocate(opts(nopt), vals(nopt), args(count(is_arg)))

      opts = pack(opts_, mask = index(opts_, achar(0)) == 0)
      vals = pack(vals_, mask = index(vals_, achar(0)) == 0)
    end if
    if(count(is_arg) > 0) then
      args = pack(args_(:narg), mask = is_arg)
    else
      deallocate(args)
    end if

	END SUBROUTINE

  SUBROUTINE MIN_CHI_SQUARE(observation, models, coeff1d, coeff2d, coeff3d, gen1d, gen2d, gen3d, chi_square)

		implicit none

		integer                        :: i, j, k, nobs, nmods
		integer, intent(out), optional :: gen1d(:), gen2d(:), gen3d(:)
		real(8), intent(in)            :: observation(:), models(:, :)
		real(8), intent(out), optional :: coeff1d(:), coeff2d(:), coeff3d(:)
		real(8), intent(out)           :: chi_square(3)
		real(8)                        :: a1d(1), a2d(2), a3d(3), current_chi_sq(3), det
		real(8), allocatable           :: comb_mat(:, :), comb_vec(:)
		real(8)                        :: m2d(2, 2), b2d(2), m3d(3, 3), b3d(3)
		logical                        :: solve_1d, solve_2d, solve_3d

		solve_1d = present(gen1d) .and. present(coeff1d)
		solve_2d = present(gen2d) .and. present(coeff2d)
		solve_3d = present(gen3d) .and. present(coeff3d)

		nobs  = size(observation)
		nmods = size(models, 2)

		if(count([solve_1d, solve_2d, solve_3d]) == 0) RETURN

		if(solve_1d) then
		if(size(gen1d) /= 1 .or. size(coeff1d) /= 1) STOP "MIN_CHI_SQUARE: gen1d (coeff1d) must have size = 1."
		end if

		if(solve_2d) then
		if(size(gen2d) /= 2 .or. size(coeff2d) /= 2) STOP "MIN_CHI_SQUARE: gen2d (coeff2d) must have size = 2."
		end if

		if(solve_3d) then
		if(size(gen3d) /= 3 .or. size(coeff3d) /= 3) STOP "MIN_CHI_SQUARE: gen3d (coeff3d) must have size = 3."
		end if

		if(nobs /= size(models, 1)) STOP "MIN_CHI_SQUARE: models and observation are not equally sampled."

		allocate(comb_mat(nmods, nmods), comb_vec(nmods))
		comb_mat   = 0.0
		comb_vec   = 0.0
		chi_square = huge(1.0)

		!$OMP PARALLEL DO
		do i = 1, nmods
			comb_mat(i, i) = sum(models(:, i) * models(:, i))
			comb_vec(i)    = sum(models(:, i) * observation)
			do j = i + 1, nmods
				comb_mat(i, j) = sum(models(:, i) * models(:, j))
				comb_mat(j, i) = comb_mat(i, j)
			end do
		end do
		!$OMP END PARALLEL DO

		if(solve_1d) then
			do i = 1, nmods

			a1d = comb_vec(i) / comb_mat(i, i)

			current_chi_sq(1) = sum((observation - a1d(1) * models(:, i)) ** 2) / (nobs - 2)
			if(current_chi_sq(1) < chi_square(1)) then
				chi_square(1) = current_chi_sq(1)
				coeff1d       = a1d
				gen1d         = i
			end if

			end do
		end if

		if(solve_2d) then
			do i = 1, nmods

			m2d(1, 1) = comb_mat(i, i)
			b2d(1)    = comb_vec(i)

			do j = i + 1, nmods

			m2d(1, 2) = comb_mat(i, j)
			m2d(2, 1) = comb_mat(j, i)
			m2d(2, 2) = comb_mat(j, j)
			b2d(2)    = comb_vec(j)

			det = m2d(1, 1) * m2d(2, 2) - m2d(1, 2) * m2d(2, 1)

			if(det /= 0.0) then
			a2d(1) = (b2d(1) * m2d(2, 2) - m2d(1, 2) * b2d(2)) / det
			a2d(2) = (m2d(1, 1) * b2d(2) - b2d(1) * m2d(2, 1)) / det

			if(all(a2d >= 0.0)) then
				current_chi_sq(2) = sum((observation - a2d(1) * models(:, i) - a2d(2) * models(:, j)) ** 2) / (nobs - 3)
				if(current_chi_sq(2) < chi_square(2)) then
					chi_square(2) = current_chi_sq(2)
					coeff2d       = a2d
					gen2d         = (/i, j/)
				end if
			end if
			end if

			end do
			end do
		end if

		if(solve_3d) then
			do i = 1, nmods

			m3d(1, 1) = comb_mat(i, i)
			b3d(1)    = comb_vec(i)

			do j = i + 1, nmods

			m3d(1, 2) = comb_mat(i, j)
			m3d(2, 1) = comb_mat(j, i)
			m3d(2, 2) = comb_mat(j, j)
			b3d(2)    = comb_vec(j)

			do k = j + 1, nmods

			m3d(1, 3) = comb_mat(i, k)
			m3d(3, 1) = comb_mat(k, i)
			m3d(2, 3) = comb_mat(j, k)
			m3d(3, 2) = comb_mat(k, j)
			m3d(3, 3) = comb_mat(k, k)
			b3d(3)    = comb_vec(k)

			det = m3d(1, 1) * m3d(2, 2) * m3d(3, 3) + &
			      m3d(2, 1) * m3d(3, 2) * m3d(1, 3) + &
			      m3d(3, 1) * m3d(1, 2) * m3d(2, 3) - &
			      m3d(1, 3) * m3d(2, 2) * m3d(3, 1) - &
			      m3d(2, 3) * m3d(3, 2) * m3d(1, 1) - &
			      m3d(3, 3) * m3d(1, 2) * m3d(2, 1)

			if(det /= 0) then
			a3d(1) = (b3d(1) * m3d(2, 2) * m3d(3, 3) + &
			          b3d(2) * m3d(3, 2) * m3d(1, 3) + &
			          b3d(3) * m3d(1, 2) * m3d(2, 3) - &
			          m3d(1, 3) * m3d(2, 2) * b3d(3) - &
			          m3d(2, 3) * m3d(3, 2) * b3d(1) - &
			          m3d(3, 3) * m3d(1, 2) * b3d(2)) / det
			a3d(2) = (m3d(1, 1) * b3d(2) * m3d(3, 3) + &
			          m3d(2, 1) * b3d(3) * m3d(1, 3) + &
			          m3d(3, 1) * b3d(1) * m3d(2, 3) - &
			          m3d(1, 3) * b3d(2) * m3d(3, 1) - &
			          m3d(2, 3) * b3d(3) * m3d(1, 1) - &
			          m3d(3, 3) * b3d(1) * m3d(2, 1)) / det
			a3d(3) = (m3d(1, 1) * m3d(2, 2) * b3d(3) + &
			          m3d(2, 1) * m3d(3, 2) * b3d(1) + &
			          m3d(3, 1) * m3d(1, 2) * b3d(2) - &
			          b3d(1) * m3d(2, 2) * m3d(3, 1) - &
			          b3d(2) * m3d(3, 2) * m3d(1, 1) - &
			          b3d(3) * m3d(1, 2) * m3d(2, 1)) / det

			if(all(a3d >= 0.0)) then
				current_chi_sq(3) = sum((observation - a3d(1) * models(:, i) - a3d(2) * models(:, j) - a3d(3) * models(:, k)) ** 2) / (nobs - 4)
				if(current_chi_sq(3) < chi_square(3)) then
					chi_square(3) = current_chi_sq(3)
					coeff3d       = a3d
					gen3d         = (/i, j, k/)
				end if
			end if
			end if

			end do
			end do
			end do
		end if

	END SUBROUTINE

	FUNCTION MEAN_VECTOR_(sample) RESULT(MEAN)

		implicit none

		real, intent(in) :: sample(:)
		real             :: mean

		mean = sum(sample) / size(sample)

	END FUNCTION

	FUNCTION MEAN_MATRIX_(sample) RESULT(MEAN)

		implicit none

		real, intent(in)  :: sample(:, :)
		real, allocatable :: mean(:)

		mean = sum(sample) / size(sample)

	END FUNCTION

	FUNCTION MEAN_MATRIX_REDUX_(sample, dim) RESULT(MEAN)

		implicit none

		integer, intent(in) :: dim
		real, intent(in)    :: sample(:, :)
		real, allocatable   :: mean(:)

		mean = sum(sample, dim) / size(sample, dim)

	END FUNCTION

	FUNCTION WEIGHTED_MEAN_VECTOR_(sample, weights) RESULT(MEAN)

		implicit none

		real, intent(in) :: sample(:), weights(:)
		real             :: mean

		if(sum(weights)==1.0) then
			mean = sum(weights * sample)
		else
			mean = sum(weights * sample) / sum(weights)
		end if

	END FUNCTION

	FUNCTION WEIGHTED_MEAN_MATRIX_(sample, weights) RESULT(MEAN)

		implicit none

		real, intent(in)  :: sample(:, :), weights(:, :)
		real, allocatable :: mean(:)

			if(sum(weights)==1.0) then
				mean = sum(weights * sample)
			else
				mean = sum(weights * sample) / sum(weights)
			end if

	END FUNCTION

	FUNCTION WEIGHTED_MEAN_MATRIX_REDUX_(sample, weights, dim) RESULT(MEAN)

		implicit none

		integer, intent(in) :: dim
		real, intent(in)    :: sample(:, :), weights(:, :)
		real, allocatable   :: mean(:)

			if(all(sum(weights, dim)==1.0)) then
				mean = sum(weights * sample, dim)
			else
				mean = sum(weights * sample, dim) / sum(weights, dim)
			end if

	END FUNCTION

  ELEMENTAL FUNCTION EVAL(string) RESULT(figure)

		implicit none

		integer                        :: figure
		character(len = *), intent(in) :: string

		read(string, *) figure

	END FUNCTION

	FUNCTION LISTING(string) RESULT(list)

		implicit none

		integer                        :: i, bracketing(2)
		integer, allocatable           :: list(:)
		character(len = *), intent(in) :: string
		character(len = 3)             :: implicit_

		if(index(string, ',...,') == 0) then
			write(6, *) 'LISTING: wrong format in string. Please enter "int1,...,int2", where int1 and'
			write(6, *) '         int2 are string representations of integer values satisfaying int1 < int2.'
			stop
		end if

		read(string, *) bracketing(1), implicit_, bracketing(2)

		if(bracketing(1) > bracketing(2)) STOP 'LISTING: in input string "int1,...,int2", integers must satisfy int1 < int2.'

		allocate(list(bracketing(2) - bracketing(1) + 1))
		list = [(i, i = bracketing(1), bracketing(2))]

	END FUNCTION

	FUNCTION SPLIT(string, substring) RESULT(string_vector)

		implicit none

		integer                              :: i, n, ascii_pos(10)
		character(*), intent(in)             :: string, substring
		character(len_trim(adjustl(string))) :: string_vector(:), copy
		character(1)                         :: replace_blank, replace_comma
		logical                              :: found_replacement, has_blanks, has_commas, last
		allocatable                          :: string_vector

		ascii_pos  = [35, 36, 37, 38, 42, 43, 45, 61, 94, 126]
		has_blanks = index(string, " ") /= 0
		has_commas = index(string, ",") /= 0
		copy       = trim(adjustl(string))
    n          = 0
    last       = .false.

		do i = 1, len(copy)
      if(copy(i:i) == substring .and. .not. last) then
        n    = n + 1
        last = .true.
      else if(copy(i:i) /= substring .and. last) then
        last = .false.
      end if
    end do
    n = n + 1
		allocate(string_vector(n))
		if(substring /= " " .and. substring /= ",") then

			if(has_blanks) then
				do i = 1, size(ascii_pos)
					replace_blank     = achar(ascii_pos(i))
					found_replacement = index(copy, replace_blank) == 0
					if(found_replacement) exit
				end do

				if(.not.found_replacement) STOP "SPLIT: replacement not found. Try including more ASCII positions."

				do while(scan(copy, " ") /= 0)
					copy = REPLACE(copy, [" ", replace_blank])
				end do
			end if

			if(has_commas) then
				do i = 1, size(ascii_pos)
					replace_comma     = achar(ascii_pos(i))
					found_replacement = index(copy, replace_comma) == 0
					if(found_replacement) exit
				end do

				if(.not.found_replacement) STOP "SPLIT: replacement not found. Try including more ASCII positions."

				do while(scan(copy, ",") /= 0)
					copy = REPLACE(copy, [",", replace_comma])
				end do
			end if

			do i = 1, n - 1
				copy = REPLACE(copy, [substring, ","//repeat(" ", len(substring) - 1)])
			end do

			read(copy, *) string_vector

			if(has_blanks) then
				do i = 1, n
					do while(index(string_vector(i), replace_blank) /= 0)
						string_vector(i) = REPLACE(string_vector(i), [replace_blank, " "])
					end do
				end do
			end if

			if(has_commas) then
				do i = 1, n
					do while(index(string_vector(i), replace_comma) /= 0)
						string_vector(i) = REPLACE(string_vector(i), [replace_comma, ","])
					end do
				end do
			end if
		else if(substring == " ") then
			if(has_commas) then
				do i = 1, size(ascii_pos)
					replace_comma     = achar(ascii_pos(i))
					found_replacement = index(copy, replace_comma) == 0
					if(found_replacement) exit
				end do

				if(.not.found_replacement) STOP "SPLIT: replacement not found. Try including more ASCII positions."

				do while(scan(copy, ",") /= 0)
					copy = REPLACE(copy, [",", replace_comma])
				end do
			end if

			read(copy, *) string_vector

			if(has_commas) then
				do i = 1, n
					do while(index(string_vector(i), replace_comma) /= 0)
						string_vector(i) = REPLACE(string_vector(i), [replace_comma, ","])
					end do
				end do
			end if
		else if(substring == ",") then
			if(has_blanks) then
				do i = 1, size(ascii_pos)
					replace_blank     = achar(ascii_pos(i))
					found_replacement = index(copy, replace_blank) == 0
					if(found_replacement) exit
				end do

				if(.not.found_replacement) STOP "SPLIT: replacement not found. Try including more ASCII positions."

				do while(scan(copy, " ") /= 0)
					copy = REPLACE(copy, [" ", replace_blank])
				end do
			end if

			read(copy, *) string_vector

			if(has_blanks) then
				do i = 1, n
					do while(index(string_vector(i), replace_blank) /= 0)
						string_vector(i) = REPLACE(string_vector(i), [replace_blank, " "])
					end do
				end do
			end if
		end if

	END FUNCTION

	FUNCTION STRING(figure)
! This routines perform the conversion between number and string.
!
! INPUT PARAMETERS:
! ----------------
!   * figure : integer number to be converted.
!
! OUTPUT PARAMETERS:
! -----------------
!   * string : character string representation of the input figure.

		implicit none

		integer, intent(in)              :: figure
		integer                          :: status, clen
		character(len = :), allocatable  :: string
		character(len = 256)             :: aux

		write(aux, *, iostat = status) figure
		if(status == 2) STOP "STRING: size of the buffer is not long enough."

		aux  = adjustl(aux)
		clen = len_trim(aux)

		string = aux(:clen)

	END FUNCTION

	PURE FUNCTION STRING_COUNT(string, substring) RESULT(scount)
! This routine returns the number of ocurrences of a substring in a string.
!
! INPUT PARAMETERS:
! ----------------
!   * string    : character string with the string in which to search for the substring.
!   * substring : character string that will be searched for in the gicen string. If is not
!                 contained in string, the result is cero.
!
! OUTPUT PARAMETERS:
! -----------------
!   * scount    : integer value of the number of occurrences of substring in string.

		implicit none

		integer                           :: i
		integer                           :: scount
		character(len = *), intent(in)    :: string, substring
		character(len = len_trim(string)) :: string_

		scount  = 0
		string_ = string
		i       = index(string_, substring)
		do while(i /= 0)
			scount  = scount + 1
			string_ = string_(i + len(substring):)
			i       = index(trim(string_), substring)
		end do

	END FUNCTION

	SUBROUTINE STRING_LOC(string, substring, position)
! This routine perform the location of a subtring within a vector string.
!
! INPUT PARAMETERS:
! ----------------
!   * string    : character string vector.
!   * substring : character string which will be searched for within the given vector.
!
! OUTPUT PARAMETERS:
! -----------------
!   * position  : integer number with the position of substring in the given vector.

		implicit none

		integer                        :: i
		integer, intent(out)           :: position
		character(len = *), intent(in) :: string(:), substring

		position = 0
		do i = 1, size(string)
			if(index(trim(string(i)), trim(substring)) /= 0) then
				position = i
				return
			end if
		end do

	END SUBROUTINE

  FUNCTION REPLACE(string, substrings)

    implicit none

    integer                         :: i
    character(len = *), intent(in)  :: string, substrings(2)
    character(len = :), allocatable :: replace

    i = index(string, trim(substrings(1)))
    if(i == 1) then
      replace = trim(substrings(2))//string(i + len_trim(substrings(1)):)
    else if(i == len(string)) then
      replace = string(:i - 1)//trim(substrings(2))
    else
      replace = string(:i - 1)//trim(substrings(2))//string(i + len_trim(substrings(1)):)
    end if

  END FUNCTION

! =================================================================================================!
! ======================================== APPEND ROUTINES ========================================!
! =================================================================================================!
! This routines perform the append action using the following possibles inputs:
!   matrix1, matrix2
!   matrix, vector
!   vector1, vector2
!
! INPUT PARAMETERS:
! ----------------
!   * matrix(1)/vector(1) : it could be a one/two-dimensional array of integer/real type. If is not
!                           allocated, then this is equivalent to copy information from the second
!                           parameter to the first.
!   * matrix(2)/vector(2) : if first parameter is a matrix this could be another matrix or a vector.
!                           Otherwise, this must be a vector, matching data type with first
!                           parameter in both cases.
!
! OUTPUT PARAMETERS:
! -----------------
!   * matrix(1)/vector(1) : first parameter with second parameter appended to it. Note that first
!                           parameter must be a dynamic array, while the second could be a constant
!                           array.

	SUBROUTINE APPEND_MATRIX_MATRIX_CHARACTER(matrix1, matrix2)

		implicit none

		integer                                  :: nrow_a, ncol_a, nrow_b, ncol_b
		character(*), allocatable, intent(inout) :: matrix1(:, :)
		character(*), intent(in)                 :: matrix2(:, :)
		character(len(matrix1)), allocatable     :: aux(:, :)

		if(len(matrix1)<len(matrix2)) STOP "APPEND: lengths must fulfill: len(matrix1) >= len(matrix2)."

		nrow_b = size(matrix2, 1)
		ncol_b = size(matrix2, 2)
		if(allocated(matrix1)) then
			nrow_a = size(matrix1, 1)
			ncol_a = size(matrix1, 2)
			if(ncol_a /= ncol_b) STOP "APPEND: arrays must have same number of columns."
		else
			allocate(matrix1(nrow_b, ncol_b))
			matrix1 = matrix2
			return
		end if

		allocate(aux(nrow_a + nrow_b, ncol_a))
		aux(:nrow_a, :)     = matrix1
		aux(nrow_a + 1:, :) = matrix2
		deallocate(matrix1)
		allocate(matrix1(nrow_a + nrow_b, ncol_a))
		matrix1 = aux

	END SUBROUTINE

	SUBROUTINE APPEND_MATRIX_VECTOR_CHARACTER(matrix, vector)

		implicit none

		integer                                  :: nrow_a, ncol_a, nrow_b, ncol_b
		character(*), allocatable, intent(inout) :: matrix(:, :)
		character(*), intent(in)                 :: vector(:)
		character(len(matrix)), allocatable      :: aux(:, :)

		if(len(matrix)<len(vector)) STOP "APPEND: lengths must fulfill: len(matrix) >= len(vector)."

		nrow_b = 1
		ncol_b = size(vector)
		if(allocated(matrix)) then
			nrow_a = size(matrix, 1)
			ncol_a = size(matrix, 2)
			if(ncol_a /= ncol_b) STOP "APPEND: arrays must have same number of columns."
		else
			allocate(matrix(nrow_b, ncol_b))
			matrix(1, :) = vector
			return
		end if

		allocate(aux(nrow_a + nrow_b, ncol_a))
		aux(:nrow_a, :)    = matrix
		aux(nrow_a + 1, :) = vector
		deallocate(matrix)
		allocate(matrix(nrow_a + nrow_b, ncol_a))
		matrix = aux

	END SUBROUTINE

	SUBROUTINE APPEND_VECTOR_VECTOR_CHARACTER(vector1, vector2)

		implicit none

		integer                                  :: nrow_a, nrow_b
		character(*), allocatable, intent(inout) :: vector1(:)
		character(*), intent(in)                 :: vector2(:)
		character(len(vector1)), allocatable     :: aux(:)

		if(len(vector1)<len(vector2)) STOP "APPEND: lengths must fulfill: len(vector1) >= len(vector2)."

		nrow_b = size(vector2)
		if(allocated(vector1)) then
			nrow_a = size(vector1)
		else
			allocate(vector1(nrow_b))
			vector1 = vector2
			return
		end if

		allocate(aux(nrow_a + nrow_b))
		aux(:nrow_a)     = vector1
		aux(nrow_a + 1:) = vector2
		deallocate(vector1)
		allocate(vector1(nrow_a + nrow_b))
		vector1 = aux

	END SUBROUTINE

	SUBROUTINE APPEND_MATRIX_MATRIX_INTEGER(matrix1, matrix2)

		implicit none

		integer                             :: nrow_a, ncol_a, nrow_b, ncol_b
		integer, allocatable, intent(inout) :: matrix1(:, :)
		integer, intent(in)                 :: matrix2(:, :)
		integer, allocatable                :: aux(:, :)

		nrow_b = size(matrix2, 1)
		ncol_b = size(matrix2, 2)
		if(allocated(matrix1)) then
			nrow_a = size(matrix1, 1)
			ncol_a = size(matrix1, 2)
			if(ncol_a/=ncol_b) STOP "APPEND: arrays must have same number of columns."
		else
			allocate(matrix1(nrow_b, ncol_b))
			matrix1 = matrix2
			return
		end if

		allocate(aux(nrow_a + nrow_b, ncol_a))
		aux(:nrow_a, :)     = matrix1
		aux(nrow_a + 1:, :) = matrix2
		deallocate(matrix1)
		allocate(matrix1(nrow_a + nrow_b, ncol_a))
		matrix1 = aux

	END SUBROUTINE

	SUBROUTINE APPEND_MATRIX_VECTOR_INTEGER(matrix, vector)

		implicit none

		integer                             :: nrow_a, ncol_a, nrow_b, ncol_b
		integer, allocatable, intent(inout) :: matrix(:, :)
		integer, intent(in)                 :: vector(:)
		integer, allocatable                :: aux(:, :)

		nrow_b = 1
		ncol_b = size(vector)
		if(allocated(matrix)) then
			nrow_a = size(matrix, 1)
			ncol_a = size(matrix, 2)
			if(ncol_a/=ncol_b) STOP "APPEND: arrays must have same number of columns."
		else
			allocate(matrix(nrow_b, ncol_b))
			matrix(1, :) = vector
			return
		end if

		allocate(aux(nrow_a + nrow_b, ncol_a))
		aux(:nrow_a, :)     = matrix
		aux(nrow_a + 1, :) = vector
		deallocate(matrix)
		allocate(matrix(nrow_a + nrow_b, ncol_a))
		matrix = aux

	END SUBROUTINE

	SUBROUTINE APPEND_VECTOR_VECTOR_INTEGER(vector1, vector2)

		implicit none

		integer                             :: nrow_a, nrow_b
		integer, allocatable, intent(inout) :: vector1(:)
		integer, intent(in)                 :: vector2(:)
		integer, allocatable                :: aux(:)

		nrow_b = size(vector2)
		if(allocated(vector1)) then
			nrow_a = size(vector1)
		else
			allocate(vector1(nrow_b))
			vector1 = vector2
			return
		end if

		allocate(aux(nrow_a + nrow_b))
		aux(:nrow_a)     = vector1
		aux(nrow_a + 1:) = vector2
		deallocate(vector1)
		allocate(vector1(nrow_a + nrow_b))
		vector1 = aux

	END SUBROUTINE

	SUBROUTINE APPEND_MATRIX_MATRIX_REAL(matrix1, matrix2)

		implicit none

		INTEGER                          :: nrow_a, ncol_a, nrow_b, ncol_b
		real, allocatable, intent(inout) :: matrix1(:, :)
		real, intent(in)                 :: matrix2(:, :)
		real, allocatable                :: aux(:, :)


		nrow_b = size(matrix2, 1)
		ncol_b = size(matrix2, 2)
		if(allocated(matrix1)) then
			nrow_a = size(matrix1, 1)
			ncol_a = size(matrix1, 2)
			if(ncol_a/=ncol_b) stop "APPEND ERROR: arrays must have same number of columns."
		else
			allocate(matrix1(nrow_b, ncol_b))
			matrix1 = matrix2
			return
		end if

		allocate(aux(nrow_a + nrow_b, ncol_a))
		aux(:nrow_a, :)     = matrix1
		aux(nrow_a + 1:, :) = matrix2
		deallocate(matrix1)
		allocate(matrix1(nrow_a + nrow_b, ncol_a))
		matrix1 = aux

	END SUBROUTINE

	SUBROUTINE APPEND_MATRIX_VECTOR_REAL(matrix, vector)

		implicit none

		INTEGER                          :: nrow_a, ncol_a, nrow_b, ncol_b
		real, allocatable, intent(inout) :: matrix(:, :)
		real, intent(in)                 :: vector(:)
		real, allocatable                :: aux(:, :)


		nrow_b = 1
		ncol_b = size(vector)
		if(allocated(matrix)) then
			nrow_a = size(matrix, 1)
			ncol_a = size(matrix, 2)
			if(ncol_a/=ncol_b) stop "APPEND ERROR: arrays must have same number of columns."
		else
			allocate(matrix(nrow_b, ncol_b))
			matrix(1, :) = vector
			return
		end if

		allocate(aux(nrow_a + nrow_b, ncol_a))
		aux(:nrow_a, :)    = matrix
		aux(nrow_a + 1, :) = vector
		deallocate(matrix)
		allocate(matrix(nrow_a + nrow_b, ncol_a))
		matrix = aux

	END SUBROUTINE

	SUBROUTINE APPEND_VECTOR_VECTOR_REAL(vector1, vector2)

		implicit none

		INTEGER                          :: nrow_a, nrow_b
		real, allocatable, intent(inout) :: vector1(:)
		real, intent(in)                 :: vector2(:)
		real, allocatable                :: aux(:)

		nrow_b = size(vector2)
		if(allocated(vector1)) then
			nrow_a = size(vector1)
		else
			allocate(vector1(nrow_b))
			vector1 = vector2
			return
		end if

		allocate(aux(nrow_a + nrow_b))
		aux(:nrow_a)     = vector1
		aux(nrow_a + 1:) = vector2
		deallocate(vector1)
		allocate(vector1(nrow_a + nrow_b))
		vector1 = aux

	END SUBROUTINE

!	=================================================================================================!
!	====================================== SUBMATRIX ROUTINES =======================================!
!	=================================================================================================!
! This routines extract selected rows/columns from a given matrix.
!
!	INPUT PARAMETERS:
!	----------------
!   * matrix(1) : two-dimensional array of integer/real type.
!   * dim       : integer number. It could take 1 (for row) or 2 (for column) as value.
!   * indexes   : scalar integer index of one row/column to extract. Is possible to extract
!                 several rows/columns if an integer vector of indexes is given.
!
! OUTPUT PARAMETERS:
! -----------------
!   * matrix2/vector : one/two-dimensional array with the corresponding submatrix extracted from the
!                      given matrix. Data type must match that of the first parameter.

  SUBROUTINE SUBMATRIX_SCALAR_INTEGER(matrix, dim, indexes, vector)

  	implicit none

  	integer, intent(in)               :: dim, indexes
  	integer, intent(in)               :: matrix(:, :)
  	integer, allocatable, intent(out) :: vector(:)

  	select case(dim)
  	case(1)

  		allocate(vector(size(matrix, 2)))
  		vector = matrix(indexes, :)

  	case(2)

  		allocate(vector(size(matrix, 1)))
  		vector = matrix(:, indexes)

  	case default
  		STOP "GET_VECTOR: invalid dim value. Allowed values are 1 and 2."
  	end select

  END SUBROUTINE

  SUBROUTINE SUBMATRIX_SCALAR_REAL(matrix, dim, indexes, vector)

  	implicit none

  	integer, intent(in)            :: dim, indexes
  	real, intent(in)               :: matrix(:, :)
  	real, allocatable, intent(out) :: vector(:)

  	select case(dim)
  	case(1)

  		allocate(vector(size(matrix, 2)))
  		vector = matrix(indexes, :)

  	case(2)

  		allocate(vector(size(matrix, 1)))
  		vector = matrix(:, indexes)

  	case default
  		STOP "GET_VECTOR: invalid dim value. Allowed values are 1 and 2."
  	end select

  END SUBROUTINE

  SUBROUTINE SUBMATRIX_VECTOR_INTEGER(matrix1, dim, indexes, matrix2)

  	implicit none

  	integer                           :: i, n, ndim
  	integer, intent(in)               :: dim, indexes(:)
  	integer, intent(in)               :: matrix1(:, :)
  	integer, allocatable, intent(out) :: matrix2(:, :)

  	n = size(indexes)

  	select case(dim)
  	case(1)

  		ndim = size(matrix1, 2)
  		allocate(matrix2(n, ndim))
  		forall(i = 1:n) matrix2(i, :) = matrix1(indexes(i), :)

  	case(2)

  		ndim = size(matrix1, 1)
  		allocate(matrix2(ndim, n))
  		forall(i = 1:n) matrix2(:, i) = matrix1(:, indexes(i))

  	case default
  		STOP "GET_VECTOR: invalid dim value. Allowed values are 1 and 2."
  	end select

  END SUBROUTINE

  SUBROUTINE SUBMATRIX_VECTOR_REAL(matrix1, dim, indexes, matrix2)

  	implicit none

  	integer                        :: i, n, ndim
  	integer, intent(in)            :: dim, indexes(:)
  	real, intent(in)               :: matrix1(:, :)
  	real, allocatable, intent(out) :: matrix2(:, :)

  	n = size(indexes)

  	select case(dim)
  	case(1)

  		ndim = size(matrix1, 2)
  		allocate(matrix2(n, ndim))
  		forall(i = 1:n) matrix2(i, :) = matrix1(indexes(i), :)

  	case(2)

  		ndim = size(matrix1, 1)
  		allocate(matrix2(ndim, n))
  		forall(i = 1:n) matrix2(:, i) = matrix1(:, indexes(i))

  	case default
  		STOP "GET_VECTOR: invalid dim value. Allowed values are 1 and 2."
  	end select

  END SUBROUTINE

  SUBROUTINE SUBMATRIX_INPLACE_SCALAR_INTEGER(matrix, dim, indexes)

  	implicit none

  	integer, intent(in)                 :: dim, indexes
  	integer, allocatable, intent(inout) :: matrix(:, :)
  	integer, allocatable                :: aux(:)

  	select case(dim)
  	case(1)

  		allocate(aux(size(matrix, 2)))
  		aux = matrix(indexes, :)

  		deallocate(matrix)
  		allocate(matrix(1, size(aux)))
  		matrix(1, :) = aux

  	case(2)

  		allocate(aux(size(matrix, 1)))
  		aux = matrix(:, indexes)

  		deallocate(matrix)
  		allocate(matrix(size(aux), 1))
  		matrix(:, 1) = aux

  	case default
  		STOP "GET_VECTOR: invalid dim value. Allowed values are 1 and 2."
  	end select

  END SUBROUTINE

  SUBROUTINE SUBMATRIX_INPLACE_SCALAR_REAL(matrix, dim, indexes)

  	implicit none

  	integer, intent(in)              :: dim, indexes
  	real, allocatable, intent(inout) :: matrix(:, :)
  	real, allocatable                :: aux(:)

  	select case(dim)
  	case(1)

  		allocate(aux(size(matrix, 2)))
  		aux = matrix(indexes, :)

  		deallocate(matrix)
  		allocate(matrix(1, size(aux)))
  		matrix(1, :) = aux

  	case(2)

  		allocate(aux(size(matrix, 1)))
  		aux = matrix(:, indexes)

  		deallocate(matrix)
  		allocate(matrix(size(aux), 1))
  		matrix(:, 1) = aux

  	case default
  		STOP "GET_VECTOR: invalid dim value. Allowed values are 1 and 2."
  	end select

  END SUBROUTINE

  SUBROUTINE SUBMATRIX_INPLACE_VECTOR_INTEGER(matrix, dim, indexes)

  	implicit none

  	integer                             :: i, n, ndim
  	integer, intent(in)                 :: dim, indexes(:)
  	integer, allocatable, intent(inout) :: matrix(:, :)
  	integer, allocatable                :: aux(:, :)

  	n = size(indexes)

  	select case(dim)
  	case(1)

  		ndim = size(matrix, 2)
  		allocate(aux(n, ndim))
  		forall(i=1:n) aux(i, :) = matrix(indexes(i), :)

  		deallocate(matrix)
  		allocate(matrix(n, ndim))
  		matrix = aux

  	case(2)

  		ndim = size(matrix, 1)
  		allocate(aux(ndim, n))
  		forall(i=1:n) aux(:, i) = matrix(:, indexes(i))

  		deallocate(matrix)
  		allocate(matrix(ndim, n))
  		matrix = aux

  	case default
  		STOP "GET_VECTOR: invalid dim value. Allowed values are 1 and 2."
  	end select

  END SUBROUTINE

  SUBROUTINE SUBMATRIX_INPLACE_VECTOR_REAL(matrix, dim, indexes)

  	implicit none

  	integer                          :: i, n, ndim
  	integer, intent(in)              :: dim, indexes(:)
  	real, allocatable, intent(inout) :: matrix(:, :)
  	real, allocatable                :: aux(:, :)

  	n = size(indexes)

  	select case(dim)
  	case(1)

  		ndim = size(matrix, 2)
  		allocate(aux(n, ndim))
  		forall(i=1:n) aux(i, :) = matrix(indexes(i), :)

  		deallocate(matrix)
  		allocate(matrix(n, ndim))
  		matrix = aux

  	case(2)

  		ndim = size(matrix, 1)
  		allocate(aux(ndim, n))
  		forall(i=1:n) aux(:, i) = matrix(:, indexes(i))

  		deallocate(matrix)
  		allocate(matrix(ndim, n))
  		matrix = aux

  	case default
  		STOP "GET_VECTOR: invalid dim value. Allowed values are 1 and 2."
  	end select

  END SUBROUTINE

  SUBROUTINE SUBVECTOR_INPLACE_VECTOR_REAL(vector, indexes)

  	implicit none

  	integer                          :: i, n
  	integer, intent(in)              :: indexes(:)
  	real, allocatable, intent(inout) :: vector(:)
  	real, allocatable                :: aux(:)

  	n = size(indexes)

  	allocate(aux(n))
  	forall(i=1:n) aux(i) = vector(indexes(i))

  	deallocate(vector)
  	allocate(vector(n))
  	vector = aux

  END SUBROUTINE

  !	=================================================================================================!
  !	======================================== SLICE ROUTINES =========================================!
  !	=================================================================================================!
  ! This routines extract a submatrix of contiguous rows/columns from a given matrix or contiguous
  ! values from a given vector.
  !
  !	INPUT PARAMETERS:
  !	----------------
  !   * matrix/vector         : one/two-dimensional array of integer/real type.
  !   * dim             : integer number. It could take 1 (for row) or 2 (for column) as value.
  !   * starts_with/ends_with : integer number corresponding to the index of the first/last row or
  !                             column (if a matrix is given) or the first/last value (if a vector is
  !                             given) to extract.
  !
  ! OUTPUT PARAMETERS:
  ! -----------------
  !   * matrix/vector         : one/two-dimensional array with the corresponding submatrix extracted
  !                             from the given matrix. Data type must match that of the first
  !                             parameter.

  SUBROUTINE SLICE_MATRIX_INTEGER(matrix, dim, starts_with, ends_with)

  	implicit none

  	integer                             :: nrow, ncol
  	integer,intent(in)                  :: dim, starts_with, ends_with
  	integer, allocatable, intent(inout) :: matrix(:, :)
  	integer, allocatable                :: aux(:, :)

  	select case(dim)
  	case(1)

  	nrow = ends_with - starts_with + 1
  	ncol = size(matrix,2)

  	allocate(aux(nrow,ncol))
  	aux = matrix(starts_with:ends_with, :)
  	deallocate(matrix)
  	allocate(matrix(nrow, ncol))

  	case(2)
  		nrow = size(matrix, 1)
  		ncol = ends_with - starts_with + 1

  		allocate(aux(nrow, ncol))
  		aux = matrix(:, starts_with:ends_with)
  		deallocate(matrix)
  		allocate(matrix(nrow, ncol))

  	case default
  		STOP "SLICE: invalid dim value. Allowed values are 1 and 2."
  	end select

  	matrix = aux

  END SUBROUTINE

  SUBROUTINE SLICE_MATRIX_REAL(matrix, dim, starts_with, ends_with)

  	implicit none

  	integer                          :: nrow, ncol
  	integer,intent(in)               :: dim, starts_with, ends_with
  	real, allocatable, intent(inout) :: matrix(:, :)
  	real, allocatable                :: aux(:, :)

  	select case(dim)
  	case(1)

  	nrow = ends_with - starts_with + 1
  	ncol = size(matrix,2)

  	allocate(aux(nrow,ncol))
  	aux = matrix(starts_with:ends_with, :)
  	deallocate(matrix)
  	allocate(matrix(nrow, ncol))

  	case(2)
  		nrow = size(matrix, 1)
  		ncol = ends_with - starts_with + 1

  		allocate(aux(nrow, ncol))
  		aux = matrix(:, starts_with:ends_with)
  		deallocate(matrix)
  		allocate(matrix(nrow, ncol))

  	case default
  		STOP "SLICE: invalid dim value. Allowed values are 1 and 2."
  	end select

  	matrix = aux

  END SUBROUTINE

  SUBROUTINE SLICE_VECTOR_INTEGER(vector, starts_with, ends_with)

  	implicit none

  	integer                             :: n
  	integer,intent(in)                  :: starts_with, ends_with
  	integer, allocatable                :: aux(:)
  	integer, allocatable, intent(inout) :: vector(:)

  	n = ends_with - starts_with + 1

  	allocate(aux(n))
  	aux = vector(starts_with:ends_with)
  	deallocate(vector)
  	allocate(vector(n))

  	vector = aux

  END SUBROUTINE

  SUBROUTINE SLICE_VECTOR_REAL(vector, starts_with, ends_with)

  	implicit none

  	integer                          :: n
  	integer,intent(in)               :: starts_with, ends_with
  	real, allocatable                :: aux(:)
  	real, allocatable, intent(inout) :: vector(:)

  	n = ends_with - starts_with + 1

  	allocate(aux(n))
  	aux = vector(starts_with:ends_with)
  	deallocate(vector)
  	allocate(vector(n))

  	vector = aux

  END SUBROUTINE

!	=================================================================================================!
!	=================================== NEAREST_ELEMENT ROUTINES ====================================!
!	=================================================================================================!
! This routines locates the element in vector nearest to a given value.
!
!	INPUT PARAMETERS:
!	----------------
!   * vector   : one-dimensional array of integer/real type.
!   * value    : scalar (matching type of the first parameter) to locate within the given vector.
!
! OUTPUT PARAMETERS:
! -----------------
!   * position : this is the position of the element in vector nearest to value.

	SUBROUTINE NEAREST_ELEMENT_INTEGER(vector, value, position)

		implicit none

		integer, intent(out) :: position
		integer, intent(in)  :: vector(:), value

		position = minloc(abs(vector-value), dim=1)

	END SUBROUTINE

	SUBROUTINE NEAREST_ELEMENT_REAL(vector, value, position)

		implicit none

		integer, intent(out) :: position
		real, intent(in)     :: vector(:), value

		position = minloc(abs(vector-value), dim=1)

	END SUBROUTINE

!	=================================================================================================!
!	===================================== VECTOR_LOC ROUTINES =======================================!
!	=================================================================================================!
! This routines locates a value/vector within a vector using the indexes of the nearest value(s).
!
!	INPUT PARAMETERS:
!	----------------
!   * domain   : one-dimensional array of integer/real type.
!   * span    : scalar or one-dimensional array (matching type of the first parameter) to locate
!                within the given domain.
!
! OUTPUT PARAMETERS:
! -----------------
!   * position : if scalar is given as second parameter, this is a scalar with the index of
!                the closest value to scalar in domain. If array is given instead, this is a vector
!                containing two indexes i, j, such that:
!
!                     domain(i) >= span(1) and domain(j) <= span(n)
!
!                where n = size(span). If any boundary of span is outside of domain, the
!                corresponding domain boundary is given instead.

	SUBROUTINE VECTOR_LOC_INTEGER(domain, span, position)

		implicit none

		integer              :: n
		integer, intent(out) :: position(2)
		integer, intent(in)  :: domain(:), span(:)

		n = size(span)
		position(1) = maxloc(domain, dim=1, mask=domain-span(1)<=0)
		position(2) = minloc(domain, dim=1, mask=domain-span(n)>=0)

	END SUBROUTINE

	SUBROUTINE VECTOR_LOC_REAL(domain, span, position)

		implicit none

		integer              :: n
		integer, intent(out) :: position(2)
		real, intent(in)     :: domain(:), span(:)

		n = size(span)
		position(1) = maxloc(domain, dim=1, mask=domain-span(1)<=0.0)
		position(2) = minloc(domain, dim=1, mask=domain-span(n)>=0.0)

	END SUBROUTINE

  SUBROUTINE VELOCITY_DISP(wlength, fssp, fluxo, sigma, wl_range, losvd_range, vel_disp)
! This routine computes LOSVD given an observed sed and a set of models. The method consist on
! fitting a given range of the sed using the set of models via MLE using the NNLS for computing the
! minimum chi square and then using an implementation of the Brent's method for finding the value of
! LOSVD correspoding to the best fit.
!
! DEPENDENCES:
! -----------
!   * BRENT    : This routine of NR(F77) for finding minima of a given one-variable function.
!   * NNLS     : This routine can be found in the GASPEX's codes. It finds the Least Square solution
!                for fitting a model to a given data. The model is a linear combination of some
!                given ingredients and the coefficients are restricted to be > 0, therefore the name
!                NonNegative Least Square.
!
! INPUT PARAMETERS:
! ----------------
!   * lambda   : real vector array containing the wavelengths of both, the problem sed and the
!                models.
!   * fssp     : real matrix array containing the fluxes of the models.
!   * fluxo    : real vector array containing the fluxes of the problem sed.
!   * sigma    : real vector array containing the dispersions in the fluxes of the problem sed.
!   * wl_range : real vector array containing the range of wavelength for computing the LOSVD.
!
! OUTPUT PARAMETERS:
! -----------------
!   * vel_disp : real value of the computed LOSVD.

    implicit none

    real, intent(in)  :: wlength(:), fssp(:, :), fluxo(:), sigma(:)
    real, intent(in)  :: wl_range(2), losvd_range(2)
    real, intent(out) :: vel_disp
    real              :: best_chi
    real, external    :: BRENT

    if(size(fluxo) /= size(fssp, 1)) then
      STOP "VELOCITY_DISP: models and observation must be equally sampled."
    end if

    best_chi = BRENT(losvd_range(1), sum(losvd_range)/2, losvd_range(2), CHI_SQUARE, 1e-4, vel_disp)

    CONTAINS
    FUNCTION CHI_SQUARE(losvd) RESULT(CHISQ)

      implicit none

      integer                            :: i, k(2), nnls_flag
      integer                            :: nages, nwlo, nwlm
      integer                            :: ni
      integer, allocatable               :: nj(:), iwlo(:)
      real, intent(in)                   :: losvd
      real, dimension(:, :), allocatable :: fssp_vd, fssp_vdsig
      real, dimension(:), allocatable    :: fssp_j, fluxo_s, ap, vect1, vect2, sed_fit
      real                               :: res
      real                               :: CHISQ
      logical, allocatable               :: mask(:)

      CALL VECTOR_LOC(wlength, wl_range, k)

      nwlo  = k(2) - k(1) + 1
      nages = size(fssp, 2)

      allocate(mask(nwlo))

      mask = fluxo(k(1):k(2)) > 0.0
      iwlo = [(i, i=1, nwlo)]
      iwlo = pack(iwlo, mask)
      nwlm = size(iwlo)

      allocate(nj(nwlo), ap(nages), vect1(nwlo), vect2(nwlo))
      allocate(fssp_j(nwlo), fssp_vd(nwlo, nages))
      allocate(fssp_vdsig(nwlm, nages), fluxo_s(nwlm), sed_fit(nwlo))

      fssp_vd = fssp(k(1):k(2), :)
      CALL APPLY_VD(wlength(k(1):k(2)), fssp_vd, losvd)
      fssp_vdsig = fssp_vd(iwlo, :) / spread(pack(sigma(k(1):k(2)), mask), 2, nages)
      fluxo_s    = pack(fluxo(k(1):k(2)) / sigma(k(1):k(2)), mask)

      CALL NNLS(fssp_vdsig, nwlm, nwlm, nages, fluxo_s, ap, res, vect1, vect2, nj, nnls_flag, ni)

      sed_fit = sum(spread(ap, 1, nwlo) * fssp_vd, dim=2)
      CHISQ = sum(((fluxo(k(1):k(2))-sed_fit)/sigma(k(1):k(2)))**2, mask=mask)
      CHISQ = CHISQ / (nwlm-count(ap>0.0)-1)

    END FUNCTION

  END SUBROUTINE

  SUBROUTINE FILTER_REDUNDANT(observation, models, tolerance, not_redundant)

    implicit none

    integer             :: i, k, nmod
    real(8), intent(in) :: observation(:), models(:, :), tolerance
    real(8)             :: F_i, F_j
    integer,allocatable :: not_redundant(:)

    nmod = size(models, 2)

    allocate(not_redundant(nmod))
    not_redundant = 0

    F_j = dot_product(observation, models(:, 1))

    k = 1
    not_redundant(k) = 1
    do i=2, nmod
      F_i = dot_product(observation, models(:, i))

      if(abs(F_i-F_j)>=tolerance*F_i) then
        F_j = F_i
        k = k + 1
      end if
      not_redundant(k) = i
    end do

    not_redundant = pack(not_redundant, mask=not_redundant>0)

  END SUBROUTINE

END MODULE
