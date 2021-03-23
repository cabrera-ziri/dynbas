PROGRAM DynBaS

  USE DYNBAS_PIECES_MODULE

  IMPLICIT NONE

  integer                     :: i, j, k, n, m
  integer                     :: unit, iostatus, clstatus
  integer                     :: nwef, nwlm, nwlo, nage, step, nmet, nmod, nbas, navs, nrsh, ind(3)
  integer                     :: nhdr, ncol
  integer                     :: it1d(1), it2d(2), it3d(3), passband
  integer, allocatable        :: filter_ids(:), iz(:), gen1d(:, :), gen2d(:, :), gen3d(:, :), ibasis(:)
  real, allocatable           :: redshifts(:), uni_ages(:)
  real, allocatable           :: t(:), ages(:), wlength(:), mods(:, :, :)
  real, allocatable           :: basis(:, :), fssps(:, :), Zs(:), Zs_(:)
  real, allocatable           :: weff(:), wlo(:), fluxo_s(:), sigma_s(:), dyn_flux(:, :)
  real, allocatable           :: sed(:, :), rflux(:), Av_grid(:), losvd_feature(:)
  real, allocatable           :: models(:, :, :), cached_models(:, :)
  real(8), allocatable        :: fluxo(:), sigma(:), chi_list(:, :)
  real(8), allocatable        :: coeff1d(:, :), coeff2d(:, :), coeff3d(:, :)
  real                        :: t1d(1), t2d(2), t3d(3), m1d(1), m2d(2), m3d(3)
  real                        :: chi(3), z1d(1), z2d(2), z3d(3)
  real                        :: mass_mod(3), t_m_mod(3), t_fr_mod(3), z_m_mod(3), z_fr_mod(3)
  real                        :: Av_mod(3), losvd_mod, losvd_wl(2), passband_weff
  real                        :: z, losvd_range(2), age_range(2), Av_range(2), basis_ages(12), Zsun
  character(:), allocatable   :: cl, dynbas_path, models_path, include, ignore, models_list
  character(256)              :: record, input_file, output_file, Av_strings(3), age_strings(3)
  character(256)              :: losvd_strings(4), passband_name*64
  character(100), allocatable :: long_opts(:), options(:), values(:), arguments(:), fname(:)*64
  character(256), allocatable :: kwords(:), files_ised(:), sel_files_ised(:)
  character(:), allocatable   :: dir_out, ID, ext_id
  logical, allocatable        :: mask(:)
  logical                     :: verbose, given_losvd, use_photon

! INITIALIZE VARIABLES WITH DEFAULT VALUES ---------------------------------------------------------
  CALL GET_PATH("DYNBAS", path=dynbas_path)

  allocate(ext_id, source="CCM")
  allocate(dir_out, source="./")
  allocate(ID, source="")
  allocate(include, source="")
  allocate(ignore, source="")
  allocate(cl, source="")

  passband     = 122
  nbas         = size(basis_ages)
  verbose      = .false.
  given_losvd  = .false.
  use_photon   = .false.
  z            = 0.0
  losvd_wl     = [3900.0, 4000.0]
  losvd_range  = [40.0, 450.0]
  age_range    = [1.05e6, 14.00e9]
  step         = 1
  Zsun         = 0.02
  Av_grid      = [0.000000000000,   1.72102769e-03,   2.96193630e-03,&
                  5.09757437e-03,   8.77306662e-03,   1.50986905e-02,&
                  2.59852645e-02,   4.47213595e-02,   7.69666979e-02,&
                  1.32461818e-01,   2.27970456e-01,   3.92343467e-01,&
                  6.75233969e-01,   1.16209635e+00,   1.500000000000]
  basis_ages   = [1.78e5, 4.17e6, 8.71e6, 1.91e7, 3.60e7, 9.05e7,&
                  3.60e8, 1.28e9, 2.75e9, 4.75e9, 9.25e9, 13.75e9]
  long_opts    = ["models=     ",&
                  "models!=    ",&
                  "help        ",&
                  "photon-count",&
                  "ID=         ",&
                  "output-dir= ",&
                  "passbands=  ",&
                  "LOSVD-wl=   ",&
                  "LOSVD-grid= ",&
                  "ext-curve=  ",&
                  "Av-grid=    ",&
                  "age-grid=   ",&
                  "L-passband= ",&
                  "Zsun=       "]
  kwords       = [character(len(kwords))::"redshift", "LOSVD", "Av"]
! --------------------------------------------------------------------------------------------------

! PARSE OPTIONS & ARGUMENTS ------------------------------------------------------------------------
  CALL PARSE_ARGS("hivz:", long_opts, options, values, arguments)

  if(.not.(allocated(options) .and. allocated(values) .and. allocated(arguments))) CALL HELP

  do m=1, size(options)
    select case(options(m))
    case("v")
      verbose = .true.
    case("L-passband")
      read(values(m), *) passband
    case("models")
      include = trim(values(m))
      n = len(include)
      if(include(1:1)/="*") include = "*"//include
      if(include(n:n)/="*") include = include//"*"
    case("models!")
      ignore = trim(values(m))
      n = len(ignore)
      if(ignore(1:1)/="*") ignore = "*"//ignore
      if(ignore(n:n)/="*") ignore = ignore//"*"
    case("output-dir")
      n       = len_trim(values(m))
      dir_out = values(m)(:n)

      if(dir_out(n:n) /= "/") dir_out = dir_out//"/"
    case("photon-count")
      use_photon = .true.
    case("ID")
      n  = len_trim(values(m))
      ID = values(m)(:n)

      if(ID(n:n)/= "_") ID = ID//"_"
    case("LOSVD-wl")
      n = STRING_COUNT(values(m), ",") + 1
      if(n==2) then
        losvd_strings(:2) = SPLIT(values(m), ",")
        read(losvd_strings(1), *) losvd_wl(1)
        read(losvd_strings(2), *) losvd_wl(2)
      end if
    case("LOSVD-grid")
      n = STRING_COUNT(values(m), ",") + 1
      if(n==2) then
        losvd_strings(:2) = SPLIT(values(m), ",")
        read(losvd_strings(1), *) losvd_range(1)
        read(losvd_strings(2), *) losvd_range(2)
      else if(n==4) then
        losvd_strings = SPLIT(values(m), ",")
        read(losvd_strings(1), *) losvd_range(1)
        read(losvd_strings(2), *) losvd_range(2)
        read(losvd_strings(3), *) losvd_wl(1)
        read(losvd_strings(4), *) losvd_wl(2)
      end if
    case("age-grid")
      n = STRING_COUNT(values(m), ",") + 1
      if(n==2) then
        age_strings(:2) = SPLIT(values(m), ",")
        read(age_strings(1), *) age_range(1)
        read(age_strings(2), *) age_range(2)
      else if(n==3) then
        age_strings = SPLIT(values(m), ",")
        read(age_strings(1), *) age_range(1)
        read(age_strings(2), *) age_range(2)
        step = EVAL(age_strings(3))
      end if
    case("ext-curve")
      n = len_trim(values(m))
      ext_id = values(m)(:n)
    case("Av-grid")
      n = STRING_COUNT(values(m), ",") + 1
      if(n==2) then
        Av_strings(:2) = SPLIT(values(m), ",")
        read(Av_strings(1), *) Av_range(1)
        read(Av_strings(2), *) Av_range(2)
        Av_grid = [(Av_range(1) + (i-1)*(Av_range(2)-Av_range(1))/9, i=1, 10)]
      else if(n==3) then
        Av_strings = SPLIT(values(m), ",")
        read(Av_strings(1), *) Av_range(1)
        read(Av_strings(2), *) Av_range(2)
        navs = EVAL(Av_strings(3))
        Av_grid = [(Av_range(1) + (i-1)*(Av_range(2)-Av_range(1))/(navs-1), i=1, navs)]
      end if
    case("Zsun")
      read(values(m), *) Zsun
    case("passbands")
      n = index(values(m), ",...,")
      if(n==0) then
        CALL APPEND(filter_ids, EVAL(SPLIT(values(m), ",")))
      else
        do k=1, STRING_COUNT(values(m), ",...,")
          i = index(values(m)(:n-1), ",", back=.true.)
          j = index(values(m)(n+5:), ",")
          if(i+j==0) then
            CALL APPEND(filter_ids, LISTING(values(m)))
          else if(i==0) then
            CALL APPEND(filter_ids, LISTING(values(m)(:j+n+3)))
          else if(j==0) then
            CALL APPEND(filter_ids, EVAL(SPLIT(values(m)(:i-1), ",")))
            CALL APPEND(filter_ids, LISTING(values(m)(i+1:)))
          else
            CALL APPEND(filter_ids, EVAL(SPLIT(values(m)(:i-1), ",")))
            CALL APPEND(filter_ids, LISTING(values(m)(i+1:j+n+3)))
          end if
          values(m) = values(m)(j+n+5:)
          n = index(values(m), ",...,")
        end do
        if(EVAL(values(m))/=filter_ids(size(filter_ids))) CALL APPEND(filter_ids, EVAL(SPLIT(values(m), ",")))
      end if

      nwef = size(filter_ids)
      allocate(fname(nwef), weff(nwef))

      do i=1, nwef
        CALL FILTER_INFO(filter_ids(i), mean_wl=weff(i), filter_name=fname(i))
      end do
    case("h", "help")
      CALL HELP
    end select
  end do

  if(size(arguments)/=1) STOP "DynBaS: number of arguments must be 1. Use --help option to see the documentation."

  input_file = arguments(1)
  CALL INQUIRE_FILE(input_file, nhdr, nwlo, ncol, status=iostatus)
  CALL IOERR(iostatus, input_file)

  unit = AVAILABLE_UNIT()
  open(unit, file=input_file, action="read")

  if(nhdr>0) then
    do i=1, nhdr
      read(unit, "(A)") record

      do j=1, size(kwords)
        if(index(record, trim(kwords(j)))==0) cycle

        select case(j)
        case(1)
          read(record(index(record, "=")+1:), *) z
        case(2)
          read(record(index(record, "=")+1:), *) losvd_mod
          given_losvd = .true.
        case(3)
          read(record(index(record, "=")+1:), *) Av_mod(1)
          Av_mod  = Av_mod(1)
          Av_grid = [Av_mod(1)]
        end select
        exit
      end do
    end do
  end if

  allocate(sed(nwlo, ncol))
  do i=1, nwlo
    read(unit, *) sed(i, :)
  end do
  close(unit)

  allocate(wlo(nwlo), fluxo(nwlo), sigma(nwlo))
  wlo   = sed(:, 1)
  fluxo = sed(:, 2)
  sigma = sed(:, 3)
  where(fluxo<=0.0) sigma = huge(1.0)

  if(.not. FILTER_EXIST(passband)) then
    write(6, *) "STOP DynBaS: selected passband for luminosity-weighted age/metallicity is not defined in $DB_FILTERS."
    stop
  else
    CALL FILTER_INFO(passband, mean_wl=passband_weff, filter_name=passband_name)
  end if
! --------------------------------------------------------------------------------------------------

! BUILD SSP MODEL MATRICES -------------------------------------------------------------------------
  CALL GET_PATH("DB_MODELS", path=models_path)
  CALL INQUIRE_FILE(models_path//"models_list.log", nrow=nmod, status=iostatus)
  if(iostatus/=0) then
    CALL SYSTEM("ls "//models_path//"*.ised|xargs -n1 basename > "//models_path//"models_list.log")
    CALL INQUIRE_FILE(models_path//"models_list.log", nrow=nmod, status=iostatus)
  end if
  if(nmod==0) then
    write(6, *) "STOP DynBaS: no .ised file was found in '"//models_path//"'. Please change the value of $DB_MODELS."
    stop
  end if
  unit = AVAILABLE_UNIT()
  open(unit, file=models_path//"models_list.log", action="read")

  allocate(files_ised(nmod), mask(nmod))
  mask = .false.
  do i=1, nmod
    read(unit, *) files_ised(i)
  end do
  close(unit)

  if(include=="" .and. ignore=="") then
    nmet = nmod
    allocate(sel_files_ised(nmet))
    sel_files_ised = files_ised
    mask = .true.
  else if(ignore=="") then
    cl = "ls "//models_path//include//"|xargs -n1 basename|grep .ised > "
  else if(include=="") then
    cl = "ls -I"//ignore//" "//models_path//"|xargs -n1 basename|grep .ised > "
  else
    cl = "ls -I"//ignore//" "//models_path//include//"|xargs -n1 basename|grep .ised > "
  endif

  if(cl/="") then
    allocate(models_list, source=models_path//"models_"//STRING(TIME())//".temp")
    CALL SYSTEM(cl//models_list, status=clstatus)
    if(clstatus/=0) then
      write(6, *) "STOP DynBaS: the include ", include, " does not match any .ised file."
      stop
    end if

    CALL INQUIRE_FILE(models_list, nrow=nmet, status=iostatus)
    CALL IOERR(iostatus, models_list)

    unit = AVAILABLE_UNIT()
    open(unit, file=models_list, action="read")

    allocate(sel_files_ised(nmet), iz(nmet))
    do i=1, nmet
      read(unit, "(A)") sel_files_ised(i)
      CALL STRING_LOC(files_ised, sel_files_ised(i), iz(i))
      mask(iz(i)) = .true.
    end do
    close(unit, status="delete")
  end if

  CALL READ_SED(models_path//sel_files_ised, wlength, mods, ages, Zs_, .true., .false.)

  if(z>0.0) then
    CALL INQUIRE_FILE(dynbas_path//"/data/universe_age.dat", nhdr, nrsh, ncol, iostatus)
    CALL IOERR(iostatus, dynbas_path//"/data/universe_age.dat")

    allocate(redshifts(nrsh), uni_ages(nrsh))

    open(unit, file=dynbas_path//"/data/universe_age.dat", action="read", iostat=iostatus)

    do i=1, nhdr
      read(unit, *)
    end do

    do i=1, nrsh
      read(unit, *) redshifts(i), uni_ages(i)
    end do
    close(unit)

    age_range(2) = minval([age_range(2), LINEAR_INTER(redshifts, uni_ages*1e9, [z])], dim=1)
  end if

  ages = ages(2:)
  mods = mods(:, 2:, :)

  m = minloc(abs(ages-age_range(1)), dim=1)
  n = minloc(abs(ages-age_range(2)), dim=1)

  nage = size(ages(m:n:step))
  nwlm = size(wlength)
  nmod = nage * nmet
  navs = size(Av_grid)

  allocate(fssps(nwlm, nmod), t(nmod), Zs(nmod))
  t     = reshape(spread(ages(m:n:step), 2, nmet), shape=[nmod])
  Zs    = reshape(spread(Zs_, 1, nage), shape=[nmod])
  fssps = reshape(mods(:, m:n:step, :), shape(fssps))

  allocate(chi_list(3, navs))
  allocate(coeff1d(1, navs), coeff2d(2, navs), coeff3d(3, navs))
  allocate(gen1d(1, navs), gen2d(2, navs), gen3d(3, navs))
  allocate(dyn_flux(nwlo, 3))
! --------------------------------------------------------------------------------------------------

! DISPLAY SED FITTING CONFIGURATION ----------------------------------------------------------------
  if(verbose) then
    write(6, *)
    write(6, *) "Available ISED models:"
    write(6, *) "---------------------"
    do i=1, size(files_ised)
    if(any(trim(files_ised(i))==sel_files_ised)) then
    write(6,             "(2A)") "[x]", trim(files_ised(i))
    else
    write(6,             "(2A)") "[ ]", trim(files_ised(i))
    end if
    end do
    write(6, *)
    write(6, *) "SED fitting parameters:"
    write(6, *) "----------------------"
    write(6,              "(A)") " * Model grid dimensions:"
    write(6,           "(A,I7)") "     > SSPs               = ", nage*nmet
    write(6,           "(A,I7)") "     > metallicities      = ", nmet
    write(6,           "(A,I7)") "     > ages               = ", nage
    write(6,           "(A,I7)") "     > ages step          = ", step
    write(6,           "(A,I7)") "     > extinctions        = ", navs
    write(6,           "(A,I7)") "     > wavelengths        = ", nwlo
    write(6,           "(A,I7)") "     > eff. wavelengths   = ", count(fluxo>0.0)
    write(6,              "(A)") " * Physical properties:"
    write(6,           "(A,A7)") "     > ext. curve ID      = ", ext_id
    write(6,        "(A,I7,3A)") "     > L-weight passband  = ", passband, " (", trim(adjustl(passband_name)), ")"
    write(6,         "(A,F7.3)") "     > L-weight wl [nm]   = ", passband_weff/10.0
    write(6,         "(A,F7.3)") "     > z (redshift)       = ", z
    write(6,         "(A,F7.3)") "     > Zsun               = ", Zsun
    write(6, "(A,F7.3,A3,F7.3)") "     > age range [log/yr] = ", log10(t(1)),      "--", log10(t(nage))
    write(6, "(A,F7.3,A3,F7.3)") "     > Z   range [Zsun]   = ", minval(Zs_)/Zsun, "--", maxval(Zs_)/Zsun
    write(6, "(A,F7.3,A3,F7.3)") "     > Av  range [mag]    = ", minval(Av_grid),  "--", maxval(Av_grid)

    if(allocated(filter_ids)) then
    write(6,              "(A)") " * Passbands:"
    do i=1, nwef
    write(6,             "(2A)") "     > ", trim(fname(i))
    end do
    end if
    write(6, *)
  end if
! --------------------------------------------------------------------------------------------------

! START SED FITTING PROCEDURE ----------------------------------------------------------------------
  write(6, "(A)", advance="no") "Computing model grid..."
  allocate(models(nwlo, nmod, navs))
  if(allocated(filter_ids)) then
    models(:, :, 1) = EFFECTIVE_FLUX(filter_ids, wlength, fssps, z, photon_count=use_photon)
  else
    models(:, :, 1) = LINEAR_INTER(wlength*(1+z), fssps, wlo)

    if(.not.given_losvd .or. losvd_mod<0.0) then
      allocate(ibasis(nbas*nmet), basis(nwlo, nbas*nmet))
      basis = LINEAR_INTER(wlength*(1+z), reshape(mods, [nwlm, size(ages)*nmet]), wlo)
      do j = 1, nbas
        CALL NEAREST_ELEMENT(ages, basis_ages(j), ibasis(j))

        forall(i=1:nmet-1) ibasis(j+i*nbas) = ibasis(j) + i*size(ages)
      end do

      basis = basis(:, ibasis)

      allocate(fluxo_s(nwlo), sigma_s(nwlo))
      fluxo_s = fluxo
      sigma_s = sigma

      if(losvd_wl(1)<wlo(1) .or. losvd_wl(2)>wlo(nwlo)) then
        write(6, *) "STOP DynBaS: the chosen feature to compute the LOSVD is not spanned by the input SED."
        write(6, *) "             Use --losvd-feature to change the deafult wavelength range."
        stop
      end if
      losvd_feature = pack(fluxo_s, wlo>=losvd_wl(1).and.wlo<=losvd_wl(2))
      if(count(losvd_feature>0.0)/size(losvd_feature)<0.6) then
        write(6, *) "WARNING DynBaS: the chosen feature to compute the LOSVD has more than 40% fluxes masked out."
      else if(count(losvd_feature>0.0)/size(losvd_feature)<0.4) then
        STOP "DynBaS: the chosen feature to compute the LOSVD has more than 60% fluxes masked out."
      end if

      CALL VELOCITY_DISP(wlo, basis, fluxo_s, sigma_s, losvd_wl, losvd_range, losvd_mod)
    end if

    if(losvd_mod>losvd_range(1)) then
      CALL APPLY_VD(wlo, models(:, :, 1), losvd_mod)
    else
      losvd_mod = 0.0
    end if
  end if

  cached_models = models(:, :, 1)

  !$OMP PARALLEL DO
  do m=1, navs
    models(:, :, m) = cached_models * spread(EXT(ext_id, wlo/(z+1), Av_grid(m)), 2, nmod)
  end do
  !$OMP END PARALLEL DO

  deallocate(cached_models)

  write(6, *) "         done."

  write(6, "(A)", advance="no") "Computing minimum chi square..."

  !$OMP PARALLEL DO
  do m=1, navs
    CALL MIN_CHI_SQUARE(fluxo/sigma, dble(models(:, :, m))/spread(sigma, 2, nmod), &
                        coeff1d(:, m), coeff2d(:, m), coeff3d(:, m), &
                        gen1d(:, m), gen2d(:, m), gen3d(:, m), chi_list(:, m))
  end do
  !$OMP END PARALLEL DO

  write(6, *) " done."
  write(6, *)
! --------------------------------------------------------------------------------------------------

! START POST-PROCESSING & WRITE OUTPUT FILE --------------------------------------------------------
  ind = minloc(chi_list, dim=2)

  Av_mod(1) = Av_grid(ind(1))
  Av_mod(2) = Av_grid(ind(2))
  Av_mod(3) = Av_grid(ind(3))

  m1d  = coeff1d(:, ind(1))
  m2d  = coeff2d(:, ind(2))
  m3d  = coeff3d(:, ind(3))
  it1d = gen1d(:, ind(1))
  it2d = gen2d(:, ind(2))
  it3d = gen3d(:, ind(3))

  chi  = [(chi_list(j, ind(j)), j=1, 3)]

  if(any(it1d<1) .or. any(it1d>nmod)) then
    it1d           = 0
    dyn_flux(:, 1) = 0.0
    mass_mod(1)    = 0.0
    t_m_mod(1)     = 0.0
    t_fr_mod(1)    = 0.0
    z_m_mod(1)     = 0.0
    z_fr_mod(1)    = 0.0
    Av_mod(1)      = 0.0
    chi(1)         = 0.0
  else
    t1d = t(it1d(1))
    z1d = Zs(it1d(1))

    dyn_flux(:, 1) = m1d(1) * models(:, it1d(1), ind(1))
    mass_mod(1)    = m1d(1)
    t_m_mod(1)     = log10(t1d(1))
    t_fr_mod(1)    = log10(t1d(1))
    z_m_mod(1)     = log10(Zs(it1d(1)) / Zsun)
    z_fr_mod(1)    = log10(Zs(it1d(1)) / Zsun)
  end if
  if(any(it2d<1) .or. any(it2d>nmod)) then
    it2d           = 0
    dyn_flux(:, 2) = 0.0
    mass_mod(2)    = 0.0
    t_m_mod(2)     = 0.0
    t_fr_mod(2)    = 0.0
    z_m_mod(2)     = 0.0
    z_fr_mod(2)    = 0.0
    Av_mod(2)      = 0.0
    chi(2)         = 0.0
  else
    t2d = [(t(it2d(j)),  j=1, 2)]
    z2d = [(Zs(it2d(j)), j=1, 2)]

    CALL MONOCHROMATIC_SEDS(wlength, fssps*spread(EXT(ext_id, wlength, Av_mod(2)), 2, nmod), passband, 0.0, it2d, rflux)

    dyn_flux(:, 2) = m2d(1) * models(:, it2d(1), ind(2)) + m2d(2) * models(:, it2d(2), ind(2))
    mass_mod(2)    = sum(m2d)
    t_m_mod(2)     = MEAN(log10(t2d), m2d)
    t_fr_mod(2)    = MEAN(log10(t2d), m2d*rflux)
    z_m_mod(2)     = MEAN(log10(z2d/Zsun), m2d)
    z_fr_mod(2)    = MEAN(log10(z2d/Zsun), m2d*rflux)
  end if
  if(any(it3d<1) .or. any(it3d>nmod)) then
    it3d           = 0
    dyn_flux(:, 3) = 0.0
    mass_mod(3)    = 0.0
    t_m_mod(3)     = 0.0
    t_fr_mod(3)    = 0.0
    z_m_mod(3)     = 0.0
    z_fr_mod(3)    = 0.0
    Av_mod(3)      = 0.0
    chi(3)         = 0.0
  else
    t3d = [(t(it3d(j)),  j=1, 3)]
    z3d = [(Zs(it3d(j)), j=1, 3)]

    CALL MONOCHROMATIC_SEDS(wlength, fssps*spread(EXT(ext_id, wlength, Av_mod(3)), 2, nmod), passband, 0.0, it3d, rflux)

    dyn_flux(:, 3) = m3d(1) * models(:, it3d(1), ind(3)) + m3d(2) * models(:, it3d(2), ind(3)) + m3d(3) * models(:, it3d(3), ind(3))
    mass_mod(3)    = sum(m3d)
    t_m_mod(3)     = MEAN(log10(t3d), m3d)
    t_fr_mod(3)    = MEAN(log10(t3d), m3d*rflux)
    z_m_mod(3)     = MEAN(log10(z3d/Zsun), m3d)
    z_fr_mod(3)    = MEAN(log10(z3d/Zsun), m3d*rflux)
  end if

  where(fluxo<=0.0) sigma = 0.0

  m = index(input_file, "/", back=.true.)
  n = index(input_file, ".", back=.true.)
  if(n<m) n = len_trim(input_file) + 1

  if(ID/="") then
    output_file = dir_out//"dynbasfit_"//ID//input_file(m+1:n-1)//".log"
  else
    output_file = dir_out//"dynbasfit_"//input_file(m+1:n-1)//".log"
  end if

  unit = AVAILABLE_UNIT()
  open(unit, file=trim(output_file), action="write")

  write(unit,                                 "(2A)") "# SED file          = ", trim(input_file)
  write(unit,                                  "(A)") "#"
  write(unit,         "(A,"//STRING(nage)//"EN13.2)") "# t grid            = ", t(:nage)
  write(unit, "(A,"//STRING(nmet)//"G10.3,A,F5.3,A)") "# Z/Zsun grid       = ", Zs_ / Zsun, "(Zsun=", Zsun, ")"
  write(unit,       "(A,"//STRING(navs)//"G10.3,3A)") "# Av grid           = ", Av_grid, "(", ext_id, ")"
  write(unit,                                  "(A)") "#"
  write(unit,                              "(A,1I7)") "# generator 1d      = ", it1d
  write(unit,                              "(A,2I7)") "# generator 2d      = ", it2d
  write(unit,                              "(A,3I7)") "# generator 3d      = ", it3d
  write(unit,                          "(A,1EN13.2)") "# t 1d              = ", t1d
  write(unit,                          "(A,2EN13.2)") "# t 2d              = ", t2d
  write(unit,                          "(A,3EN13.2)") "# t 3d              = ", t3d
  write(unit,                           "(A,1G10.3)") "# Z/Zsun 1d         = ", z1d / Zsun
  write(unit,                           "(A,2G10.3)") "# Z/Zsun 2d         = ", z2d / Zsun
  write(unit,                           "(A,3G10.3)") "# Z/Zsun 3d         = ", z3d / Zsun
  write(unit,                           "(A,1E10.2)") "# coefficients 1d   = ", m1d
  write(unit,                           "(A,2E10.2)") "# coefficients 2d   = ", m2d
  write(unit,                           "(A,3E10.2)") "# coefficients 3d   = ", m3d
  write(unit,                                  "(A)") "#"
  write(unit,                            "(A,F10.5)") "# redshift          = ", z
  if(allocated(filter_ids)) then
  write(unit,                                 "(2A)") "# LOSVD             = ", "None"
  else
  write(unit,                             "(A,F8.2)") "# LOSVD             = ", losvd_mod
  end if
  write(unit,                                 "(4A)") "# prop.\dimension   = ", "             1d",&
                                                                                "             2d",&
                                                                                "             3d"
  write(unit,                           "(A,3E15.6)") "# sum(coefficients) = ", mass_mod
  write(unit,                           "(A,3F15.6)") "# <log t>_M         = ", t_m_mod
  write(unit,                           "(A,3F15.6)") "# <log t>_Lr        = ", t_fr_mod
  write(unit,                           "(A,3F15.6)") "# <log Z/Zo>_M      = ", z_m_mod
  write(unit,                           "(A,3F15.6)") "# <log Z/Zo>_Lr     = ", z_fr_mod
  write(unit,                           "(A,3F15.6)") "# Av[mag]           = ", Av_mod
  write(unit,                           "(A,3F15.6)") "# chi square        = ", chi
  write(unit,                                  "(A)") "#"
  write(unit,                                 "(6A)") "#    wavelength", "       obs_flux",&
                                                      "          sigma", "        1d_flux",&
                                                      "        2d_flux", "        3d_flux"
  do i = 1, nwlo
  write(unit, "(F15.2,5E15.6)") wlo(i), fluxo(i), sigma(i), dyn_flux(i, 1), dyn_flux(i, 2), dyn_flux(i, 3)
  end do
! --------------------------------------------------------------------------------------------------

END PROGRAM
