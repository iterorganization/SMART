include 'mathematical_constants.f90'

module mod_smart

contains

   subroutine smart(eq_in, cp_in, pellets_in, cp_out, codeparam, &
                    error_flag, error_message)

      ! ---------------------------------------
      ! PELLET ABLATION MODEL FROM ASTRA: SMART
      ! IDS INPUT: EQUILIBRIUM, CORE_PROFILES
      ! IDS OUTPUT: CORE_PROFILES
      ! ---------------------------------------

      use ids_schemas,  only : ids_equilibrium, ids_core_profiles, ids_pellets
      use ids_schemas,  only : ids_parameters_input, ids_is_valid
      use ids_routines, only : ids_copy
      use mod_codeparam_smart
      use mathematical_constants, only : M_PI
      implicit none

      ! Arguments
      type(ids_equilibrium)        :: eq_in
      type(ids_core_profiles)      :: cp_in, cp_out
      type(ids_pellets)            :: pellets_in
      type(ids_parameters_input)   :: codeparam
      integer, intent(out)         :: error_flag
      character(len=:), pointer, intent(out) :: error_message

      ! Local derived type
      type(type_smart_data) :: smart_in

      ! Flags
      logical :: from_pellets_ids = .false.

      ! Integer scalars
      integer :: i, j, i_time, j_time, n_xcp, n_xeq, n_ion, nrd, NA1, jtime, jend
      integer :: iH, iD, iT, na, nb1, Khydr, JABS, istep, NABEG, NEBEG, NE1
      integer :: jheat, jdens, jpsi, jpel, jprint, nbnd1, nbnd2, IMIX, jmix
      integer :: nkey4control

      ! Double precision scalars
      double precision :: denA2D, temA2D, presA2D, cuA2D
      double precision :: HRO, ROC, BTOR, SHIFT, ABC, RTOR, ALFA
      double precision :: YDABL, YDDEP, ai, zi, Dtpel
      double precision :: NEB, TEB, TIB, F0B, F1B, F2B, F3B, F4B, F5B, F6B, F7B, F8B, F9B
      double precision :: F01B, F02B, F03B, GN2E, GN2I, RHOEC, RHODR
      double precision :: QECR, YEFFec, y, y1, y2, y3, y4, pdt, QNB, IPL
      double precision :: PDD1, PDD2, PDD3, TAUSp
      double precision :: elon, trian, sign_psi
      double precision :: TAU, dtau, TIME, TIMBEG
      double precision, save :: TIMPEL = 0.d0

      ! Arrays
      integer, allocatable :: ispec(:), key4control(:)
      double precision, allocatable :: ne(:), ni(:), Te(:), Ti(:), TN(:)
      double precision, allocatable :: F0(:), F1(:), F2(:), F3(:), F4(:), F5(:)
      double precision, allocatable :: F6(:), F7(:), F8(:), F9(:), FP(:)
      double precision, allocatable :: neo(:), nio(:), Teo(:), Tio(:)
      double precision, allocatable :: F0o(:), F1o(:), F2o(:), F3o(:), F4o(:), F5o(:)
      double precision, allocatable :: F6o(:), F7o(:), F8o(:), F9o(:), FPo(:)
      double precision, allocatable :: neX(:), niX(:), TeX(:), TiX(:)
      double precision, allocatable :: F0x(:), F1x(:), F2x(:), F3x(:), F4x(:), F5x(:)
      double precision, allocatable :: F6x(:), F7x(:), F8x(:), F9x(:), NN(:)
      double precision, allocatable :: cu(:), cutor(:), cd(:), cubs(:), UPL(:), ULON(:), EZ(:), FV(:)
      double precision, allocatable :: ZEF(:), AMAIN(:), Z2NdA(:), ZMAIN(:), SQEPS(:)
      double precision, allocatable :: F1fast(:), F2fast(:), F3fast(:), F4fast(:), F5fast(:)
      double precision, allocatable :: NHYDR(:), NDEUT(:), NTRIT(:), NALF(:), NHE3(:)
      double precision, allocatable :: PBLON(:), PBPER(:), PFAST(:), VTOR(:)
      double precision, allocatable :: SN(:), SNN(:), SNTOT(:), QN(:), GN(:), GNX(:), QE(:), QI(:)
      double precision, allocatable :: SF0(:), SFF0(:), SF0TOT(:), QF0(:), GF0(:), GF0X(:)
      double precision, allocatable :: SF1(:), SFF1(:), SF1TOT(:), QF1(:), GF1(:), GF1X(:)
      double precision, allocatable :: SF2(:), SFF2(:), SF2TOT(:), QF2(:), GF2(:), GF2X(:)
      double precision, allocatable :: SF3(:), SFF3(:), SF3TOT(:), QF3(:), GF3(:), GF3X(:)
      double precision, allocatable :: SF4(:), SFF4(:), SF4TOT(:), QF4(:), GF4(:), GF4X(:)
      double precision, allocatable :: SF5(:), SFF5(:), SF5TOT(:), QF5(:), GF5(:), GF5X(:)
      double precision, allocatable :: SF6(:), SFF6(:), SF6TOT(:), QF6(:), GF6(:), GF6X(:)
      double precision, allocatable :: SF7(:), SFF7(:), SF7TOT(:), QF7(:), GF7(:), GF7X(:)
      double precision, allocatable :: SF8(:), SFF8(:), SF8TOT(:), QF8(:), GF8(:), GF8X(:)
      double precision, allocatable :: SF9(:), SFF9(:), SF9TOT(:), QF9(:), GF9(:), GF9X(:)
      double precision, allocatable :: PE(:), PET(:), PETOT(:), PI(:), PIT(:), PITOT(:), PEI(:)
      double precision, allocatable :: PEECR(:), CUECR(:), YPELSRS(:)
      double precision, allocatable :: PEFUS(:), PIFUS(:), PEAUX(:), PIAUX(:)
      double precision, allocatable :: PEN(:), PIN(:), PJOUL(:), PRAD(:)
      double precision, allocatable :: Sn14(:), Sn245(:), SCUBM(:), PEICR(:), PIICR(:), PEBM(:), PIBM(:)
      double precision, allocatable :: DF0(:), VF0(:), DF1(:), VF1(:), DF2(:), VF2(:), DF3(:), VF3(:)
      double precision, allocatable :: DF4(:), VF4(:), DF5(:), VF5(:), DF6(:), VF6(:), DF7(:), VF7(:)
      double precision, allocatable :: DF8(:), VF8(:), DF9(:), VF9(:), DN(:), CN(:), HE(:), XI(:)
      double precision, allocatable :: cc(:), DSI(:), DSE(:), DSN(:)
      double precision, allocatable :: ametr(:), shif(:), mu(:), VOL(:), vr(:), VOLo(:), vro(:)
      double precision, allocatable :: IPOL(:), G11(:), G33(:), SLAT(:)
      double precision, allocatable :: BMINT(:), BMAXT(:), BDB0(:), BDB02(:), B0DB2(:)
      double precision, allocatable :: FOFB(:), G22(:), EQFF(:), EQPF(:)
      double precision, allocatable :: ametre(:), shife(:), VOLe(:), XCP(:), XEQ(:), RHO(:), FPe(:)
      double precision, allocatable :: IPOLe(:), G11e(:), G33e(:), SLATe(:)
      double precision, allocatable :: BMINTe(:), BMAXTe(:), BDB0e(:), BDB02e(:), B0DB2e(:), FOFBe(:), G22e(:)

      character(len=14) :: ARRNAME(30), VARNAME(30)

      double precision, external :: VINTa

      data ARRNAME / ' x ',' Te ',' Ti ',' ne ',' ni ', &
                     ' nHth ',' nHf ',' nDth ',' nTth ',' nTf ',' n4Heth ',' n4Hef ', ' n3Heth ',' n3Hef ', &
                     ' U ',' q ',' Fp ',' Jtot',' Jbs ',' Jec ',' Zeff ',' Zmain ',' Amain ', &
                     ' Pe ',' Pi ',' Pefus ',' Pifus ',' Poh ',' Pec ',' Pei ' /

      data VARNAME / ' time ',' Te0 ',' Ti0 ',' ne0 ',' ni0 ',' Tea ',' Tia ',' nea ',' nia ', &
                     ' n00 ',' n0a ',' <ne> ', ' Ibs ',' Itot ',' Icd ',' U0 ',' Ua ',' Fp0 ',' Fp ', &
                     ' Pe ',' Pi ',' Pefus ',' Pifus ',' Poh ',' Pec ',' Pei ', &
                     ' q(0) ',' q(a) ',' <Shdt> ',' <Sn0> ' /

      ! Initialization of error flag
      error_flag = 0
      nullify(error_message)

      ! Unit conversions
      denA2D  = 1.d19
      temA2D  = 1.d3
      presA2D = 1.6d3
      cuA2D   = 1.d6

      ! Pellet IDS present?
      if (pellets_in%ids_properties%homogeneous_time >= 0) from_pellets_ids = .true.

      ! Validate input IDS
      if (ids_is_valid(eq_in%ids_properties%homogeneous_time) .and. &
          size(eq_in%time) > 0 .and. &
          ids_is_valid(cp_in%ids_properties%homogeneous_time) .and. &
          size(cp_in%time) > 0) then

         call assign_codeparam(codeparam%parameters_value, smart_in)

         ! Copy control keys
         if (allocated(key4control)) deallocate(key4control)
         if (allocated(smart_in%key4control) .and. size(smart_in%key4control) >= 15) then
            nkey4control = size(smart_in%key4control)
            allocate(key4control(nkey4control))
            key4control = smart_in%key4control
         else
            nkey4control = 15
            allocate(key4control(nkey4control))
            key4control = 0
            if (allocated(smart_in%key4control) .and. size(smart_in%key4control) > 0) &
               key4control(1:size(smart_in%key4control)) = smart_in%key4control
            write(*,*) 'Notice: Control keys missing or too short in XML, using defaults where needed'
         end if
         write(*,'(A,*(I2,1X))') 'key4control=', key4control(1:min(15,nkey4control))

         if (from_pellets_ids) then
            write(*,*) 'Input pellets IDS detected'
            smart_in%YAM  = pellets_in%time_slice(1)%pellet(1)%species(1)%a
            smart_in%YVP  = pellets_in%time_slice(1)%pellet(1)%velocity_initial * 1.e-3
            smart_in%YVOL = pellets_in%time_slice(1)%pellet(1)%shape%size(1)**2 * &
                            pellets_in%time_slice(1)%pellet(1)%shape%size(2) * M_PI * 1.e+9
         else
            if (smart_in%sw_stdout /= 0) write(*,*) 'Input pellets IDS NOT detected'
         end if

         call ids_copy(cp_in, cp_out)
      else
         error_flag = -1
         allocate(character(50) :: error_message)
         error_message = 'Error in SMART: input IDS not valid'
         return
      end if

      !======================
      ! Initialization
      !======================
      i_time = size(cp_in%time)
      TIMBEG = cp_in%time(i_time)
      if (smart_in%sw_stdout /= 0) write(*,'(/,A17, I5, F10.4)') 'i_time, TIMBEG = ', i_time, TIMBEG

      j_time = size(eq_in%time_slice)
      n_xcp  = size(cp_in%profiles_1d(i_time)%grid%rho_tor_norm)
      n_ion  = size(cp_in%profiles_1d(i_time)%ion)
      n_xeq  = size(eq_in%time_slice(j_time)%profiles_1d%psi)
      nrd    = 2 * max(n_xcp, n_xeq)

      ! Temporary block
      nbnd1 = size(eq_in%time_slice(j_time)%boundary%outline%r)
      write(*,*) 'nbnd1=', nbnd1
      open(20, file='IMAS.dat')
      write(20,*) i_time, nbnd1
      write(20,*) j_time, TIMBEG
      do j = 1, nbnd1
         write(20,*) eq_in%time_slice(j_time)%boundary%outline%r(j), &
                     eq_in%time_slice(j_time)%boundary%outline%z(j)
      end do
      elon  = eq_in%time_slice(j_time)%boundary%elongation
      trian = eq_in%time_slice(j_time)%boundary%triangularity
      close(20)

      ! Elimination of fake -1/2 point
      if (cp_in%profiles_1d(i_time)%grid%rho_tor_norm(1) <= 0.d0) then
         NABEG = 2; NA1 = n_xcp - 1
      else
         NABEG = 1; NA1 = n_xcp
      end if
      if (eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(1) <= 0.d0) then
         NEBEG = 2; NE1 = n_xeq - 1
      else
         NEBEG = 1; NE1 = n_xeq
      end if

      NB1 = NA1; NA = NA1 - 1

      allocate(ispec(n_ion + 3))

      ! Allocate large profile arrays
      allocate(ne(n_xcp), ni(n_xcp), Te(n_xcp), Ti(n_xcp), &
               neX(n_xcp), niX(n_xcp), TeX(n_xcp), TiX(n_xcp), &
               F1x(n_xcp), F2x(n_xcp), F3x(n_xcp), F4x(n_xcp), F5x(n_xcp), &
               F6x(n_xcp), F7x(n_xcp), F8x(n_xcp), F9x(n_xcp), &
               F1(n_xcp), F2(n_xcp), F3(n_xcp), F4(n_xcp), F5(n_xcp), &
               F6(n_xcp), F7(n_xcp), F8(n_xcp), F9(n_xcp), FP(n_xcp), FV(n_xcp), &
               neo(n_xcp), nio(n_xcp), Teo(n_xcp), Tio(n_xcp), &
               F1o(n_xcp), F2o(n_xcp), F3o(n_xcp), F4o(n_xcp), F5o(n_xcp), &
               F6o(n_xcp), F7o(n_xcp), F8o(n_xcp), F9o(n_xcp), FPo(n_xcp), VTOR(n_xcp), &
               F0(n_xcp), F0o(n_xcp), F0x(n_xcp), NN(n_xcp), TN(n_xcp), &
               cu(n_xcp), cutor(n_xcp), cd(n_xcp), cubs(n_xcp), &
               UPL(n_xcp), ULON(n_xcp), EZ(n_xcp), ZEF(n_xcp), AMAIN(n_xcp), Z2NdA(n_xcp), ZMAIN(n_xcp), &
               PBLON(n_xcp), PBPER(n_xcp), PFAST(n_xcp), &
               F1fast(n_xcp), F2fast(n_xcp), F3fast(n_xcp), F4fast(n_xcp), F5fast(n_xcp), &
               NHYDR(n_xcp), NDEUT(n_xcp), NTRIT(n_xcp), NALF(n_xcp), NHE3(n_xcp) )

      ! Zero initialize
      ne = 0.d0; ni = 0.d0; Te = 0.d0; Ti = 0.d0
      neX = 0.d0; niX = 0.d0; TeX = 0.d0; TiX = 0.d0
      F1x = 0.d0; F2x = 0.d0; F3x = 0.d0; F4x = 0.d0; F5x = 0.d0
      F6x = 0.d0; F7x = 0.d0; F8x = 0.d0; F9x = 0.d0
      F1  = 0.d0; F2  = 0.d0; F3  = 0.d0; F4  = 0.d0; F5  = 0.d0
      F6  = 0.d0; F7  = 0.d0; F8  = 0.d0; F9  = 0.d0
      FP  = 0.d0; FV  = 0.d0
      neo = 0.d0; nio = 0.d0; Teo = 0.d0; Tio = 0.d0
      F1o = 0.d0; F2o = 0.d0; F3o = 0.d0; F4o = 0.d0; F5o = 0.d0
      F6o = 0.d0; F7o = 0.d0; F8o = 0.d0; F9o = 0.d0; FPo = 0.d0
      F0  = 0.d0; F0o = 0.d0; F0x = 0.d0
      NN  = 0.d0; TN  = 0.d0
      VTOR = 0.d0; cu = 0.d0; cutor = 0.d0; cd = 0.d0; cubs = 0.d0
      UPL = 0.d0; ULON = 0.d0; EZ = 0.d0
      ZEF = 0.d0; AMAIN = 1.d0; Z2NdA = 1.d0; ZMAIN = 1.d0
      NHYDR = 0.d0; NDEUT = 0.d0; NALF = 0.d0; NHE3 = 0.d0
      F1fast = 0.d0; F2fast = 0.d0; F3fast = 0.d0; F4fast = 0.d0; F5fast = 0.d0

      allocate(ametr(n_xcp), shif(n_xcp), vr(n_xcp), vro(n_xcp), &
               mu(n_xcp), VOL(n_xcp), VOLo(n_xcp), XCP(n_xcp), &
               IPOL(n_xcp), G11(n_xcp), G33(n_xcp), SLAT(n_xcp), RHO(n_xcp), SQEPS(n_xcp), &
               BMINT(n_xcp), BMAXT(n_xcp), BDB0(n_xcp), BDB02(n_xcp), B0DB2(n_xcp), &
               FOFB(n_xcp), G22(n_xcp), EQFF(n_xcp), EQPF(n_xcp) )
      ametr = 0.d0; shif = 0.d0; vr = 0.d0; vro = 0.d0; mu = 0.d0
      VOL = 0.d0; VOLo = 0.d0; XCP = 0.d0; IPOL = 0.d0; G11 = 0.d0; G33 = 0.d0
      SLAT = 0.d0; BMINT = 0.d0; BMAXT = 0.d0; BDB0 = 0.d0; BDB02 = 0.d0; B0DB2 = 0.d0
      FOFB = 0.d0; G22 = 0.d0; EQFF = 0.d0; EQPF = 0.d0

      allocate(ametre(n_xeq), shife(n_xeq), VOLe(n_xeq), &
               FPe(n_xeq), XEQ(n_xeq), &
               IPOLe(n_xeq), G11e(n_xeq), G33e(n_xeq), SLATe(n_xeq), &
               BMINTe(n_xeq), BMAXTe(n_xeq), BDB0e(n_xeq), BDB02e(n_xeq), B0DB2e(n_xeq), &
               FOFBe(n_xeq), G22e(n_xeq) )
      ametre = 0.d0; shife = 0.d0; VOLe = 0.d0; FPe = 0.d0; XEQ = 0.d0; RHO = 0.d0
      IPOLe = 0.d0; G11e = 0.d0; G33e = 0.d0; SLATe = 0.d0; BMINTe = 0.d0; BMAXTe = 0.d0
      BDB0e = 0.d0; BDB02e = 0.d0; B0DB2e = 0.d0; FOFBe = 0.d0

      allocate(SN(n_xcp), SNN(n_xcp), SNTOT(n_xcp), QN(n_xcp), GN(n_xcp), GNX(n_xcp), &
               SF0(n_xcp), SFF0(n_xcp), SF0TOT(n_xcp), QF0(n_xcp), GF0(n_xcp), GF0X(n_xcp), &
               SF1(n_xcp), SFF1(n_xcp), SF1TOT(n_xcp), QF1(n_xcp), GF1(n_xcp), GF1X(n_xcp), &
               SF2(n_xcp), SFF2(n_xcp), SF2TOT(n_xcp), QF2(n_xcp), GF2(n_xcp), GF2X(n_xcp), &
               SF3(n_xcp), SFF3(n_xcp), SF3TOT(n_xcp), QF3(n_xcp), GF3(n_xcp), GF3X(n_xcp), &
               SF4(n_xcp), SFF4(n_xcp), SF4TOT(n_xcp), QF4(n_xcp), GF4(n_xcp), GF4X(n_xcp), &
               SF5(n_xcp), SFF5(n_xcp), SF5TOT(n_xcp), QF5(n_xcp), GF5(n_xcp), GF5X(n_xcp), &
               SF6(n_xcp), SFF6(n_xcp), SF6TOT(n_xcp), QF6(n_xcp), GF6(n_xcp), GF6X(n_xcp), &
               SF7(n_xcp), SFF7(n_xcp), SF7TOT(n_xcp), QF7(n_xcp), GF7(n_xcp), GF7X(n_xcp), &
               SF8(n_xcp), SFF8(n_xcp), SF8TOT(n_xcp), QF8(n_xcp), GF8(n_xcp), GF8X(n_xcp), &
               SF9(n_xcp), SFF9(n_xcp), SF9TOT(n_xcp), QF9(n_xcp), GF9(n_xcp), GF9X(n_xcp), &
               DF0(n_xcp), VF0(n_xcp), DF1(n_xcp), VF1(n_xcp), DF2(n_xcp), VF2(n_xcp), DF3(n_xcp), VF3(n_xcp), &
               DF4(n_xcp), VF4(n_xcp), DF5(n_xcp), VF5(n_xcp), DF6(n_xcp), VF6(n_xcp), DF7(n_xcp), VF7(n_xcp), &
               DF8(n_xcp), VF8(n_xcp), DF9(n_xcp), VF9(n_xcp), &
               PE(n_xcp), PET(n_xcp), PETOT(n_xcp), QE(n_xcp), QI(n_xcp), &
               PI(n_xcp), PIT(n_xcp), PITOT(n_xcp), PEI(n_xcp), &
               PEFUS(n_xcp), PIFUS(n_xcp), PEAUX(n_xcp), PIAUX(n_xcp), &
               PEN(n_xcp), PIN(n_xcp), PJOUL(n_xcp), PRAD(n_xcp), &
               Sn14(n_xcp), Sn245(n_xcp), SCUBM(n_xcp), &
               PEICR(n_xcp), PIICR(n_xcp), PEBM(n_xcp), PIBM(n_xcp), &
               YPELSRS(n_xcp) )

      IMIX = 0
      SN = 0.d0; SNN = 0.d0; SNTOT = 0.d0; QN = 0.d0; GN = 0.d0; GNX = 0.d0
      SF0 = 0.d0; SFF0 = 0.d0; SF0TOT = 0.d0; QF0 = 0.d0; GF0 = 0.d0; GF0X = 0.d0
      SF1 = 0.d0; SFF1 = 0.d0; SF1TOT = 0.d0; QF1 = 0.d0; GF1 = 0.d0; GF1X = 0.d0
      SF2 = 0.d0; SFF2 = 0.d0; SF2TOT = 0.d0; QF2 = 0.d0; GF2 = 0.d0; GF2X = 0.d0
      SF3 = 0.d0; SFF3 = 0.d0; SF3TOT = 0.d0; QF3 = 0.d0; GF3 = 0.d0; GF3X = 0.d0
      SF4 = 0.d0; SFF4 = 0.d0; SF4TOT = 0.d0; QF4 = 0.d0; GF4 = 0.d0; GF4X = 0.d0
      SF5 = 0.d0; SFF5 = 0.d0; SF5TOT = 0.d0; QF5 = 0.d0; GF5 = 0.d0; GF5X = 0.d0
      SF6 = 0.d0; SFF6 = 0.d0; SF6TOT = 0.d0; QF6 = 0.d0; GF6 = 0.d0; GF6X = 0.d0
      SF7 = 0.d0; SFF7 = 0.d0; SF7TOT = 0.d0; QF7 = 0.d0; GF7 = 0.d0; GF7X = 0.d0
      SF8 = 0.d0; SFF8 = 0.d0; SF8TOT = 0.d0; QF8 = 0.d0; GF8 = 0.d0; GF8X = 0.d0
      SF9 = 0.d0; SFF9 = 0.d0; SF9TOT = 0.d0; QF9 = 0.d0; GF9 = 0.d0; GF9X = 0.d0
      DF0 = 0.d0; VF0 = 0.d0; DF1 = 0.d0; VF1 = 0.d0
      DF2 = 0.d0; VF2 = 0.d0; DF3 = 0.d0; VF3 = 0.d0
      DF4 = 0.d0; VF4 = 0.d0; DF5 = 0.d0; VF5 = 0.d0
      DF6 = 0.d0; VF6 = 0.d0; DF7 = 0.d0; VF7 = 0.d0
      DF8 = 0.d0; VF8 = 0.d0; DF9 = 0.d0; VF9 = 0.d0
      PE = 0.d0; PET = 0.d0; PETOT = 0.d0; QE = 0.d0; QI = 0.d0
      PI = 0.d0; PIT = 0.d0; PITOT = 0.d0; PEI = 0.d0
      YPELSRS = 0.d0; PEFUS = 0.d0; PIFUS = 0.d0; PEAUX = 0.d0; PIAUX = 0.d0
      PEN = 0.d0; PIN = 0.d0; PJOUL = 0.d0; PRAD = 0.d0
      Sn14 = 0.d0; Sn245 = 0.d0; SCUBM = 0.d0
      PEBM = 0.d0; PIBM = 0.d0; PEICR = 0.d0; PIICR = 0.d0

      allocate(PEECR(n_xcp), CUECR(n_xcp))
      PEECR = 0.d0; CUECR = 0.d0

      allocate(DN(n_xcp), CN(n_xcp), HE(n_xcp), XI(n_xcp), CC(n_xcp), &
               DSI(n_xcp), DSE(n_xcp), DSN(n_xcp))
      DN = 0.d0; CN = 0.d0; HE = 0.d0; XI = 0.d0; CC = 0.d0
      DSI = 0.d0; DSE = 0.d0; DSN = 0.d0

      ! Populate initial profiles
      ne(1:NA1) = cp_in%profiles_1d(i_time)%electrons%density(NABEG:n_xcp) / denA2D
      Te(1:NA1) = cp_in%profiles_1d(i_time)%electrons%temperature(NABEG:n_xcp) / temA2D
      Ti(1:NA1) = cp_in%profiles_1d(i_time)%ion(1)%temperature(NABEG:n_xcp) / temA2D
      MU(1:NA1) = 1. / cp_in%profiles_1d(i_time)%q(NABEG:n_xcp)

      if (cp_in%profiles_1d(i_time)%grid%psi(1) < cp_in%profiles_1d(i_time)%grid%psi(n_xcp)) then
         sign_psi = 1.d0
      else
         sign_psi = -1.d0
      end if

      FP(1:NA1)   = sign_psi * cp_in%profiles_1d(i_time)%grid%psi(NABEG:n_xcp)
      CU(1:NA1)   = abs(cp_in%profiles_1d(i_time)%j_total(NABEG:n_xcp) / cuA2D)
      CUbs(1:NA1) = abs(cp_in%profiles_1d(i_time)%j_bootstrap(NABEG:n_xcp) / cuA2D)
      CD(1:NA1)   = abs(cp_in%profiles_1d(i_time)%j_non_inductive(NABEG:n_xcp) / cuA2D)
      CC(1:NA1)   = abs(cp_in%profiles_1d(i_time)%conductivity_parallel(NABEG:n_xcp) / cuA2D)

      F0 = 0.d0; F1 = 0.d0; F2 = 0.d0; F3 = 0.d0
      ispec(1:n_ion + 3) = 0; Khydr = 0; istep = 0; iH = 0; iD = 0; iT = 0

      do i = 1, n_ion
         ai = cp_in%profiles_1d(i_time)%ion(i)%element(1)%a
         zi = cp_in%profiles_1d(i_time)%ion(i)%element(1)%z_n
         if (zi == 1.d0) then
            if (ai < 1.5d0) then
               ispec(1) = i; Khydr = Khydr + 1
            end if
            if (ai > 1.d0 .and. ai < 3.d0) then
               ispec(2) = i; Khydr = Khydr + 1
            end if
            if (ai > 2.5d0) then
               ispec(3) = i; Khydr = Khydr + 1
            end if
         else
            istep = istep + 1
            ispec(3 + istep) = i
         end if
      end do

      if (Khydr < 1) then
         write(*,*) 'No hydrogen isotopes in the ion list'
         write(*,*) 'Hydrogen species: iH,iD,iT Khydr', iH, iD, iT, Khydr
         stop
      end if
      if (smart_in%sw_stdout /= 0) write(*,200) 'ispec          = ', ispec

      if (ispec(1) /= 0) then
         iH = ispec(1)
         F1(1:NA1) = cp_in%profiles_1d(i_time)%ion(iH)%density(NABEG:n_xcp) / denA2D
         F0(1:NA1) = F0(1:n_xcp) + cp_in%profiles_1d(i_time)%neutral(iH)%density(NABEG:n_xcp) / denA2D
      end if
      if (ispec(2) /= 0) then
         iD = ispec(2)
         F2(1:NA1) = cp_in%profiles_1d(i_time)%ion(iD)%density(NABEG:n_xcp) / denA2D
         F0(1:NA1) = F0(1:n_xcp) + cp_in%profiles_1d(i_time)%neutral(iD)%density(1:n_xcp) / denA2D
      end if
      if (ispec(3) /= 0) then
         iT = ispec(3)
         F3(1:NA1) = cp_in%profiles_1d(i_time)%ion(iT)%density(NABEG:n_xcp) / denA2D
         F0(1:NA1) = F0(NABEG:n_xcp) + cp_in%profiles_1d(i_time)%neutral(iT)%density(NABEG:n_xcp) / denA2D
      end if

      AMAIN(1:NA1) = (F1(1:NA1) + 2.*F2(1:NA1) + 3.*F3(1:NA1)) / &
                     (F1(1:NA1) + F2(1:NA1) + F3(1:NA1))

      do j = 1, NA1
         ni(j)    = 0.d0
         Z2NdA(j) = 0.d0
         ZEF(j)   = 0.d0
         ZMAIN(j) = 0.d0
         do i = 1, n_ion
            ai = cp_in%profiles_1d(i_time)%ion(i)%element(1)%a
            zi = cp_in%profiles_1d(i_time)%ion(i)%z_ion_1D(j)
            ni(j) = ni(j) + cp_in%profiles_1d(i_time)%ion(i)%density(j) / denA2D
            if (ai > 0.) then
               Z2NdA(j) = Z2NdA(j) + cp_in%profiles_1d(i_time)%ion(i)%density(j) / denA2D * zi**2 / ai
            else
               write(*,*) 'wrong ion mass, i=', i
            end if
            ZEF(j) = ZEF(j) + cp_in%profiles_1d(i_time)%ion(i)%density(j) / denA2D * &
                     cp_in%profiles_1d(i_time)%ion(i)%z_ion_1D(j)**2 / ne(j)
            ZMAIN(j) = ZMAIN(j) + &
                       cp_in%profiles_1d(i_time)%ion(i)%density(j) / denA2D * &
                       cp_in%profiles_1d(i_time)%ion(i)%z_ion_1D(j) / ne(j)
         end do
         if (ne(j) >= 0. .and. Te(j) > 0.) then
            PEI(j) = 0.00246 * (15.9 - .5 * dlog(ne(j)) + dlog(Te(j))) * &
                     ne(j) * Z2NdA(j) * (Te(j) - Ti(j)) / Te(j) / dsqrt(Te(j))
         end if
      end do

      RTOR = eq_in%vacuum_toroidal_field%r0
      BTOR = dabs(eq_in%vacuum_toroidal_field%b0(j_time))
      XEQ(1:NE1)  = eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(NEBEG:n_xeq)
      VOLe(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%VOLUME(NEBEG:n_xeq)
      XCP(1:NA1)  = cp_in%profiles_1d(i_time)%grid%rho_tor_norm(NABEG:n_xcp)
      ROC = eq_in%time_slice(j_time)%profiles_1d%phi(n_xeq)
      ROC = dsqrt(dabs(ROC / BTOR / M_PI))
      RHO(1:NA1) = ROC * XCP(1:NA1)

      if (eq_in%time_slice(j_time)%profiles_1d%psi(n_xeq) > eq_in%time_slice(j_time)%profiles_1d%psi(NEBEG)) then
         sign_psi = 1.d0
      else
         sign_psi = -1.d0
      end if
      FPe(1:NE1) = sign_psi * eq_in%time_slice(j_time)%profiles_1d%psi(NEBEG:n_xeq)

      IPOLe(1:NE1) = dabs(eq_in%time_slice(j_time)%profiles_1d%f(NEBEG:n_xeq) / RTOR / BTOR)
      SLATe(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%surface(NEBEG:n_xeq)
      G11e(1:NE1)  = eq_in%time_slice(j_time)%profiles_1d%gm3(NEBEG:n_xeq)
      G33e(1:NE1)  = eq_in%time_slice(j_time)%profiles_1d%gm1(NEBEG:n_xeq) * RTOR**2
      BDB02e(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%gm5(NEBEG:n_xeq) / BTOR**2
      B0DB2e(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%gm4(NEBEG:n_xeq) * BTOR**2
      G22e(1:NE1)   = eq_in%time_slice(j_time)%profiles_1d%gm2(NEBEG:n_xeq)

      if (associated(eq_in%time_slice(j_time)%profiles_1d%b_field_min)) then
         BMINTe(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%b_field_min(NEBEG:n_xeq)
      else
         BMINTe(1:NE1) = BTOR * RTOR / (RTOR + SHIFe(1:NE1) + AMETRe(1:NE1))
      end if
      if (associated(eq_in%time_slice(j_time)%profiles_1d%b_field_max)) then
         BMAXTe(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%b_field_max(NEBEG:n_xeq)
      else
         BMAXTe(1:NE1) = BTOR * RTOR / (RTOR + SHIFe(1:NE1) - AMETRe(1:NE1))
      end if
      if (associated(eq_in%time_slice(j_time)%profiles_1d%b_field_max)) then
         BDB0e(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%b_field_max(NEBEG:n_xeq)
      else
         BDB0e(1:NE1) = sqrt(BDB02e(1:NE1))
      end if

      IPL  = dabs(eq_in%time_slice(j_time)%global_quantities%Ip / cuA2D)
      ABC  = ametre(NE1)
      SHIFT = shife(NE1)
      HRO  = (XCP(3) - XCP(2)) * ROC
      ALFA = 0.001d0

      if (n_xeq /= n_xcp) then
         call SMTH(ALFA, NE1, VOLe, XEQ, NA1, VOL, XCP, nrd)
         call SMTH(ALFA, NE1, ametre, XEQ, NA1, ametr, XCP, nrd)
         call SMTH(ALFA, NE1, shife, XEQ, NA1, shif, XCP, nrd)
         call SMTH(ALFA, NE1, IPOLe, XEQ, NA1, IPOL, XCP, nrd)
         call SMTH(ALFA, NE1, G11e, XEQ, NA1, G11, XCP, nrd)
         call SMTH(ALFA, NE1, G33e, XEQ, NA1, G33, XCP, nrd)
         call SMTH(ALFA, NE1, SLATe, XEQ, NA1, SLAT, XCP, nrd)
         call SMTH(ALFA, n_xeq, BMINTe, XEQ, n_xcp, BMINT, XCP, nrd)
         call SMTH(ALFA, n_xeq, BMAXTe, XEQ, n_xcp, BMAXT, XCP, nrd)
         call SMTH(ALFA, n_xeq, BDB0e, XEQ, n_xcp, BDB0, XCP, nrd)
         call SMTH(ALFA, NE1, BDB02e, XEQ, NA1, BDB02, XCP, nrd)
         call SMTH(ALFA, NE1, B0DB2e, XEQ, NA1, B0DB2, XCP, nrd)
         call SMTH(ALFA, NE1, G22e, XEQ, NA1, G22, XCP, nrd)
      else
         IPOL = IPOLe; VOL = VOLe; shif = shife; ametr = ametre; G11 = G11e; G33 = G33e
         G22 = G22e; SLAT = SLATe; BMINT = BMINTe; BMAXT = BMAXTe; BDB0 = BDB0e; BDB02 = BDB02e
         B0DB2 = B0DB2e; FP = FPe
      end if

      do j = 2, NA1
         VR(j)  = (VOL(j) - VOL(j - 1)) / HRO
         G11(j) = VR(j) * G11(j)
      end do
      VR(1)  = VR(2) / 3.
      G11(1) = VR(1) * G11(1)

      G22(1:NA1) = VR(1:NA1) * G22(1:NA1) / IPOL(1:NA1) * RTOR / 4. / M_PI**2
      SQEPS(1:NA1) = sqrt(AMETR(1:NA1) / (RTOR + SHIF(1:NA1)))
      if (SQEPS(1) <= 1.d-2) SQEPS(1) = SQEPS(2) / 1.4

      do j = 1, NA1
         y = min(.99d0, BTOR * RTOR / BMAXT(j) / (RTOR + SHIF(j)))
         FOFB(j) = B0DB2(j) * (1. - sqrt(1.d0 - y) * (1. + 0.5 * y))
      end do

      TAU  = smart_in%TAU
      dtau = smart_in%dtau

      write(*,*) 'NA1', NA1, NE1

      ! (Diagnostics prints retained as in original code) ...

      ! Conversion to IMAS units (conditional on NABEG)
      cp_out%time(i_time) = TIME
      if (NABEG == 1) then
         cp_out%profiles_1d(i_time)%electrons%density(1:n_xcp)     = ne(1:n_xcp) * denA2D
         cp_out%profiles_1d(i_time)%electrons%temperature(1:n_xcp) = Te(1:n_xcp) * temA2D
         cp_out%profiles_1d(i_time)%q(1:n_xcp)                     = 1. / MU(1:n_xcp)
         cp_out%profiles_1d(i_time)%grid%psi(1:n_xcp)              = sign_psi * FP(1:n_xcp)
         cp_out%profiles_1d(i_time)%j_total(1:n_xcp)               = sign_psi * CU(1:n_xcp) * cuA2D
         cp_out%profiles_1d(i_time)%j_bootstrap(1:n_xcp)           = sign_psi * CUBS(1:n_xcp) * cuA2D
         cp_out%profiles_1d(i_time)%j_non_inductive(NABEG:n_xcp)   = sign_psi * CD(1:n_xcp) * cuA2D
         if (ispec(1) /= 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(1))%density(1:n_xcp)    = F1(1:n_xcp) * denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(1))%density(1:n_xcp) = F01B / F0B * F0(1:n_xcp) * denA2D
         end if
         if (ispec(2) /= 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(2))%density(1:n_xcp)    = F2(1:n_xcp) * denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(2))%density(1:n_xcp) = F02B / F0B * F0(1:n_xcp) * denA2D
         end if
         if (ispec(3) /= 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(3))%density(1:n_xcp)    = F3(1:n_xcp) * denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(3))%density(1:n_xcp) = F03B / F0B * F0(1:n_xcp) * denA2D
         end if
         do i = 1, n_ion
            cp_out%profiles_1d(i_time)%ion(i)%temperature(1:n_xcp) = Ti(1:n_xcp) * temA2D
         end do
      else
         cp_out%profiles_1d(i_time)%electrons%density(NABEG:n_xcp)     = ne(1:NA1) * denA2D
         cp_out%profiles_1d(i_time)%electrons%temperature(NABEG:n_xcp) = Te(1:NA1) * temA2D
         cp_out%profiles_1d(i_time)%q(NABEG:n_xcp)                    = 1. / MU(1:NA1)
         cp_out%profiles_1d(i_time)%grid%psi(NABEG:n_xcp)             = sign_psi * FP(1:NA1)
         cp_out%profiles_1d(i_time)%j_total(NABEG:n_xcp)              = sign_psi * CU(1:NA1) * cuA2D
         cp_out%profiles_1d(i_time)%j_bootstrap(NABEG:n_xcp)          = sign_psi * CUBS(1:NA1) * cuA2D
         cp_out%profiles_1d(i_time)%j_non_inductive(NABEG:n_xcp)      = sign_psi * CD(1:NA1) * cuA2D
         if (ispec(1) /= 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(1))%density(NABEG:n_xcp)    = F1(1:NA1) * denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(1))%density(NABEG:n_xcp) = F01B / F0B * F0(1:NA1) * denA2D
         end if
         if (ispec(2) /= 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(2))%density(NABEG:n_xcp)    = F2(1:NA1) * denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(2))%density(NABEG:n_xcp) = F02B / F0B * F0(1:NA1) * denA2D
         end if
         if (ispec(3) /= 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(3))%density(NEBEG:n_xcp)    = F3(1:NA1) * denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(3))%density(NABEG:n_xcp) = F03B / F0B * F0(1:NA1) * denA2D
         end if
         do i = 1, n_ion
            cp_out%profiles_1d(i_time)%ion(i)%temperature(NABEG:n_xcp) = Ti(1:NA1) * temA2D
         end do
         do j = 1, NABEG
            cp_out%profiles_1d(i_time)%electrons%density(j)     = ne(1) * denA2D
            cp_out%profiles_1d(i_time)%electrons%temperature(j) = Te(1) * temA2D
            cp_out%profiles_1d(i_time)%q(j)                     = 1. / MU(1)
            cp_out%profiles_1d(i_time)%grid%psi(j)              = sign_psi * FP(1)
            cp_out%profiles_1d(i_time)%j_total(j)               = sign_psi * CU(1) * cuA2D
            cp_out%profiles_1d(i_time)%j_bootstrap(j)           = sign_psi * CUBS(1) * cuA2D
            cp_out%profiles_1d(i_time)%j_non_inductive(j)       = sign_psi * CD(1) * cuA2D
            if (ispec(1) /= 0) then
               cp_out%profiles_1d(i_time)%ion(ispec(1))%density(j)    = F1(1) * denA2D
               cp_out%profiles_1d(i_time)%neutral(ispec(1))%density(j) = F01B / F0B * F0(1) * denA2D
            end if
            if (ispec(2) /= 0) then
               cp_out%profiles_1d(i_time)%ion(ispec(2))%density(j)    = F2(1) * denA2D
               cp_out%profiles_1d(i_time)%neutral(ispec(2))%density(j) = F02B / F0B * F0(1) * denA2D
            end if
            if (ispec(3) /= 0) then
               cp_out%profiles_1d(i_time)%ion(ispec(3))%density(j)    = F3(1) * denA2D
               cp_out%profiles_1d(i_time)%neutral(ispec(3))%density(j) = F03B / F0B * F0(1) * denA2D
            end if
            do i = 1, n_ion
               cp_out%profiles_1d(i_time)%ion(i)%temperature(j) = Ti(1) * temA2D
            end do
         end do
      end if

      if (smart_in%sw_stdout /= 0) then
         if (jprint == 1) then
            write(*,100) 'YDABL, YDDEP   = ', YDABL, YDDEP
            write(*,100) 'RHOEC, RHODR   = ', RHOEC, RHODR
            write(*,100) 'ROC, QECR      = ', ROC, QECR
            write(*,100) 'YEFFec         = ', YEFFec
            write(*,100) 'Pecr, Pe, Pi   = ', VINTa(PEECR, ROC, RHO, VR, NA1), &
                                            VINTa(PEECR, ROC, RHO, VR, NA1), &
                                            VINTa(PI, ROC, RHO, VR, NA1)
         end if
      end if

100   format(A17,5(1PE15.6))
200   format(A17,10I5)

      ! Deallocate
      deallocate(key4control)
      deallocate(ispec)
      deallocate(ne, ni, Te, Ti, neX, niX, TeX, TiX, TN, NN, &
                 F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, &
                 F0x, F1x, F2x, F3x, F4x, F5x, F6x, F7x, F8x, F9x, &
                 VTOR, FP, neo, nio, Teo, Tio, F0o, F1o, F2o, F3o, &
                 F4o, F5o, F6o, F7o, F8o, F9o, FPo, cu, cutor, cd, &
                 cubs, UPL, ULON, EZ, ZEF, AMAIN, Z2NdA, ZMAIN, &
                 NHYDR, NDEUT, NTRIT, NALF, NHE3 )

      deallocate(SN, SNN, SNTOT, QN, GN, GNX, PE, PET, PETOT, &
                 SF0, SFF0, SF0TOT, QF0, GF0, GF0X, &
                 SF1, SFF1, SF1TOT, QF1, GF1, GF1X, &
                 SF2, SFF2, SF2TOT, QF2, GF2, GF2X, &
                 SF3, SFF3, SF3TOT, QF3, GF3, GF3X, &
                 SF4, SFF4, SF4TOT, QF4, GF4, GF4X, &
                 SF5, SFF5, SF5TOT, QF5, GF5, GF5X, &
                 SF6, SFF6, SF6TOT, QF6, GF6, GF6X, &
                 SF7, SFF7, SF7TOT, QF7, GF7, GF7X, &
                 SF8, SFF8, SF8TOT, QF8, GF8, GF8X, &
                 SF9, SFF9, SF9TOT, QF9, GF9, GF9X, &
                 PI, PIT, PITOT, PEI, QE, QI, PEECR, CUECR, YPELSRS, &
                 PEFUS, PIFUS, PEAUX, PIAUX, PEN, PIN, PJOUL, PRAD, &
                 PEBM, PIBM, PEICR, PIICR, Sn14, Sn245, SCUBM )

      deallocate(DF0, VF0, DF1, VF1, DF2, VF2, DF3, VF3, DF4, VF4, &
                 DF5, VF5, DF6, VF6, DF7, VF7, DF8, VF8, DF9, VF9, &
                 DN, CN, HE, XI, CC, DSE, DSI, DSN)

      deallocate(IPOL, G11, G33, SLAT, ametr, shif, vr, vro, SQEPS, &
                 mu, VOL, XCP, ametre, shife, VOLe, FPe, XEQ, RHO, &
                 IPOLe, G11e, G33e, SLATe, BMINT, BMAXT, BDB0, &
                 BDB02, B0DB2, FOFB, G22, EQFF, EQPF, FV, &
                 BMINTe, BMAXTe, BDB0e, BDB02e, B0DB2e, FOFBe, G22e)

      call deallocate_smart_data(smart_in)
   end subroutine smart

end module mod_smart

!======================================================================|
