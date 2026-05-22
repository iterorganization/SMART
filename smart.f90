include 'mathematical_constants.f90'

module mod_smart

contains

   subroutine smart(eq_in, cp_in, pellets_in, cp_out, codeparam, &
      error_flag, error_message)

      ! --------------------------------------- VERSION 08-OCT-2025
      ! PELLET ABLATION MODEL FROM ASTRA: SMART
      ! IDS INPUT: EQUILIBRIUM, CORE_PROFILES
      ! IDS OUTPUT: CORE_PROFILES
      ! ---------------------------------------

      use ids_schemas, only: ids_equilibrium, ids_core_profiles, ids_pellets
      use ids_schemas, only: ids_parameters_input, ids_is_valid
      use ids_routines, only: ids_copy
      use mod_codeparam_smart
      use mathematical_constants, only: M_PI

      implicit none

      type(ids_equilibrium) :: eq_in
      type(ids_core_profiles) :: cp_in, cp_out
      type(ids_pellets) :: pellets_in
      type(ids_parameters_input) :: codeparam
      type(type_smart_data) :: smart_in
      integer, intent(out) :: error_flag
      logical:: from_pellets_ids = .False.
      character(len=:), pointer, intent(out) :: error_message

      integer :: i, j, i_time, j_time, n_xcp, n_xeq, n_ion, nrd, NA1, jtime, jend
      integer:: iH, iD, iT, na, nb1, Khydr, JABS, istep, NABEG, NEBEG, NE1
      integer:: jheat, jdens, jpsi, jpel, jprint, nbnd1, nbnd2, IMIX, jmix
      double precision &
         denA2D, temA2D, presA2D, cuA2D, &
         HRO, ROC, BTOR, SHIFT, ABC, RTOR, ALFA, &
         YDABL, YDDEP, ai, zi, Dtpel, &
         NEB, TEB, TIB, F0B, F1B, F2B, F3B, F4B, F5B, F6B, F7B, F8B, F9B, F01B, F02B, F03B, &
         GN2E, GN2I, SVRCy, SVCXy, SVIEy, RHOEC, RHODR, QECR, YEFFec, &
         y, y1, y2, y3, y4, pdt, PAIONy, VINTa, QNB, IINTa, IPL, &
         SVD1y, SVD2y, SVDBHy, stbrny, PENLIy, PDD1, PDD2, PDD3, TAUSp, FTLLMRy
      double precision &
         elon, trian, sign_psi
      !
      double precision &
         TAU, dtau, TIME, TIMBEG
      double precision, save :: TIMPEL = 0.d0

      integer, allocatable :: ispec(:)
      !
      double precision, allocatable :: &
         ne(:), ni(:), Te(:), Ti(:), TN(:), &
         F0(:), F1(:), F2(:), F3(:), F4(:), F5(:), F6(:), F7(:), F8(:), F9(:), FP(:), &
         neo(:), nio(:), Teo(:), Tio(:), &
         F0o(:), F1o(:), F2o(:), F3o(:), F4o(:), F5o(:), F6o(:), F7o(:), F8o(:), F9o(:), FPo(:), &
         neX(:), niX(:), TeX(:), TiX(:), &
         F0x(:), F1x(:), F2x(:), F3x(:), F4x(:), F5x(:), F6x(:), F7x(:), F8x(:), F9x(:), NN(:), &
         cu(:), cutor(:), cd(:), cubs(:), UPL(:), ULON(:), EZ(:), FV(:), &
         ZEF(:), AMAIN(:), Z2NdA(:), ZMAIN(:), SQEPS(:), PBLON(:), PBPER(:), PFAST(:), &
         F1fast(:),  F2fast(:),  F3fast(:),  F4fast(:), F5fast(:), VTOR(:), &
         NHYDR(:),NDEUT(:),NTRIT(:),NALF(:),NHE3(:)

      double precision, allocatable :: &
         SN(:), SNN(:), SNTOT(:), QN(:), GN(:), GNX(:), QE(:), QI(:), &
         SF0(:), SFF0(:), SF0TOT(:), QF0(:), GF0(:), GF0X(:), &
         SF1(:), SFF1(:), SF1TOT(:), QF1(:), GF1(:), GF1X(:), &
         SF2(:), SFF2(:), SF2TOT(:), QF2(:), GF2(:), GF2X(:), &
         SF3(:), SFF3(:), SF3TOT(:), QF3(:), GF3(:), GF3X(:), &
         SF4(:), SFF4(:), SF4TOT(:), QF4(:), GF4(:), GF4X(:), &
         SF5(:), SFF5(:), SF5TOT(:), QF5(:), GF5(:), GF5X(:), &
         SF6(:), SFF6(:), SF6TOT(:), QF6(:), GF6(:), GF6X(:), &
         SF7(:), SFF7(:), SF7TOT(:), QF7(:), GF7(:), GF7X(:), &
         SF8(:), SFF8(:), SF8TOT(:), QF8(:), GF8(:), GF8X(:), &
         SF9(:), SFF9(:), SF9TOT(:), QF9(:), GF9(:), GF9X(:), &
         PE(:), PET(:), PETOT(:), PI(:), PIT(:), PITOT(:), PEI(:), &
         PEECR(:), CUECR(:), YPELSRS(:), &
         PEFUS(:), PIFUS(:), PEAUX(:), PIAUX(:), PEN(:), PIN(:), PJOUL(:), PRAD(:), &
         Sn14(:), Sn245(:), SCUBM(:), PEICR(:), PIICR(:), PEBM(:), PIBM(:)

      double precision, allocatable :: &
         DF0(:), VF0(:), DF1(:), VF1(:), DF2(:), VF2(:), DF3(:), VF3(:), &
         DF4(:), VF4(:), DF5(:), VF5(:), DF6(:), VF6(:), DF7(:), VF7(:), &
         DF8(:), VF8(:), DF9(:), VF9(:), DN(:), CN(:), HE(:), XI(:), &
         cc(:), DSI(:), DSE(:), DSN(:)

      double precision, allocatable :: &
         ametr(:), shif(:), mu(:), VOL(:), vr(:), VOLo(:), vro(:), &
         IPOL(:), G11(:), G33(:), SLAT(:), &
         BMINT(:), BMAXT(:), BDB0(:), BDB02(:), B0DB2(:), &
         FOFB(:), G22(:), EQFF(:), EQPF(:)
      double precision, allocatable :: &
         ametre(:), shife(:), VOLe(:), XCP(:), XEQ(:), RHO(:), FPe(:), &
         IPOLe(:), G11e(:), G33e(:), SLATe(:), &
         BMINTe(:), BMAXTe(:), BDB0e(:), BDB02e(:), B0DB2e(:), FOFBe(:), G22e(:)
      character*14 ARRNAME(30), VARNAME(30)

      external SVRCy, SVCXy, SVIEy, PAIONy, VINTa, IINTa, SVD1y, SVD2y, SVDBHy, stbrny, PENLIy
      external SMTH, ECH2a, STEPUPN, STEPUPN0, STEPUPT, STEPUPF, CUBSy, RHSEQy, ARR, PelIMAS1, FTLLMRy, MIXALL

      data ARRNAME/' x ',' Te ',' Ti ',' ne ',' ni ', &
         ' nHth ',' nHf ',' nDth ',' nTth ',' nTf ',' n4Heth ',' n4Hef ', ' n3Heth ',' n3Hef ', &
         ' U ',' q ',' Fp ',' Jtot',' Jbs ',' Jec ',' Zeff ',' Zmain ',' Amain ', &
         ' Pe ',' Pi ',' Pefus ',' Pifus ',' Poh ',' Pec ',' Pei ' &
         /
      data VARNAME/' time ',' Te0 ',' Ti0 ',' ne0 ',' ni0 ',' Tea ',' Tia ',' nea ',' nia ', &
         ' n00 ',' n0a ',' <ne> ', ' Ibs ',' Itot ',' Icd ',' U0 ',' Ua ',' Fp0 ',' Fp ', &
         ' Pe ',' Pi ',' Pefus ',' Pifus ',' Poh ',' Pec ',' Pei ', &
         ' q(0) ',' q(a) ',' <Shdt> ',' <Sn0> '  &
         /

      ! INITIALISATION OF ERROR FLAG
      error_flag = 0
      nullify (error_message) ! do this otherwise gfortran behaviour is undefined

      ! ASTRA- IMAS units transfer
      !       include 'declar.units'
      denA2D = 1.d19
      temA2D = 1.d3
      presA2D = 1.6d3 ! in ASTRA pressure is keV 10^19m-3, have to change Pa
      cuA2D = 1.d6

      ! Code specific parameters
      !       include 'init.pel'
      ! CHECK IF INPUT PELLETS IDS IS FILLED: IF SO, USE IT INSTEAD OF CODEPARAM
      ! FOR PELLET DESCRIPTION
      if (pellets_in%ids_properties%homogeneous_time .ge. 0) from_pellets_ids = .True.

      ! CHECK IF INPUT IDS IS VALID
      if (ids_is_valid(eq_in%ids_properties%homogeneous_time) .and. &
         size(eq_in%time) > 0 .and. &
         ids_is_valid(cp_in%ids_properties%homogeneous_time) .and. &
         size(cp_in%time) > 0) then

         call assign_codeparam(codeparam%parameters_value, smart_in)

         if (from_pellets_ids .eqv. .true.) then ! REPLACE PELLET INFORMATION
            write (*, *) 'Input pellets IDS detected'
            smart_in%YAM = pellets_in%time_slice(1)%pellet(1)%species(1)%a
            smart_in%YVP = pellets_in%time_slice(1)%pellet(1)%velocity_initial*1.e-3
            smart_in%YVOL = pellets_in%time_slice(1)%pellet(1)%shape%size(1)**2* &
               pellets_in%time_slice(1)%pellet(1)%shape%size(2)*M_PI*1.e+9
         else
            if (smart_in%sw_stdout .ne. 0) then
               write (*, *) 'Input pellets IDS NOT detected'
            end if
         end if

         ! COPY THE INPUT IDS IN THE OUTPUT IDS
         call ids_copy(cp_in, cp_out)
      else
         error_flag = -1
         allocate (character(50):: error_message)
         error_message = 'Error in SMART: input IDS not valid'
         return
      end if

      !======================
      ! initialization
      i_time = size(cp_in%time)
      TIMBEG = cp_in%time(i_time)
      if (smart_in%sw_stdout .ne. 0) then
         write(*, '(/,A17, I5, F10.4)')'i_time, TIMBEG = ', i_time, TIMBEG
      endif
      !write (*, *) 'i_time,time=', i_time, TIMBEG
      !              i_time=size(cp_in%profiles_1d(:)%time)
      j_time = size(eq_in%time_slice)
      n_xcp = size(cp_in%profiles_1d(i_time)%grid%rho_tor_norm)
      n_ion = size(cp_in%profiles_1d(i_time)%ion)
      n_xeq = size(eq_in%time_slice(j_time)%profiles_1d%psi)
      nrd = 2*max(n_xcp,n_xeq)
! temporary vvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      nbnd1 = size(eq_in%time_slice(j_time)%boundary%outline%r)
      jprint = smart_in%key4control(16)
   if (jprint.gt.0) then
   write(*,*), 'nbnd1=',nbnd1
!
      open(20,file='IMAS.dat')
      write(20,*) i_time, nbnd1
      write(20,*) j_time, timbeg
      do j=1,nbnd1
         write(20,*) eq_in%time_slice(j_time)%boundary%outline%r(j), &
            eq_in%time_slice(j_time)%boundary%outline%z(j)
      enddo
      close(20)
      endif
!
      elon  =eq_in%time_slice(j_time)%boundary%elongation
      trian =eq_in%time_slice(j_time)%boundary%triangularity
! temporary  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! elimination of the fake -1/2 point for JINTRAC output
      if(cp_in%profiles_1d(i_time)%grid%rho_tor_norm(1).le.0.d0) then
         NABEG=2
         NA1=n_xcp-1
      else
         NABEG=1
         NA1=n_xcp
      endif
      if(eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(1).le.0.d0) then
         NEBEG=2
         NE1=n_xeq-1
      else
         NEBEG=1
         NE1=n_xeq
      endif

      NB1 = NA1
      NA = NA1 - 1
      !write (*, *) 'NA1,nrd,n_xeq', NA1, nrd, n_xeq
      !==========================
      allocate (ispec(n_ion + 3))
      !==========================
      !        include 'alloc.corprf'
      allocate ( &
         ne(n_xcp), ni(n_xcp), Te(n_xcp), Ti(n_xcp), &
         nex(n_xcp), nix(n_xcp), TEX(n_xcp), TIX(n_xcp), &
         F1x(n_xcp), F2x(n_xcp), F3x(n_xcp), F4x(n_xcp), F5x(n_xcp),&
         F6x(n_xcp), F7x(n_xcp), F8x(n_xcp), F9x(n_xcp), &
         F1(n_xcp), F2(n_xcp), F3(n_xcp), F4(n_xcp), F5(n_xcp), &
         F6(n_xcp), F7(n_xcp), F8(n_xcp), F9(n_xcp), FP(n_xcp), FV(n_xcp),&
         neo(n_xcp), nio(n_xcp), Teo(n_xcp), Tio(n_xcp), &
         F1o(n_xcp), F2o(n_xcp), F3o(n_xcp), F4o(n_xcp), F5o(n_xcp), &
         F6o(n_xcp), F7o(n_xcp), F8o(n_xcp), F9o(n_xcp),  FPo(n_xcp), VTOR(n_xcp),&
         F0(n_xcp), F0o(n_xcp), F0x(n_xcp), NN(n_xcp), TN(n_xcp), &
         cu(n_xcp), cutor(n_xcp), cd(n_xcp), cubs(n_xcp), &
         UPL(n_xcp), ULON(n_xcp), EZ(n_xcp), ZEF(n_xcp), AMAIN(n_xcp), Z2NdA(n_xcp), ZMAIN(n_xcp), &
         PBLON(n_xcp), PBPER(n_xcp), PFAST(n_xcp), &
         F1fast(n_xcp),  F2fast(n_xcp),  F3fast(n_xcp),  F4fast(n_xcp), F5fast(n_xcp), &
         NHYDR(n_xcp), NDEUT(n_xcp),NTRIT(n_xcp), NALF(n_xcp), NHE3(n_xcp) &
         )
      !        include 'alloc.eq'
      ne = 0.d0
      ni = 0.d0
      Te = 0.d0
      Ti = 0.d0
      nex = 0.d0
      nix = 0.d0
      TEX = 0.d0
      TIX = 0.d0
      F1x = 0.d0
      F2x = 0.d0
      F3x = 0.d0
      F4x = 0.d0
      F5x = 0.d0
      F6x = 0.d0
      F7x = 0.d0
      F8x = 0.d0
      F9x = 0.d0
      F1 = 0.d0
      F2 = 0.d0
      F3 = 0.d0
      F4 = 0.d0
      F5 = 0.d0
      F6 = 0.d0
      F7 = 0.d0
      F8 = 0.d0
      F9 = 0.d0
      FP = 0.d0
      neo = 0.d0
      nio = 0.d0
      Teo = 0.d0
      Tio = 0.d0
      F1o = 0.d0
      F2o = 0.d0
      F3o = 0.d0
      F4o = 0.d0
      F5o = 0.d0
      F6o = 0.d0
      F7o = 0.d0
      F8o = 0.d0
      F9o = 0.d0
      FPo = 0.d0
      F0 = 0.d0
      F0o = 0.d0
      F0x = 0.d0
      NN = 0.d0
      TN = 0.d0
      VTOR = 0.d0
      cu = 0.d0
      cutor = 0.d0
      cd = 0.d0
      cubs = 0.d0
      UPL = 0.d0
      ULON = 0.d0
      EZ = 0.d0
      ZEF = 0.d0
      AMAIN = 1.d0
      Z2NdA = 1.d0
      FV=0.d0
      ZMAIN = 1.d0
      NHYDR =0.d0
      NDEUT =0.d0
      NALF =0.d0
      NHE3 =0.d0
      F1fast=0.d0
      F2fast=0.d0
      F3fast=0.d0
      F4fast=0.d0
      F5fast=0.d0

      allocate ( &
         ametr(n_xcp), shif(n_xcp), vr(n_xcp), vro(n_xcp), &
         mu(n_xcp), VOL(n_xcp), VOLo(n_xcp), XCP(n_xcp), &
         IPOL(n_xcp), G11(n_xcp), G33(n_xcp), SLAT(n_xcp), RHO(n_xcp), SQEPS(n_xcp), &
         BMINT(n_xcp), BMAXT(n_xcp), BDB0(n_xcp), BDB02(n_xcp), B0DB2(n_xcp), &
         FOFB(n_xcp), G22(n_xcp), EQFF(n_xcp), EQPF(n_xcp))
      ametr = 0.d0
      shif = 0.d0
      vr = 0.d0
      vro = 0.d0
      mu = 0.d0
      VOL = 0.d0
      VOLo = 0.d0
      XCP = 0.d0
      IPOL = 0.d0
      G11 = 0.d0
      G33 = 0.d0
      SLAT = 0.d0
      BMINT = 0.d0
      BMAXT = 0.d0
      BDB0 = 0.d0
      BDB02 = 0.d0
      B0DB2 = 0.d0
      FOFB = 0.d0
      G22 = 0.d0
      EQFF = 0.d0
      EQPF = 0.d0

      allocate ( &
         ametre(n_xeq), shife(n_xeq), VOLe(n_xeq), &
         FPe(n_xeq), XEQ(n_xeq),  &
         IPOLe(n_xeq), G11e(n_xeq), G33e(n_xeq), SLATe(n_xeq), &
         BMINTe(n_xeq), BMAXTe(n_xeq), BDB0e(n_xeq), BDB02e(n_xeq), B0DB2e(n_xeq), &
         FOFBe(n_xeq), G22e(n_xeq))
      !        include 'alloc.corsrs'
      ametre = 0.d0
      shife = 0.d0
      VOLe = 0.d0
      FPe = 0.d0
      XEQ = 0.d0
      RHO = 0.d0
      IPOLe = 0.d0
      G11e = 0.d0
      G33e = 0.d0
      SLATe = 0.d0
      BMINTe = 0.d0
      BMAXTe = 0.d0
      BDB0e = 0.d0
      BDB02e = 0.d0
      B0DB2e = 0.d0
      FOFBe = 0.d0

      allocate ( &
         SN(n_xcp), SNN(n_xcp), SNTOT(n_xcp), &
         QN(n_xcp), GN(n_xcp), GNX(n_xcp), &
         SF0(n_xcp), SFF0(n_xcp), SF0TOT(n_xcp), &
         QF0(n_xcp), GF0(n_xcp), GF0X(n_xcp), &
         SF1(n_xcp), SFF1(n_xcp), SF1TOT(n_xcp), &
         QF1(n_xcp), GF1(n_xcp), GF1X(n_xcp), &
         SF2(n_xcp), SFF2(n_xcp), SF2TOT(n_xcp), &
         QF2(n_xcp), GF2(n_xcp), GF2X(n_xcp), &
         SF3(n_xcp), SFF3(n_xcp), SF3TOT(n_xcp), &
         QF3(n_xcp), GF3(n_xcp), GF3X(n_xcp), &
         SF4(n_xcp), SFF4(n_xcp), SF4TOT(n_xcp), &
         QF4(n_xcp), GF4(n_xcp), GF4X(n_xcp), &
         SF5(n_xcp), SFF5(n_xcp), SF5TOT(n_xcp), &
         QF5(n_xcp), GF5(n_xcp), GF5X(n_xcp), &
         SF6(n_xcp), SFF6(n_xcp), SF6TOT(n_xcp), &
         QF6(n_xcp), GF6(n_xcp), GF6X(n_xcp), &
         SF7(n_xcp), SFF7(n_xcp), SF7TOT(n_xcp), &
         QF7(n_xcp), GF7(n_xcp), GF7X(n_xcp), &
         SF8(n_xcp), SFF8(n_xcp), SF8TOT(n_xcp), &
         QF8(n_xcp), GF8(n_xcp), GF8X(n_xcp), &
         SF9(n_xcp), SFF9(n_xcp), SF9TOT(n_xcp), &
         QF9(n_xcp), GF9(n_xcp), GF9X(n_xcp), &
         DF0(n_xcp), VF0(n_xcp), DF1(n_xcp), VF1(n_xcp), &
         DF2(n_xcp), VF2(n_xcp), DF3(n_xcp), VF3(n_xcp), &
         DF4(n_xcp), VF4(n_xcp), DF5(n_xcp), VF5(n_xcp), &
         DF6(n_xcp), VF6(n_xcp), DF7(n_xcp), VF7(n_xcp), &
         DF8(n_xcp), VF8(n_xcp), DF9(n_xcp), VF9(n_xcp), &
         PE(n_xcp), PET(n_xcp), PETOT(n_xcp), QE(n_xcp), QI(n_xcp), &
         PI(n_xcp), PIT(n_xcp), PITOT(n_xcp), PEI(n_xcp), &
         PEFUS(n_xcp), PIFUS(n_xcp), PEAUX(n_xcp), PIAUX(n_xcp), &
         PEN(n_xcp), PIN(n_xcp), PJOUL(n_xcp), PRAD(n_xcp), &
         Sn14(n_xcp), Sn245(n_xcp), SCUBM(n_xcp), &
         PEICR(n_xcp), PIICR(n_xcp), PEBM(n_xcp), PIBM(n_xcp), &
         YPELSRS(n_xcp) )
      IMIX =0  !
      SN = 0.d0
      SNN = 0.d0
      SNTOT = 0.d0
      QN = 0.d0
      GN = 0.d0
      GNX = 0.d0
      SF0 = 0.d0
      SFF0 = 0.d0
      SF0TOT = 0.d0
      QF0 = 0.d0
      GF0 = 0.d0
      GF0X = 0.d0
      SF1 = 0.d0
      SFF1 = 0.d0
      SF1TOT = 0.d0
      QF1 = 0.d0
      GF1 = 0.d0
      GF1X = 0.d0
      SF2 = 0.d0
      SFF2 = 0.d0
      SF2TOT = 0.d0
      QF2 = 0.d0
      GF2 = 0.d0
      GF2X = 0.d0
      SF3 = 0.d0
      SFF3 = 0.d0
      SF3TOT = 0.d0
      QF3 = 0.d0
      GF3 = 0.d0
      GF3X = 0.d0
      SF4 = 0.d0
      SFF4 = 0.d0
      SF4TOT = 0.d0
      QF4 = 0.d0
      GF4 = 0.d0
      GF4X = 0.d0
      SF5 = 0.d0
      SFF5 = 0.d0
      SF5TOT = 0.d0
      QF5 = 0.d0
      GF5 = 0.d0
      GF5X = 0.d0
      SF6 = 0.d0
      SFF6 = 0.d0
      SF6TOT = 0.d0
      QF6 = 0.d0
      GF6 = 0.d0
      GF6X = 0.d0
      SF7 = 0.d0
      SFF7 = 0.d0
      SF7TOT = 0.d0
      QF7 = 0.d0
      GF7 = 0.d0
      GF7X = 0.d0
      SF8 = 0.d0
      SFF8 = 0.d0
      SF8TOT = 0.d0
      QF8 = 0.d0
      GF8 = 0.d0
      GF8X = 0.d0
      SF9 = 0.d0
      SFF9 = 0.d0
      SF9TOT = 0.d0
      QF9 = 0.d0
      GF9 = 0.d0
      GF9X = 0.d0
      DF0 = 0.d0
      VF0 = 0.d0
      DF1 = 0.d0
      VF1 = 0.d0
      DF2 = 0.d0
      VF2 = 0.d0
      DF3 = 0.d0
      VF3 = 0.d0
      DF4 = 0.d0
      VF4 = 0.d0
      DF5 = 0.d0
      VF5 = 0.d0
      DF6 = 0.d0
      VF6 = 0.d0
      DF7 = 0.d0
      VF7 = 0.d0
      DF8 = 0.d0
      VF8 = 0.d0
      DF9 = 0.d0
      VF9 = 0.d0
      PE = 0.d0
      PET = 0.d0
      PETOT = 0.d0
      QE = 0.d0
      QI = 0.d0
      PI = 0.d0
      PIT = 0.d0
      PITOT = 0.d0
      PEI = 0.d0
      YPELSRS = 0.d0
      PEFUS = 0.d0
      PIFUS = 0.d0
      PEAUX = 0.d0
      PIAUX = 0.d0
      PEN = 0.d0
      PIN = 0.d0
      PJOUL = 0.d0
      PRAD = 0.d0
      Sn14 = 0.d0
      Sn245 = 0.d0
      SCUBM = 0.d0
      PEBM = 0.d0
      PIBM = 0.d0
      PEICR = 0.d0
      PIICR = 0.d0

      allocate ( PEECR(n_xcp), CUECR(n_xcp) )
      !        include 'alloc.cortran'
      PEECR = 0.d0
      CUECR = 0.d0

      allocate ( DN(n_xcp), CN(n_xcp), HE(n_xcp), XI(n_xcp), CC(n_xcp), &
         DSI(n_xcp), DSE(n_xcp), DSN(n_xcp) )
      DN = 0.d0
      CN = 0.d0
      HE = 0.d0
      XI = 0.d0
      CC = 0.d0
      DSI = 0.d0
      DSE = 0.d0
      DSN = 0.d0


      ne(1:NA1) = cp_in%profiles_1d(i_time)%electrons%density(NABEG:n_xcp)/denA2D
      Te(1:NA1) = cp_in%profiles_1d(i_time)%electrons%temperature(NABEG:n_xcp)/temA2D
      Ti(1:NA1) = cp_in%profiles_1d(i_time)%ion(1)%temperature(NABEG:n_xcp)/temA2D
      MU(1:NA1) = 1./cp_in%profiles_1d(i_time)%q(NABEG:n_xcp)
      if(cp_in%profiles_1d(i_time)%grid%psi(1).lt.cp_in%profiles_1d(i_time)%grid%psi(n_xcp)) then
         sign_psi=1.d0
      else
         sign_psi=-1.d0
      endif
      FP(1:NA1) = sign_psi*cp_in%profiles_1d(i_time)%grid%psi(NABEG:n_xcp)
      CU(1:NA1) = abs(cp_in%profiles_1d(i_time)%j_total(NABEG:n_xcp)/cuA2D)
      CUbs(1:NA1) = abs(cp_in%profiles_1d(i_time)%j_bootstrap(NABEG:n_xcp)/cuA2D)
      CD(1:NA1) = abs(cp_in%profiles_1d(i_time)%j_non_inductive(NABEG:n_xcp)/cuA2D)
      CC(1:NA1) = abs(cp_in%profiles_1d(i_time)%conductivity_parallel(NABEG:n_xcp)/cuA2D)
!      UPL(1:NA1) = cp_in%profiles_1d(i_time)%e_field%toroidal(NABEG:n_xcp)*RTOR*3.1416*2.
! write(*,*) 'J0, CC0 DB', cp_in%profiles_1d(i_time)%j_total(1)/cuA2D, &
!      cp_in%profiles_1d(i_time)%conductivity_parallel(1)/cuA2D
      !============================== detect hydrogen isotopes
      F0 = 0.d0
      F1 = 0.d0
      F2 = 0.d0
      F3 = 0.d0
      ispec(1:n_ion + 3) = 0
      Khydr = 0
      istep = 0
      iH = 0
      iD = 0
      iT = 0
      do i = 1, n_ion
         ai = cp_in%profiles_1d(i_time)%ion(i)%element(1)%a
         zi = cp_in%profiles_1d(i_time)%ion(i)%element(1)%z_n
         if (zi .eq. 1.d0) then
            if (ai .lt. 1.5d0) then
               ispec(1) = i
               Khydr = Khydr + 1
            end if
            if (ai .gt. 1.d0 .and. ai .lt. 3.d0) then
               ispec(2) = i
               Khydr = Khydr + 1
            end if
            if (ai .gt. 2.5d0) then
               ispec(3) = i
               Khydr = Khydr + 1
            end if
         else
            istep = istep + 1
            ispec(3 + istep) = i
         end if
      end do
      if (Khydr .lt. 1) then
         write (*, *) 'No hydrogen isotopes in the ion list'
         write (*, *) 'Hyrogen species: iH,iD,iT Khydr', iH, iD, iT, Khydr
         stop
      end if
      if (smart_in%sw_stdout .ne.0) then
         write (*, 200) 'ispec          = ', ispec
      end if
      if (ispec(1) .ne. 0) then
         iH = ispec(1)
         F1(1:NA1) = cp_in%profiles_1d(i_time)%ion(iH)%density(NABEG:n_xcp)/denA2D
         F0(1:NA1) = F0(1:n_xcp) + cp_in%profiles_1d(i_time)%neutral(iH)%density(NABEG:n_xcp)/denA2D
      end if
      if (ispec(2) .ne. 0) then
         iD = ispec(2)
         F2(1:NA1) = cp_in%profiles_1d(i_time)%ion(iD)%density(NABEG:n_xcp)/denA2D
         F0(1:NA1) = F0(1:n_xcp) + cp_in%profiles_1d(i_time)%neutral(iD)%density(1:n_xcp)/denA2D
      end if
      if (ispec(3) .ne. 0) then
         iT = ispec(3)
         F3(1:NA1) = cp_in%profiles_1d(i_time)%ion(iT)%density(NABEG:n_xcp)/denA2D
         F0(1:NA1) = F0(NABEG:n_xcp) + cp_in%profiles_1d(i_time)%neutral(iT)%density(NABEG:n_xcp)/denA2D
      end if
      AMAIN(1:NA1) = (F1(1:NA1) + 2.*F2(1:NA1) + 3.*F3(1:NA1))/ &
      & (F1(1:NA1) + F2(1:NA1) + F3(1:NA1))
      !!!!!!!!!!!!!!!!!!!!!!!!!! separate ions, av mass, zeff, pei
      do j = 1, NA1
         ni(j) = 0.d0
         Z2NdA(j) = 0.d0
         ZEF(j) = 0.d0
         ZMAIN(j) = 0.d0
         do i = 1, n_ion
            ai = cp_in%profiles_1d(i_time)%ion(i)%element(1)%a
            zi = cp_in%profiles_1d(i_time)%ion(i)%z_ion_1D(j)
            ni(j) = ni(j) + cp_in%profiles_1d(i_time)%ion(i)%density(j)/denA2D
            if (ai .le. 0.) then
               write (*, *) 'wrong ion mass, i=', i
            else
               Z2NdA(j) = Z2NdA(j) + cp_in%profiles_1d(i_time)%ion(i)%density(j)/denA2D*zi**2/ai

            end if
            ZEF(j) = ZEF(j) + cp_in%profiles_1d(i_time)%ion(i)%density(j)/denA2D*&
            & cp_in%profiles_1d(i_time)%ion(i)%z_ion_1D(j)**2/ne(j)
            ZMAIN(j) = Zmain(J)+ &
               cp_in%profiles_1d(i_time)%ion(i)%density(j)/denA2D &
               *cp_in%profiles_1d(i_time)%ion(i)%z_ion_1D(j)/ne(j)
         end do                                        ! j ion
         if (ne(j) .ge. 0. .and. te(j) .gt. 0.) then
            pei(j) = 0.00246*(15.9 - .5*dlog(NE(j)) + dlog(TE(j))) &
               *NE(j)*Z2NdA(j)*(TE(j)-TI(j))/TE(j)/dsqrt(TE(j))
         else
            !                write(*,*) 'ne,te,j',ne(j),te(j),j
         end if
      end do ! na1

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! equilibrium
!      RTOR = (eq_in%time_slice(j_time)%profiles_1d%r_outboard(n_xeq) + &
!           &  eq_in%time_slice(j_time)%profiles_1d%r_inboard(n_xeq))/2.
!      if (RTOR .gt. 0.) then
!         BTOR = dabs(eq_in%vacuum_toroidal_field%b0(i_time))*eq_in%vacuum_toroidal_field%r0/RTOR
!      else
!         write (*, *) 'wrong RTOR =', RTOR
!         stop
!      end if
      RTOR = eq_in%vacuum_toroidal_field%r0
      BTOR = dabs(eq_in%vacuum_toroidal_field%b0(j_time))
      XEQ(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(NEBEG:n_xeq)
      VOLe(1:NE1) =  eq_in%time_slice(j_time)%profiles_1d%VOLUME(NEBEG:n_xeq)
      XCP(1:NA1) = cp_in%profiles_1d(i_time)%grid%rho_tor_norm(NABEG:n_xcp)
      ROC = eq_in%time_slice(j_time)%profiles_1d%phi(n_xeq)
      ROC = dsqrt(dabs(ROC/BTOR/M_PI))
      RHO(1:NA1) = ROC*XCP(1:NA1)
 if(jprint.gt.0) then
         write(*,*) 'Xe 1,2,3,NA,NA1',eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(1), &
            eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(2), &
            eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(3), &
            eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(n_xeq-1), &
            eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(n_xeq)
         write(*,*) 'X 1,2,3,NA,NA1',cp_in%profiles_1d(i_time)%grid%rho_tor_norm(1), &
            cp_in%profiles_1d(i_time)%grid%rho_tor_norm(2), &
            cp_in%profiles_1d(i_time)%grid%rho_tor_norm(3), &
            cp_in%profiles_1d(i_time)%grid%rho_tor_norm(n_xcp-1), &
            cp_in%profiles_1d(i_time)%grid%rho_tor_norm(n_xcp)
 endif
      ametre(1:NE1) =  (eq_in%time_slice(j_time)%profiles_1d%r_outboard(NEBEG:n_xeq) - &
      & eq_in%time_slice(j_time)%profiles_1d%r_inboard(NEBEG:n_xeq))/2.
      shife(1:NE1) = (eq_in%time_slice(j_time)%profiles_1d%r_outboard(NEBEG:n_xeq) + &
      & eq_in%time_slice(j_time)%profiles_1d%r_inboard(NEBEG:n_xeq))/2.-RTOR
!         FPe(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%psi(1:n_xeq)
! temporary
      if(eq_in%time_slice(j_time)%profiles_1d%psi(n_xeq).gt.eq_in%time_slice(j_time)%profiles_1d%psi(NEBEG)) then
         sign_psi=1.d0
      else
         sign_psi=-1.d0
      endif
      FPe(1:NE1) = sign_psi* eq_in%time_slice(j_time)%profiles_1d%psi(NEBEG:n_xeq)
!+ eq_in%time_slice(j_time)%profiles_1d%psi(1)
      IPOLe(1:NE1) = dabs(eq_in%time_slice(j_time)%profiles_1d%f(NEBEG:n_xeq)/RTOR/BTOR)
      SLATe(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%surface(NEBEG:n_xeq)
      G11e(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%gm3(NEBEG:n_xeq)
      G33e(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%gm1(NEBEG:n_xeq)*RTOR**2
      BDB02e(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%gm5(NEBEG:n_xeq)/BTOR**2
      B0DB2e(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%gm4(NEBEG:n_xeq)*BTOR**2
      G22e(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%gm2(NEBEG:n_xeq)

      if(associated(eq_in%time_slice(j_time)%profiles_1d%b_field_min)) then
         BMINTe(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%b_field_min(NEBEG:n_xeq)
      else
         BMINTe(1:NE1) = BTOR*RTOR/(RTOR+SHIFe(1:NE1)+AMETRe(1:NE1))
      endif
      if(associated(eq_in%time_slice(j_time)%profiles_1d%b_field_max)) then
         BMAXTe(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%b_field_max(NEBEG:n_xeq)
      else
         BMAXTe(1:NE1) = BTOR*RTOR/(RTOR+SHIFe(1:NE1)-AMETRe(1:NE1))
      endif
      if(associated(eq_in%time_slice(j_time)%profiles_1d%b_field_max)) then
         BDB0e(1:NE1) = eq_in%time_slice(j_time)%profiles_1d%b_field_max(NEBEG:n_xeq)
      else
         BDB0e(1:NE1) = sqrt(BDB02e(1:NE1))
      endif

      IPL = dabs(eq_in%time_slice(j_time)%global_quantities%Ip/cuA2D)
      ABC = ametre(NE1)
      SHIFT = Shife(NE1) ! - RTOR
      HRO = (XCP(3) - XCP(2))*ROC
      ALFA = 0.001d0
      !SMTH(ALFA,NO,FO,XO,N,FN,XN) mapping from EQ to CP grids
      !        write(*,*) 'trace 184'
      if (n_xeq .ne. n_xcp) then
         call SMTH(ALFA, NE1, VOLe, XEQ, NA1, VOL, XCP, NRD)
         call SMTH(ALFA, NE1, ametre, XEQ, NA1, ametr, XCP, NRD)
         !        call SMTH(ALFA, NE1, FPe, XEQ, NA1, FP, XCP, NRD)
         call SMTH(ALFA, NE1, shife, XEQ, NA1, shif, XCP, NRD)
         call SMTH(ALFA, NE1, IPOLe, XEQ, NA1, IPOL, XCP, NRD)
         call SMTH(ALFA, NE1, G11e, XEQ, NA1, G11, XCP, NRD)
         call SMTH(ALFA, NE1, G33e, XEQ, NA1, G33, XCP, NRD)
         call SMTH(ALFA, NE1, SLATe, XEQ, NA1, SLAT, XCP, NRD)
         call SMTH(ALFA, n_xeq, BMINTe, XEQ, n_xcp,BMINT, XCP, NRD)
         call SMTH(ALFA, n_xeq, BMAXTe, XEQ, n_xcp,BMAXT, XCP, NRD)
         call SMTH(ALFA, n_xeq, BDB0e, XEQ, n_xcp, BDB0, XCP, NRD)
         call SMTH(ALFA, NE1, BDB02e, XEQ, NA1, BDB02, XCP, NRD)
         call SMTH(ALFA, NE1, B0DB2e, XEQ, NA1, B0DB2, XCP, NRD)
         call SMTH(ALFA, NE1, G22e, XEQ, NA1, G22, XCP, NRD)
      else
         IPOL  = IPOLe
         VOL = VOLE
         shif = shife
         AMETR = AMETRe
         G11 = G11e
         G33 = G33e
         G22 = G22e
         SLAT = SLATe
         BMINT = BMINTe
         BMAXT = BMAXTe
         BDB0 = BDB0e
         BDB02 = BDB02e
         B0DB2 = B0DB2e
         FP = FPe

      endif


!============================temporary ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      do j = 2, NA1
         VR(j) = (VOL(j) - VOL(j - 1))/HRO
         G11(j) = VR(j)*G11(j)
      end do
      VR(1) = VR(2)/3.
      G11(1) = VR(1)*G11(1)
!         write(*,*) 'VR',(VR(j),j=1,5)
      G22(1:NA1)= &
         VR(1:NA1)*G22(1:NA1)/IPOL(1:NA1)*RTOR/4./M_PI**2
      SQEPS(1:NA1)  = SQRT(AMETR(1:NA1)/(RTOR+SHIF(1:NA1)))
      if(SQEPS(1).le.1.d-2) SQEPS(1)=SQEPS(2)/1.4
 if(jprint.gt.0) then
         write(*,*) 'SQEPS',(SQEPS(j),j=1,5)
         write(*,*) 'AMETR',(AMETR(j),j=1,5)
         write(*,*) 'AMETRe',(AMETRe(j),j=1,5)
         write(*,*) 'G11',(G11(j),j=1,5)
         write(*,*) 'G22',(G22(j),j=1,5)
         write(*,*) 'G33',(G33(j),j=1,5)
         write(*,*) 'BMINT',(BMINT(j),j=1,5)
         write(*,*) 'BMAXT',(BMAXT(j),j=1,5)
         write(*,*) 'BDB0',(BDB0(j),j=1,5)
 endif

      do j=1,NA1
         y=min(.99d0,BTOR*RTOR/BMAXT(j)/(RTOR+SHIF(j)))
         FOFB(j) = B0DB2(j)*(1.-sqrt(1.d0-y)*(1.+0.5*y))
      enddo
 if(jprint.gt.0) then
         write(*,*) 'VOLe',(VOLE(j),j=1,5)
         write(*,*) 'VOL',(VOL(j),j=1,5)
         write(*,*) 'XCP',(XCP(j),j=1,5)
         write(*,*) 'XEQ',(XEQ(j),j=1,5)
         write(*,*) 'RTOR, AMETR(NA1),SHIF(1),SHIF(NA1)', RTOR, AMETR(NA1),SHIF(1),SHIF(NA1)
         write(*,*) 'q(0), q(na1), BTOR, VOL ',1./MU(1), 1./MU(NA1), BTOR, VOL(NA1)
         write(*,*) 'FP(1), FP(NA1)', FP(1), FP(NA1), FP(1)- FP(NA1)
         write(*,*) 'FPe(1), FPe(NE1)', FPe(1), FPe(NE1), FPe(1)- FPe(NE1)
 endif
      TAU = smart_in%TAU
      dtau = smart_in%dtau
! temporay vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv lines 768 -791 should be commented for external time loop control
!      open(20,file='tau.dat')
!      read(20,*) TAU
!      close(20)
   if(jprint.gt.0) then
         write(*,*) 'NA1', NA1,NE1
!         write(*,*) 'FP(1), FP(n_xcp)', FP(1), FP(n_xcp), FP(1)- FP(n_xcp)
         write(*,*) 'XCP',(XCP(j),j=1,NA1,30),XCP(NA1)
         write(*,*) 'RHO',(RHO(j),j=1,NA1,30),RHO(NA1)
         write(*,*) 'AMETR',(AMETR(j),j=1,NA1,30),AMETR(NA1)
         write(*,*) 'shif',(shif(j),j=1,NA1,30),shif(NA1)
         write(*,*) 'VOL',(VOL(j),j=1,NA1,30),VOL(NA1)
         write(*,*) 'VR',(VR(j),j=1,NA1,30),VR(NA1)
         write(*,*) 'SLAT',(SLAT(j),j=1,NA1,30),SLAT(NA1)
         write(*,*) 'SQEPS',(SQEPS(j),j=1,NA1,30),SQEPS(NA1)
         write(*,*) 'CU',(CU(j),j=1,NA1,30),CU(NA1)
         write(*,*) 'CUbs',(CUbs(j),j=1,NA1,30),CUbs(NA1)
         write(*,*) 'q',(1./MU(j),j=1,NA1,30),1./MU(NA1)
         write(*,*) 'Fp',(FP(j),j=1,NA1,30),FP(NA1)
         write(*,*) 'CC',(CC(j),j=1,NA1,30),CC(NA1)
         write(*,*) 'CD',(CD(j),j=1,NA1,30),CD(NA1)
         write(*,*) 'G33',(G33(j),j=1,NA1,30),G33(NA1)
         write(*,*) 'G22',(G22(j),j=1,NA1,30),G22(NA1)
         write(*,*) 'IPOL',(IPOL(j),j=1,NA1,30),IPOL(NA1)
         write(*,*) 'Ipl', IINTa(CU,ROC,RHO,G33,IPOL,NA1)
   endif
!=========================================================== time loop
      !open (1, file='out_Peltran.dat')
      TIME = TIMBEG
      !TIMPEL = 0.d0
      if(jprint.eq.1) then
      open(1,file='out_VAR.dat')
      write(1,997) VARNAME
      endif
997   format(30A14)
      if(jprint.lt.0) write(*,997) VARNAME
!if(jprint.gt.0) write(*,*) 'TE(1), TI(1)=', TE(1), TI(1)
!      open(20,file='jtime.dat')
!      read(20,*) jend
!      close(20)
!      open(20,file='jprint.dat')
!      read(20,*) jprint
!      close(20)
      jend=smart_in%key4control(17)
!      write(*,*) 'jprint,jend=',jprint,jend
      do jtime=1,jend   ! the line should be commented when external time loop is used
         !=============================================== OLDNEW
         TEo(1:NA1) = TE(1:NA1)
         TEx(1:NA1) = TE(1:NA1)
         TIo(1:NA1) = TI(1:NA1)
         TIx(1:NA1) = TI(1:NA1)
         NEo(1:NA1) = NE(1:NA1)
         NEx(1:NA1) = NE(1:NA1)
         NIo(1:NA1) = NI(1:NA1)
         NIx(1:NA1) = NI(1:NA1)
         F0o(1:NA1) = F0(1:NA1)
         F0x(1:NA1) = F0(1:NA1)
         F1o(1:NA1) = F1(1:NA1)
         F1x(1:NA1) = F1(1:NA1)
         F2o(1:NA1) = F2(1:NA1)
         F2x(1:NA1) = F2(1:NA1)
         F3o(1:NA1) = F3(1:NA1)
         F3x(1:NA1) = F3(1:NA1)
         F4o(1:NA1) = F4(1:NA1)
         F4x(1:NA1) = F4(1:NA1)
         F5o(1:NA1) = F5(1:NA1)
         F5x(1:NA1) = F5(1:NA1)
         FPo(1:NA1) = FP(1:NA1)
         VRo(1:NA1) = VR(1:NA1)
         do j=1,NA1
            pei(j) = 0.00246*(15.9 - .5*dlog(NE(j)) + dlog(TE(j))) &
               *NE(j)*Z2NdA(j)*(TE(J)-TI(J))/TE(j)/dsqrt(TE(j))
         enddo
         !================================================== auxilliary H&CD
         !=================================================== EC heating
   if(jprint.gt.0) write(*,*) 'Peecr1', VINTa(PEECR,ROC,RHO,VR,NA1)
         if (smart_in%sw_ech2a .ne. 0) then
            RHOEC  = smart_in%ROCEC*ROC ! EC location
            RHODR  = smart_in%ROCDR*ROC ! EC width
!
      QECR   = smart_in%QECR      ! QEC= 10 MW
            YEFFec = smart_in%YEFFec    ! IEC/QEC MA/MW

            !  ECH2a(YR0,YDR,YQ,YEFF,YP,YC,NA1,RHO,VR)
            if(smart_in%key4control(13).ne.0)        call ECH2a&
               (RHOEC, RHODR, QECR, YEFFec, PEECR, CUECR, NA1, RHO, VR, G33, IPOL, TE, NE)
!            YR0,  YDR,  YQ,   Y  EFF1,    YP,   YC,   NA1,  RHO,  VR,   G33,  IPOL, TE,NE
!         PE(1:NA1) = PEECR(1:NA1)
            if(jprint.gt.0)     write(*,*) 'Peecr2', VINTa(PEECR,ROC,RHO,VR,NA1)
   endif
!write(*,*) 'after ECH2, SQEPS', SQEPS(1)
         !=============================================== boundary conditions
         QNB = smart_in%QNB
         F01B = smart_in%F01B
         F02B = smart_in%F02B
         F03B = smart_in%F03B
         GN2E = smart_in%GN2E
         GN2I = smart_in%GN2I

         NEB = NE(NA1)
         TEB = TE(NA1)
         TIB = TI(NA1)
         F1B = F1(NA1)
         F2B = F2(NA1)
         F3B = F3(NA1)
         F4B = F4(NA1)
         F5B = F5(NA1)
         F6B = F6(NA1)
         F7B = F7(NA1)
         F8B = F8(NA1)
         F9B = F9(NA1)

         F0B = F01B + F02B + F03B
         !F0B is used only for normalization: puffing is controlled by QNB
         !fractions of neutral species: nH0B=  F01B/F0B, nD0B=  F02B/F0B, nT0B= F03B/F0B
   if(jprint.gt.0)               write(*,*) 'F0B=',F0B
         !========================================================transport coefficients
         !======================================== charged species
! to be replaced by external transport coefficients vvvvvvvvv
         do j = 1, NA1
            XI(j) = 0.5*(1.+3.*XCP(j)**2)
            HE(j) = XI(j)
            DN(j) = (HE(j) + XI(j))/10.
            CN(J) = 0.d0
            DF1(j) = (HE(j) + XI(j))/10.
            VF1(J) = 0.d0
            DF2(j) = (HE(j) + XI(j))/10.
            VF2(J) = 0.d0
            DF3(j) = (HE(j) + XI(j))/10.
            VF3(J) = 0.d0
            DF4(j) = (HE(j) + XI(j))/10.
            VF4(J) = 0.d0
            DF5(j) = (HE(j) + XI(j))/10.
            VF5(J) = 0.d0
            DF6(j) = (HE(j) + XI(j))/10.
            VF6(J) = 0.d0
            DF7(j) = (HE(j) + XI(j))/10.
            VF7(J) = 0.d0
            DF8(j) = (HE(j) + XI(j))/10.
            VF8(J) = 0.d0
            DF9(j) = (HE(j) + XI(j))/10.
            VF9(J) = 0.d0

         end do
! to be replaced by external transport coefficients ^^^^^^^^^

         ! for Pereverzev-Corrigan scheme
         do j=1,NA1
            DSE(J)=0.
            DSI(J)=0.
            DSN(J)=0.
         enddo
   if(jprint.gt.0)  then
            write(*,*) 'UPL',(UPL(j),j=1,NA1,30),UPL(NA1)
            write(*,*) 'Te',(TE(j),j=1,NA1,30),TE(NA1)
            write(*,*) 'Ti',(Ti(j),j=1,NA1,30),Ti(NA1)
            write(*,*) 'ne',(Ne(j),j=1,NA1,30),ne(NA1)
            write(*,*) 'AMAIN',(AMAIN(j),j=1,NA1,30),AMAIN(NA1)
   endif
         !======================================== for neutral transport
         do j = 1, NA1
            DF0(J) = 9.584d10*(TI(J) + 1.d-9)/(SVCXy(Ti(j),AMAIN(j)) + 1.d-10)/AMAIN(J)/NE(j)

            if (j .lt. NA) then
               VF0(J) = -7.6d-1*DF0(J)*(TI(J + 1) - TI(J))/HRO/(TI(J) + 1.d-9)
            else
               VF0(J) = -7.6d-1*DF0(J)*(TI(NA1) - TI(NA))/HRO/(TI(J) + 1.d-9)
            end if

            SF0(J) = (F1(J) + F2(J) + F3(J))*NE(J)*SVRCy(TE(J))

            SFF0(J) = -NEo(j)*SVIEy(TE(J))

         end do
         !        write(*,*) 'a1=',TEo(1:NA1)
         !==============================================density stepup
!
         !write(*,*) 'SF0(J), SFF0(J)=', SF0(1), SFF0(1)
! if(jprint.gt.0) then
!          write(*,*) 'DF0',(DF0(j),j=1,NA1,30),DF0(NA1)
!         write(*,*) 'VF0',(VF0(j),j=1,NA1,30),VF0(NA1)
!         write(*,*) 'SF0',(SF0(j),j=1,NA1,30),SF0(NA1)
!         write(*,*) 'SFF0',(SFF0(j),j=1,NA1,30),SFF0(NA1)
!         write(*,*) 'F0',(F0(j),j=1,NA1,30),F0(NA1)
!endif
         !======================================= neutrals
         if(smart_in%key4control(10).ne.0)      call STEPUPN0(&
            NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO,&
            VF0, DF0, F0o, F0, F0X, QNB, SF0, SFF0, SF0TOT, QF0, GF0, GF0X)
! write(*,*) 'VF0, DF0 F0(0) F0(a)=',VF0(1), DF0(1), F0(1), F0(NA1)
!if(jprint.eq.2) then
!          write(*,*) 'F0',(F0(j),j=1,NA1,30),F0(NA1)

!         write(*,*) 'SF0TOT', VINTa(SF0TOT,ROC,RHO,VR,NA1)
!endif
         !=================================================
         do j = 1, NA1
! particle sources/sinks  due to fusion and neutrals
            Y   = f2(J)*f3(J)*SVDBHy(TI(j))             ! n(14MeV)+4He(3.52MeV)
            Y1  = F2(j)*F2(j)*SVD1y(TI(j))/2.   ! n2.45(MeV)+3He(0.817MeV)
            Y2  = F2(j)*F2(j)*SVD2y(TI(j))/2.   ! t(1MeV)+p(3MeV)
            Y3  = F2(j)*stbrny(TE(j),NE(j),Z2NdA(j)) ! probability of burn out t(1MeV+dth=> 4He+n14MeV
            SF1(j)      = NEo(J)*(SVIEy(TE(J))*F0(J)*F01B/F0B-SVRCy(TE(J))*F1o(J)) + Y2
            SF2(j)      = NEo(J)*(SVIEy(TE(J))*F0(J)*F02B/F0B-SVRCy(TE(J))*F2o(J)) - Y - Y1 - Y2*Y3
            SF3(j)      = NEo(J)*(SVIEy(TE(J))*F0(J)*F03B/F0B-SVRCy(TE(J))*F3o(J)) -Y + Y2*(1.-Y3)
            SF4(j)      = Y + Y2*Y3
            SF5(j)      = Y1
            Sn14(j)     = Y + Y2*Y3
            Sn245(j)= Y1

! distributions of fusion products
            PDT         = 5.632*Sn14(j) ! 3.52*1.6
            PDD1        = 1.3072*Y1             ! 0.817*1.6
            PDD2        = 1.6128*Y2             !1.008*1.6
            PDD3        = 4.8*Y2                !3.*1.6

            Y4  =       (NE(j)/Z2NdA(j))**.66667
            Y   =       60.27*Y4                        ! 3520/14.6/4
            Y1  =       18.653*Y4                       ! 817./14.6/3.
            Y2  =       23.01*Y4                        ! 1008/14.6/3
            Y3  =       205.48*Y4                       ! 3000/14.6/1
!               TAUSp=2.d0*yABEAM/yNEJ*yTEJ*DSQRT(yTEJ)/YLE slowing down of proton
            TAUSp       =2.d0/NE(J)*TE(J)*DSQRT(TE(J))/(15.85d0+DLOG(TE(J)/DSQRT(NE(J))))
! nfast = S*tauSp*Afast/Zfast**2*ln(1+(Vb/Vc)**3)/3 fast ion desnity
            F1fast(j) = SF1(j)*TAUSp*DLOG(1.+(Y3/TE(J))**1.5)/3.          !fast protons
            F3fast(j) = SF3(j)*TAUSp*3.*DLOG(1.+(Y2/TE(J))**1.5)/3.       !fast t
            F4fast(j) = SF4(j)*TAUSp*DLOG(1.+(Y/TE(J))**1.5)/3.       !fast 4He
            F5fast(j) = SF5(j)*TAUSp*0.75*DLOG(1.+(Y1/TE(J))**1.5)/3.  !fast 3He
            NHYDR(j)  = F1(j) +F1fast(j)
            NDEUT(j)  = F2(j) +F2fast(j)
            NTRIT(j)  = F3(j) + F3fast(j)
            NALF(j)   = F4(j) + F4fast(j)
            NHE3(j)   = F5(j) + F5fast(j)
! fast ion pressure [keV*10^19/m3]
            PFAST(j)=0. &
               +PDT*(1.-PAIONy(TE(J),y))*625.*TAUSp/3. &
               +PDD1*(1.-PAIONy(TE(J),y1))*625.*TAUSp*0.25 &
               +PDD2*(1.-PAIONy(TE(J),y2))*625.*TAUSp &
               +PDD3*(1.-PAIONy(TE(J),y3))*625.*TAUSp/3.
! heating by fusion products
            pifus(j)  =0. &
               +PDT*PAIONy(TE(J),y) &
               +PDD1*PAIONy(TE(J),y1) &
               +PDD2*PAIONy(TE(J),y2) &
               +PDD3*PAIONy(TE(J),y3)
            pefus(j)=0. &
               +PDT*(1.-PAIONy(TE(J),y)) &
               +PDD1*(1.-PAIONy(TE(J),y1)) &
               +PDD2*(1.-PAIONy(TE(J),y2)) &
               +PDD3*(1.-PAIONy(TE(J),y3))
! temparature of neutrals in the diffisive model
            TN(j)  =Ti(j)+ &
               (Te(j)-Ti(j))*195.*PEI(j)/SVCXy(Ti(j),AMAIN(j))/(ni(j)+f0(j))/ni(j)
! heat loss due to neutrals
            PEN(j) =-SVRCy(TE(J))*TE(j)*1.5/625.*NE(j)*(f1(j)+f2(j)+f3(j)) &
               -(4.d-2*SVIEy(TE(J))/625.+PENLIy(TE(J)))*f0(j)*NE(J)
            PIN(J) =1.5/625.*( &
               +(SVCXy(Ti(j),AMAIN(j))*f0(j)*(TN(j)-TI(j)) &
               -SVRCy(TE(J))*TI(j))*(f1(j)+f2(j)+f3(j)))*NE(j)

         enddo
!write(*,*) 'Pfast, Pefus, PEN, PIN ',Pfast(1), Pefus(1), PEN(1), PIN(1)
!========================================= Bootstrap current
!write(*,*) 'AMAIN, ZMAIN, ZEF, Z2NdA ',AMAIN(1), ZMAIN(1), ZEF(1), Z2NdA(1)
!write(*,*) 'TI, TE, MU, NE, NI ',TI(1), TE(1), MU(1), NE(1), NI(1)
!goto 111
!      open(20,file='jpsi.dat')
!      read(20,*) jpsi
!      close(20)
!      if(jpsi.ne.0) then
         call CUBSy( &
            NA1, RTOR, BTOR, IPL, &
            FP, MU, ZEF, TE, TI, NE, NI, AMAIN, ZMAIN, &
            BMINT, BMAXT, BDB0, BDB02, FOFB, SQEPS, RHO, &
            CUBS, CC)   ! output: bootsrap current density and curent conductivity by Sauter
!write(*,*) 'CUBS, CC ', CUBS(1), CC(1)
!      cubs(1:na1) =0.
         call RHSEQy( &
            NA1, RTOR, BTOR, RHO, NE, NI, TE, TI, PBLON, PBPER, PFAST, &
            MU, CU, G22, G33, IPOL, AMETR, &
            CUTOR, EQFF, EQPF)  !out: toroidal current density, RHS for equilibrium equation

!write(*,*) 'CUTOR, EQFF, EQPF',  CUTOR(1), EQFF(1), EQPF(1)
!============================================ current diffusion (+equilibrium)
         if(smart_in%key4control(11).gt.0)     call STEPUPF( &
            NA1, RHO, TAU, RTOR, BTOR, IPL, CUBS, CD, CC, G22, G33, IPOL, &
            FP, FPo, MU, CU, UPL, ULON, FV)
!         endif ! psi
! ======================================= Ohmic heating
! PJOUL=CUTOR(J)*UPL(J)/(M_PI2*RTOR)
!               PJOUL(1:NA1)=CUTOR(1:NA1)*UPL(1:NA1)/(M_PI2*RTOR)
! POH [MW/m#3]: Power of Ohmic Heating
!       P=sigma*Ez**2
!               (Pereverzev 12-FEB-90)
         do j=1,NA1
            PJOUL(j) =CC(j)*(ULON(j)/(2.*M_PI*RTOR*IPOL(j)))**2/G33(j)
         enddo
! 111 continue
!write(*,*) 'PJOUL, CC, ULON ',PJOUL(1), CC(1), ULON(1)
!=======================================================
         PEAUX(1:NA1) = PEBM(1:NA1)+PEECR(1:NA1)+PEICR(1:NA1)
         PIAUX(1:NA1) = PIBM(1:NA1)+PIICR(1:NA1)
!============================= total heat sources w/o equipartition
         PE(1:NA1) = PEAUX(1:NA1) +PJOUL(1:NA1) +PEFUS(1:NA1) +PEN(1:NA1) -PRAD(1:NA1)
         PI(1:NA1) = PIAUX(1:NA1) +PIFUS(1:NA1) +PIN(1:NA1)
         CD(1:NA1) = CUECR(1:NA1)
   if(jprint.gt.0) then
            write(*,*) 'Peecr', VINTa(PEECR,ROC,RHO,VR,NA1)
            write(*,*) 'Peaux', VINTa(PEAUX,ROC,RHO,VR,NA1)
            write(*,*) 'Piaux', VINTa(PiAUX,ROC,RHO,VR,NA1)
            write(*,*) 'Pjoul', VINTa(Pjoul,ROC,RHO,VR,NA1)
            write(*,*) 'Pefus', VINTa(PEFUS,ROC,RHO,VR,NA1)
            write(*,*) 'Pifus', VINTa(Pifus,ROC,RHO,VR,NA1)
            write(*,*) 'CD', IINTa(CD,ROC,RHO,G33,IPOL,NA1)
            write(*,*) 'CUBS', IINTa(CUBS,ROC,RHO,G33,IPOL,NA1)
   endif
!     open(20,file='jdens.dat')
!      read(20,*) jdens
!      close(20)
!      if(jdens.ne.0) then
         !======================================= hydrogen species
!      if (iH .ne. 0)
         if(smart_in%key4control(1).gt.0) call STEPUPN( &
            NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
            VF1, DF1, DSN, F1o, F1, F1X, F1B, SF1, SFF1, SF1TOT, QF1, GF1, GF1X)
!      if (iD .ne. 0)
         if(smart_in%key4control(2).gt.0) call STEPUPN( &
            NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
            VF2, DF2, DSN, F2o, F2, F2X, F2B, SF2, SFF2, SF2TOT, QF2, GF2, GF2X)
!      if (iT .ne. 0)
         if(smart_in%key4control(3).gt.0) call STEPUPN( &
            NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
            VF3, DF3, DSN, F3o, F3, F3X, F3B, SF3, SFF3, SF3TOT, QF3, GF3, GF3X)


         if(smart_in%key4control(4).gt.0) call STEPUPN( &
            NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
            VF4, DF4, DSN, F4o, F4, F4X, F4B, SF4, SFF4, SF4TOT, QF4, GF4, GF4X)
!      if (iT .ne. 0)
         if(smart_in%key4control(5).gt.0) call STEPUPN( &
            NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
            VF5, DF5, DSN, F5o, F5, F5X, F5B, SF5, SFF5, SF5TOT, QF5, GF5, GF5X)

         if(smart_in%key4control(6).gt.0) call STEPUPN( &
            NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
            VF6, DF6, DSN, F6o, F6, F6X, F6B, SF6, SFF6, SF6TOT, QF6, GF6, GF6X)

         if(smart_in%key4control(7).gt.0) call STEPUPN( &
            NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
            VF7, DF7, DSN, F7o, F7, F7X, F7B, SF7, SFF7, SF7TOT, QF7, GF7, GF7X)
         if(smart_in%key4control(8).gt.0) call STEPUPN( &
            NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
            VF8, DF8, DSN, F8o, F8, F8X, F8B, SF8, SFF8, SF8TOT, QF8, GF8, GF8X)

         if(smart_in%key4control(9).gt.0) call STEPUPN( &
            NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
            VF9, DF9, DSN, F9o, F9, F9X, F9B, SF9, SFF9, SF9TOT, QF9, GF9, GF9X)
         !============================ electron density from quasineutrality
         ne(1:na1) = f1(1:na1) + f2(1:na1) + f3(1:na1) + 2.*(f4(1:na1) + f5(1:na1)) &
            + f1fast(1:na1) + f2fast(1:na1) + f3fast(1:na1) &
            + 2.*(f4fast(1:na1) + f5fast(1:na1))
         ni(1:na1) = f1(1:na1) + f2(1:na1) + f3(1:na1) + f4(1:na1) + f5(1:na1) ! themal ion density

!======
         if (ispec(4) .gt. 0) then
            do j = 4, n_ion + 3 - Khydr
               i = ispec(j)
               ne(1:na1) = ne(1:na1) + &
                  cp_in%profiles_1d(i_time)%ion(i)%density(1:na1)/denA2D* &
                  cp_in%profiles_1d(i_time)%ion(i)%z_ion_1D(1:na1)
               ni(1:na1) = ni(1:na1) + &
                  cp_in%profiles_1d(i_time)%ion(i)%density(1:na1)/denA2D
            end do
         end if
         !        write(*,*) 'QF0,GF0,QFB',QF0(NA1),GF0(NA1),QNB
         !        write(*,*) 'QF1,GF1',QF1(NA1),GF1(NA1)
         !        write(*,*) 'QF2,GF2',QF2(NA1),GF2(NA1)
         !        write(*,*) 'QF3,GF3',QF3(NA1),GF3(NA1)
         !        write(*,*) 'QF3(vint)',VINTa(SF3TOT,ROC,RHO,VR,NA1)
         !==============================================temperature stepup
         !        write(*,*) 'trace 325 before STEPUPT'
!      endif !density

         !     open(20,file='jheat.dat')
         !     read(20,*) jheat
         !     close(20)
         !     if(jheat.ne.0) &
         if(smart_in%key4control(12).ne.0) call stepupt(&
            NA1, NB1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
            XI, HE, DSI, DSE, &
            PE, PET, PETOT, PI, PIT, PITOT, Z2NdA, &
            TEX, TE, TEo, TEB, TIX, TI, TIo, TIB, &
            NEX, NEo, NE, NIX, NIo, NI, Qe, Qi, GNX, GN2E, GN2I &
            )
         !================================================================
         if (smart_in%sw_stdout .ne.0) then
            if(jprint.gt.0)     write (*, 100) 'QE, QI, Ge     = ', QE(NA1), QI(NA1), QF1(NA1)+QF2(NA1)+QF3(NA1)
         endif

         !============================================= pelshot
         !        Write(*,*) 'before pellet'
         !        time=time + TAU
         TIMPEL = TIMPEL + TAU
         YDABL = 0.d0
         YDDEP = 0.d0
      if (TIMPEL .ge. (dtau - 1.d-7)) then
!         TIMPEL = 0.d0
         !== Pellet Ablation Model: SMART
!         if (smart_in%sw_smart .ne. 0) then
!      open(20,file='jpel.dat')
!      read(20,*) jpel
!      close(20)
         if(smart_in%key4control(15).ne.0) then
            call pelIMAS1(smart_in%YAM, smart_in%YVP, smart_in%YVOL, &
               smart_in%YCOS0, smart_in%YEFF, smart_in%YDL, &
               YDABL, YDDEP, YPELSRS, smart_in%yswitch, &
               ne, ni, Te, Ti, F1, F2, F3, FP, &
               ametr, shif, vr, mu, &
               HRO, ROC, BTOR, RTOR, NA1, NRD)
           end if
            TIMPEL = 0.d0
      end if

         !        Write(*,*) 'after pellet'
         !======================= calculation of the delay between shoot and ablation
         !        write(*,*) 'ne-neo',(ne(1:NA1)-neo(1:NA1))
         !        write(*,*) 'Te-Teo',(Te(1:NA1)-Teo(1:NA1))
         jabs = idint((1.0d0-YDABL)*NA1)
         if (smart_in%YCOS0 .gt. 0.) then
            Dtpel = (smart_in%yglength + smart_in%YDL*(ametr(NA1) - ametr(jabs) + shif(jabs)))/smart_in%YVP/1000.
         else
            Dtpel = (smart_in%yglength + smart_in%YDL*(ametr(NA1) - ametr(jabs) - shif(jabs)))/smart_in%YVP/1000.
         end if
         !        write(*,*) 'Dtpel =',Dtpel
         !============================================ current diffusion (+equilibrium)
         !============================================ end of itterations

         if (smart_in%sw_stdout .ne.0) then
      if(jprint.gt.0) then
               write(*, 100)'Te(1),  Ti(1)  = ',Te(1),Ti(1)
               write(*, 100)'ne(1),  ni(1)  = ',ne(1),ni(1)
               write(*, 100)'n0(1),  n0(a)  = ',F0(1),F0(NA1)
               write(*, 100)'<ne>,   QF0B   = ',VINTa(ne, ROC, RHO, VR, NA1)/VOL(NA1),QNB
               write(*, 100)'<Shdt>, <Sn0>  = ',VINTa(SF3TOT, ROC, RHO, VR, NA1) + &
                  VINTa(SF2TOT, ROC, RHO, VR, NA1) + &
                  VINTa(SF1TOT, ROC, RHO, VR, NA1),  &
                  VINTa(SF0TOT, ROC, RHO, VR, NA1)
      endif
         end if
!      if(jmix.ne.0) &
!         if(smart_in%key4control(14).ne.0.and.time.gt.277.d0) &
        if(smart_in%key4control(14).ne.0) &
            call MIXF19(1.4d0,1.d0,NA1,RHO,VR,AMAIN,ZEF,G33,G22,IPOL, &
            F1,F2,F3,F4,F5,F6,F7,F8,F9,TE,TI,NE,NI,MU,PFAST,PBLON,PBPER, &
            CU,CUTOR,EQPF,EQFF,NHYDR,NDEUT,NTRIT,NALF,NHE3,FP, &
            RTOR,BTOR,SHIFT,IMIX)

100      format(A17,5(1PE15.6))
200      format(A17,10i5)

         TIME = TIME + TAU
         if(jprint.gt.0) then
            write (*, 100) 'time   = ', time
            write (*, 100) 'Pec,Pe,Pi,cc0,J0 = ', VINTa(PEECR, ROC, RHO, VR, NA1),&
               VINTa(PE, ROC, RHO, VR, NA1),&
               VINTa(PI, ROC, RHO, VR, NA1), CC(1), CU(1)
            write (*, 100) 'Pei, POH, FP(a) = ', VINTa(Pei, ROC, RHO, VR, NA1),&
               VINTa(PJOUL, ROC, RHO, VR, NA1),FP(NA1)
            write (*, 100) '<ne>, q(0), q(a) = ', VINTa(NE, ROC, RHO, VR, NA1)/VOL(NA1),1./mu(1),1./mu(NA1)
            write (*, 100) 'Ibs, Itot, U(0), U(a) =', IINTa(CUBS,ROC,RHO,G33,IPOL,NA1), &
               IINTa(CU,ROC,RHO,G33,IPOL,NA1),UPL(1),UPL(NA1), &
               IINTa(CD,ROC,RHO,G33,IPOL,NA1)
         endif
   write(1,998) &
   time,Te(1),Ti(1),ne(1),ni(1),Te(NA1),Ti(NA1),ne(NA1),ni(NA1), &
            F0(1),F0(NA1),VINTa(NE, ROC, RHO, VR, NA1)/VOL(NA1),IINTa(CUBS,ROC,RHO,G33,IPOL,NA1),&
            IINTa(CU, ROC, RHO, G33, IPOL, NA1), IINTa(CD, ROC, RHO, G33, IPOL, NA1),ULON(1),ULON(NA1),FP(1),FP(NA1), &
            VINTa(Pe, ROC, RHO, VR, NA1), VINTa(Pi, ROC, RHO, VR, NA1), VINTa(Pefus, ROC, RHO, VR, NA1), &
            VINTa(Pifus, ROC, RHO, VR, NA1), VINTa(PJOUL, ROC, RHO, VR, NA1), VINTa(PEECR, ROC, RHO, VR, NA1), &
            VINTa(Pei, ROC, RHO, VR, NA1), 1./mu(1), 1./mu(na1), VINTa(SF3TOT, ROC, RHO, VR, NA1) + &
            VINTa(SF2TOT, ROC, RHO, VR, NA1) + &
            VINTa(SF1TOT, ROC, RHO, VR, NA1),  &
            VINTa(SF0TOT, ROC, RHO, VR, NA1)
   if(jprint.lt.0)    write(*,998) &
   time,Te(1),Ti(1),ne(1),ni(1),Te(NA1),Ti(NA1),ne(NA1),ni(NA1), &
            F0(1),F0(NA1),VINTa(NE, ROC, RHO, VR, NA1)/VOL(NA1),IINTa(CUBS,ROC,RHO,G33,IPOL,NA1),&
            IINTa(CU, ROC, RHO, G33, IPOL, NA1), IINTa(CD, ROC, RHO, G33, IPOL, NA1),ULON(1),ULON(NA1),FP(1),FP(NA1), &
            VINTa(Pe, ROC, RHO, VR, NA1), VINTa(Pi, ROC, RHO, VR, NA1), VINTa(Pefus, ROC, RHO, VR, NA1), &
            VINTa(Pifus, ROC, RHO, VR, NA1), VINTa(PJOUL, ROC, RHO, VR, NA1), VINTa(PEECR, ROC, RHO, VR, NA1), &
            VINTa(Pei, ROC, RHO, VR, NA1), 1./mu(1), 1./mu(na1), VINTa(SF3TOT, ROC, RHO, VR, NA1) + &
            VINTa(SF2TOT, ROC, RHO, VR, NA1) + &
            VINTa(SF1TOT, ROC, RHO, VR, NA1),  &
            VINTa(SF0TOT, ROC, RHO, VR, NA1)
!         write(1,*) 'time', time
      enddo      ! end of time loop . The line should de commented for external time control
      if(jprint.eq.1) then
      close(1)
      write(*,*) 'jprint',jprint
      endif
998   format(30(1XPE13.6))
   if(jprint.gt.0) then
         write(*,*) 'UPL',(UPL(j),j=1,NA1,30),UPL(NA1)
         write(*,*) 'CU',(CU(j),j=1,NA1,30),CU(NA1)
         write(*,*) 'CUbs',(CUbs(j),j=1,NA1,30),CUbs(NA1)
         write(*,*) 'q',(1./MU(j),j=1,NA1,30),1./MU(NA1)
         write(*,*) 'Fp',(FP(j),j=1,NA1,30),FP(NA1)
         write(*,*) 'CC',(CC(j),j=1,NA1,30),CC(NA1)
         write(*,*) 'CD',(CD(j),j=1,NA1,30),CD(NA1)
         write(*,*) 'CU(1-5)',(CU(j),j=1,5)
         write(*,*) 'q(1-5)',(1./mu(j),j=1,5)
         write(*,*) 'UPL(1-5)',(UPL(j),j=1,5)
         write(*,*) 'Ipl', IINTa(CU,ROC,RHO,G33,IPOL,NA1)
! temporary for ASTRA vvvvvvvvvvvvvv
         write(*,*) 'ABC', time, abc
         write(*,*) 'BTOR',time,BTOR
         write(*,*) 'RTOR',time,rtor
         write(*,*)  'AB', time, abc*1.1
         write(*,*) 'ELON',time, elon
         write(*,*) 'TRIAN',time,trian
         write(*,*) 'ELONG',time, elon
         write(*,*) 'TRICH',time,trian
         write(*,*) 'IPL',time,IPL
         write(*,*) 'AMJ',time,AMAIN(1)
         write(*,*) 'ZMJ',time,ZMAIN(1)
   endif
!      include 'Out4ASTRA.inc'
! temporayr for ASTRA ^^^^^^^^^^^^^^
!         write(*,*) 'G33',(G33(j),j=1,NA1,30),G33(NA1)
!         write(*,*) 'G22',(G22(j),j=1,NA1,30),G22(NA1)
!         write(*,*) 'IPOL',(IPOL(j),j=1,NA1,30),IPOL(NA1)
      !========================================================
      !== conversion to IMAS units
      cp_out%time(i_time) = TIME
      if(NABEG.eq.1) then
         cp_out%profiles_1d(i_time)%electrons%density(1:n_xcp) = ne(1:n_xcp)*denA2D
         cp_out%profiles_1d(i_time)%electrons%temperature(1:n_xcp) = Te(1:n_xcp)*temA2D
         cp_out%profiles_1d(i_time)%q(1:n_xcp) = 1./ MU(1:n_xcp)
         cp_out%profiles_1d(i_time)%grid%psi(1:n_xcp) =sign_psi*FP(1:n_xcp)
         cp_out%profiles_1d(i_time)%j_total(1:n_xcp)=sign_psi*CU(1:n_xcp)*cuA2D
         cp_out%profiles_1d(i_time)%j_bootstrap(1:n_xcp)=sign_psi*CUBS(1:n_xcp)*cuA2D
         cp_out%profiles_1d(i_time)%j_non_inductive(NABEG:n_xcp)=sign_psi*CD(1:n_xcp)*cuA2D
         if (ispec(1) .ne. 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(1))%density(1:n_xcp) = F1(1:n_xcp)*denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(1))%density(1:n_xcp) = F01B/F0B*F0(1:n_xcp)*denA2D
         end if
         if (ispec(2) .ne. 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(2))%density(1:n_xcp) = F2(1:n_xcp)*denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(2))%density(1:n_xcp) = F02B/F0B*F0(1:n_xcp)*denA2D
         end if
         if (ispec(3) .ne. 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(3))%density(1:n_xcp) = F3(1:n_xcp)*denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(3))%density(1:n_xcp) = F03B/F0B*F0(1:n_xcp)*denA2D
         end if
         do i = 1, n_ion
            cp_out%profiles_1d(i_time)%ion(i)%temperature(1:n_xcp) = Ti(1:n_xcp)*temA2D
         end do
      else
         cp_out%profiles_1d(i_time)%electrons%density(NABEG:n_xcp) = ne(1:NA1)*denA2D
         cp_out%profiles_1d(i_time)%electrons%temperature(NABEG:n_xcp) = Te(1:NA1)*temA2D
         cp_out%profiles_1d(i_time)%q(NABEG:n_xcp) = 1./ MU(1:NA1)
         cp_out%profiles_1d(i_time)%grid%psi(NABEG:n_xcp) =sign_psi*FP(1:NA1)
         cp_out%profiles_1d(i_time)%j_total(NABEG:n_xcp)=sign_psi*CU(1:NA1)*cuA2D
         cp_out%profiles_1d(i_time)%j_bootstrap(NABEG:n_xcp)=sign_psi*CUBS(1:NA1)*cuA2D
         cp_out%profiles_1d(i_time)%j_non_inductive(NABEG:n_xcp)=sign_psi*CD(1:NA1)*cuA2D
         if (ispec(1) .ne. 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(1))%density(NABEG:n_xcp) = F1(1:NA1)*denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(1))%density(NABEG:n_xcp) = F01B/F0B*F0(1:NA1)*denA2D
         end if
         if (ispec(2) .ne. 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(2))%density(NABEG:n_xcp) = F2(1:NA1)*denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(2))%density(NABEG:n_xcp) = F02B/F0B*F0(1:NA1)*denA2D
         end if
         if (ispec(3) .ne. 0) then
            cp_out%profiles_1d(i_time)%ion(ispec(3))%density(NEBEG:n_xcp) = F3(1:NA1)*denA2D
            cp_out%profiles_1d(i_time)%neutral(ispec(3))%density(NABEG:n_xcp) = F03B/F0B*F0(1:NA1)*denA2D
         end if
         do i = 1, n_ion
            cp_out%profiles_1d(i_time)%ion(i)%temperature(NABEG:n_xcp) = Ti(1:NA1)*temA2D
         end do
         do j=1,NABEG
            cp_out%profiles_1d(i_time)%electrons%density(j) = ne(1)*denA2D
            cp_out%profiles_1d(i_time)%electrons%temperature(j) = Te(1)*temA2D
            cp_out%profiles_1d(i_time)%q(j) = 1./ MU(1)
            cp_out%profiles_1d(i_time)%grid%psi(j) =sign_psi*FP(1)
            cp_out%profiles_1d(i_time)%j_total(j)=sign_psi*CU(1)*cuA2D
            cp_out%profiles_1d(i_time)%j_bootstrap(j)=sign_psi*CUBS(1)*cuA2D
            cp_out%profiles_1d(i_time)%j_non_inductive(j)=sign_psi*CD(1)*cuA2D
            if (ispec(1) .ne. 0) then
               cp_out%profiles_1d(i_time)%ion(ispec(1))%density(j) = F1(1)*denA2D
               cp_out%profiles_1d(i_time)%neutral(ispec(1))%density(j) = F01B/F0B*F0(1)*denA2D
            end if
            if (ispec(2) .ne. 0) then
               cp_out%profiles_1d(i_time)%ion(ispec(2))%density(j) = F2(1)*denA2D
               cp_out%profiles_1d(i_time)%neutral(ispec(2))%density(j) = F02B/F0B*F0(1)*denA2D
            end if
            if (ispec(3) .ne. 0) then
               cp_out%profiles_1d(i_time)%ion(ispec(3))%density(j) = F3(1)*denA2D
               cp_out%profiles_1d(i_time)%neutral(ispec(3))%density(j) = F03B/F0B*F0(1)*denA2D
            end if
            do i = 1, n_ion
               cp_out%profiles_1d(i_time)%ion(i)%temperature(j) = Ti(1)*temA2D
            end do
         enddo
      endif
      !        write(*,*) 'vr, n_xcp',n_xcp, vr(1:n_xcp)
      !        write(*,*) 'vole, n_xcp',n_xcp, vr(1:n_xcp)
      !        write(*,*) 'shif, n_xcp',n_xcp, shif(1:n_xcp)
      if (smart_in%sw_stdout .ne.0) then
         if(jprint.gt.0) then
            write (*, 100) 'YDABL, YDDEP   = ', YDABL, YDDEP
            write (*, 100) 'RHOEC, RHODR   = ', RHOEC, RHODR
            write (*, 100) 'ROC, QECR      = ', ROC, QECR
            write (*, 100) 'YEFFec         = ', YEFFec
            write (*, 100) 'Pecr, Pe, Pi   = ', VINTa(PEECR, ROC, RHO, VR, NA1),&
               VINTa(PEECR, ROC, RHO, VR, NA1),&
               VINTa(PI, ROC, RHO, VR, NA1)
         endif
      endif
      deallocate (ispec)
      !        include 'dealloc.corprf
      deallocate (ne, ni, Te, Ti, nex, nix, TEX, TIX, TN, NN,&
         F0, F1, F2, F3, F4, F5, F6, F7, F8, F9, &
         F0x, F1x, F2x, F3x, F4x, F5x, F6x, F7x, F8x, F9x, VTOR, FP,&
         neo, nio, Teo, Tio, F0o, F1o, F2o, F3o, F4o, F5o, F6o, F7o, F8o, F9o, FPo,&
         cu, cutor, cd, cubs, UPL, ULON, EZ, ZEF, AMAIN, Z2NdA, ZMAIN, &
         NHYDR, NDEUT, NTRIT, NALF, NHE3 &
         )
      !        include 'dealloc.corsrs'
      deallocate (SN, SNN, SNTOT, QN, GN, GNX, PE, PET, PETOT,&
         SF0, SFF0, SF0TOT, QF0, GF0, GF0X,&
         SF1, SFF1, SF1TOT, QF1, GF1, GF1X,&
         SF2, SFF2, SF2TOT, QF2, GF2, GF2X,&
         SF3, SFF3, SF3TOT, QF3, GF3, GF3X,&
         SF4, SFF4, SF4TOT, QF4, GF4, GF4X,&
         SF5, SFF5, SF5TOT, QF5, GF5, GF5X,&
         SF6, SFF6, SF6TOT, QF6, GF6, GF6X,&
         SF7, SFF7, SF7TOT, QF7, GF7, GF7X,&
         SF8, SFF8, SF8TOT, QF8, GF8, GF8X,&
         SF9, SFF9, SF9TOT, QF9, GF9, GF9X,&
         PI, PIT, PITOT, PEI, QE, QI, PEECR, CUECR, YPELSRS, &
         PEFUS, PIFUS, PEAUX, PIAUX, PEN, PIN, PJOUL, PRAD, &
         PEBM, PIBM, PEICR, PIICR, &
         Sn14, Sn245, SCUBM )
      !        include 'dealloc.cortran'
      deallocate (DF0, VF0, DF1, VF1, DF2, VF2, DF3, VF3, DF4, VF4, DF5, VF5, &
         DF6, VF6, DF7, VF7, DF8, VF8, DF9, VF9, &
         DN, CN, HE, XI, CC, DSE, DSI, DSN)
      !        write(*,*) '308'
      !        include 'dealloc.eq'
      deallocate (IPOL, G11, G33, SLAT, ametr, shif, vr, vro, SQEPS, &
         mu, VOL, XCP, ametre, shife, VOLe,&
         FPe, XEQ, RHO, IPOLe, G11e, G33e, SLATe, &
         BMINT, BMAXT, BDB0, BDB02, B0DB2, FOFB, G22, EQFF, EQPF, FV, &
         BMINTe, BMAXTe, BDB0e, BDB02e, B0DB2e, FOFBe, G22e)
      !        write(*,*) '313'

   end subroutine smart

end module mod_smart

!======================================================================|
