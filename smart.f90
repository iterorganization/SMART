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

      integer :: i, j, i_time, j_time, n_xcp, n_xeq, n_ion, nrd, NA1, jtime

      double precision &
         denA2D, temA2D, presA2D, cuA2D, &
         HRO, ROC, BTOR, SHIFT, ABC, RTOR, ALFA, &
         YDABL, YDDEP, ai, zi, Dtpel, &
         NEB, TEB, TIB, F0B, F1B, F2B, F3B, F4B, F5B, F01B, F02B, F03B, &
         GN2E, GN2I, SVRCy, SVCXy, SVIEy, RHOEC, RHODR, QECR, YEFFec, &
         y, y1, y2, y3, y4, pdt, PAIONy, VINTa, QNB, IINTa, IPL, &
         SVD1y, SVD2y, SVDBHy, stbrny, PENLIy, PDD1, PDD2, PDD3, TAUSp

      !
      double precision &
         TAU, dtau, TIME, TIMBEG
      double precision, save :: TIMPEL = 0.d0

      integer, allocatable :: ispec(:)
      !
      integer:: iH, iD, iT, na, nb1, Nhydr, JABS, istep

      double precision, allocatable :: &
         ne(:), ni(:), Te(:), Ti(:), TN(:), F0(:), F1(:), F2(:), F3(:), F4(:), F5(:), FP(:), &
         neo(:), nio(:), Teo(:), Tio(:), F0o(:), F1o(:), F2o(:), F3o(:), F4o(:), F5o(:), FPo(:), &
         neX(:), niX(:), TeX(:), TiX(:), F0x(:), F1x(:), F2x(:), F3x(:), F4x(:), F5x(:), NN(:), &
         cu(:), cutor(:), cd(:), cubs(:), UPL(:), ULON(:), EZ(:), FV(:), &
         ZEF(:), AMAIN(:), Z2NdA(:), ZMAIN(:), SQEPS(:), PBLON(:), PBPER(:), PFAST(:), &
         F1fast(:),  F2fast(:),  F3fast(:),  F4fast(:), F5fast(:)

      double precision, allocatable :: &
         SN(:), SNN(:), SNTOT(:), QN(:), GN(:), GNX(:), QE(:), QI(:), &
         SF0(:), SFF0(:), SF0TOT(:), QF0(:), GF0(:), GF0X(:), &
         SF1(:), SFF1(:), SF1TOT(:), QF1(:), GF1(:), GF1X(:), &
         SF2(:), SFF2(:), SF2TOT(:), QF2(:), GF2(:), GF2X(:), &
         SF3(:), SFF3(:), SF3TOT(:), QF3(:), GF3(:), GF3X(:), &
         SF4(:), SFF4(:), SF4TOT(:), QF4(:), GF4(:), GF4X(:), &
         SF5(:), SFF5(:), SF5TOT(:), QF5(:), GF5(:), GF5X(:), &
         PE(:), PET(:), PETOT(:), PI(:), PIT(:), PITOT(:), PEI(:), &
         PEECR(:), CUECR(:), YPELSRS(:), &
         PEFUS(:), PIFUS(:), PEAUX(:), PIAUX(:), PEN(:), PIN(:), PJOUL(:), PRAD(:), &
         Sn14(:), Sn245(:)

      double precision, allocatable :: &
         DF0(:), VF0(:), DF1(:), VF1(:), DF2(:), VF2(:), DF3(:), VF3(:), &
         DF4(:), VF4(:), DF5(:), VF5(:), &
         DN(:), CN(:), HE(:), XI(:), cc(:), DSI(:), DSE(:), DSN(:)

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
      external SMTH, ECH2a, STEPUPN, STEPUPN0, STEPUPT, STEPUPF, CUBSy, RHSEQy, ARR, PelIMAS1

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
      NA1 = n_xcp
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
         F1(n_xcp), F2(n_xcp), F3(n_xcp), F4(n_xcp), F5(n_xcp), FP(n_xcp), FV(n_xcp),&
         neo(n_xcp), nio(n_xcp), Teo(n_xcp), Tio(n_xcp), &
         F1o(n_xcp), F2o(n_xcp), F3o(n_xcp), F4o(n_xcp), F5o(n_xcp), FPo(n_xcp), &
         F0(n_xcp), F0o(n_xcp), F0x(n_xcp), NN(n_xcp), TN(n_xcp), &
         cu(n_xcp), cutor(n_xcp), cd(n_xcp), cubs(n_xcp), &
         UPL(n_xcp), ULON(n_xcp), EZ(n_xcp), ZEF(n_xcp), AMAIN(n_xcp), Z2NdA(n_xcp), ZMAIN(n_xcp), &
         PBLON(n_xcp), PBPER(n_xcp), PFAST(n_xcp), &
         F1fast(n_xcp),  F2fast(n_xcp),  F3fast(n_xcp),  F4fast(n_xcp), F5fast(n_xcp))
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
      F1 = 0.d0
      F2 = 0.d0
      F3 = 0.d0
      F4 = 0.d0
      F5 = 0.d0
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
      FPo = 0.d0
      F0 = 0.d0
      F0o = 0.d0
      F0x = 0.d0
      NN = 0.d0
      TN = 0.d0
      cu = 0.d0
      cutor = 0.d0
      cd = 0.d0
      cubs = 0.d0
      UPL = 0.d0
      ULON = 0.d0
      EZ = 0.d0
      ZEF = 0.d0
      AMAIN = 0.d0
      Z2NdA = 0.d0
      FV=0.d0
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
         DF0(n_xcp), VF0(n_xcp), DF1(n_xcp), VF1(n_xcp), &
         DF2(n_xcp), VF2(n_xcp), DF3(n_xcp), VF3(n_xcp), &
         DF4(n_xcp), VF4(n_xcp), DF5(n_xcp), VF5(n_xcp), &
         PE(n_xcp), PET(n_xcp), PETOT(n_xcp), QE(n_xcp), QI(n_xcp), &
         PI(n_xcp), PIT(n_xcp), PITOT(n_xcp), PEI(n_xcp), &
         PEFUS(n_xcp), PIFUS(n_xcp), PEAUX(n_xcp), PIAUX(n_xcp), &
         PEN(n_xcp), PIN(n_xcp), PJOUL(n_xcp), PRAD(n_xcp), &
         Sn14(n_xcp), Sn245(n_xcp), &
         YPELSRS(n_xcp) )

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
      DF0 = 0.d0
      VF0 = 0.d0
      DF1 = 0.d0
      VF1 = 0.d0
      DF2 = 0.d0
      VF2 = 0.d0
      DF3 = 0.d0
      VF3 = 0.d0
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

      !====================================================== core profiles
      ne(1:n_xcp) = cp_in%profiles_1d(i_time)%electrons%density(1:n_xcp)/denA2D
!      Te(1:n_xcp) = cp_in%profiles_1d(i_time)%electrons%temperature(1:n_xcp)/temA2D
      Ti(1:n_xcp) = cp_in%profiles_1d(i_time)%ion(1)%temperature(1:n_xcp)/temA2D
      Te(1:n_xcp) = Ti(1:n_xcp)
      write(*,*), 'TE(1)', TE(1)
!cp_in%profiles_1d(i_time)%electrons%temperature(1:n_xcp)/temA2D
      MU(1:n_xcp) = 1./cp_in%profiles_1d(i_time)%q(1:n_xcp)
 write(*,*) 'J0, CC0 DB', cp_in%profiles_1d(i_time)%j_total(1)/cuA2D, &
      cp_in%profiles_1d(i_time)%conductivity_parallel(1)/cuA2D
      !============================== detect hydrogen isotopes
      F0(1:n_xcp) = 0.d0
      F1(1:n_xcp) = 0.d0
      F2(1:n_xcp) = 0.d0
      F3(1:n_xcp) = 0.d0
      ispec(1:n_ion + 3) = 0
      Nhydr = 0
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
               Nhydr = Nhydr + 1
            end if
            if (ai .gt. 1.d0 .and. ai .lt. 3.d0) then
               ispec(2) = i
               Nhydr = Nhydr + 1
            end if
            if (ai .gt. 2.5d0) then
               ispec(3) = i
               Nhydr = Nhydr + 1
            end if
         else
            istep = istep + 1
            ispec(3 + istep) = i
         end if
      end do
      if (Nhydr .lt. 1) then
         write (*, *) 'No hydrogen isotopes in the ion list'
         write (*, *) 'Hyrogen species: iH,iD,iT Nhydr', iH, iD, iT, Nhydr
         stop
      end if
      if (smart_in%sw_stdout .ne.0) then
         write (*, 200) 'ispec          = ', ispec
      end if
      if (ispec(1) .ne. 0) then
         iH = ispec(1)
         F1(1:n_xcp) = cp_in%profiles_1d(i_time)%ion(iH)%density(1:n_xcp)/denA2D
         F0(1:n_xcp) = F0(1:n_xcp) + cp_in%profiles_1d(i_time)%neutral(iH)%density(1:n_xcp)/denA2D
      end if
      if (ispec(2) .ne. 0) then
         iD = ispec(2)
         F2(1:n_xcp) = cp_in%profiles_1d(i_time)%ion(iD)%density(1:n_xcp)/denA2D
         F0(1:n_xcp) = F0(1:n_xcp) + cp_in%profiles_1d(i_time)%neutral(iD)%density(1:n_xcp)/denA2D
      end if
      if (ispec(3) .ne. 0) then
         iT = ispec(3)
         F3(1:n_xcp) = cp_in%profiles_1d(i_time)%ion(iT)%density(1:n_xcp)/denA2D
         F0(1:n_xcp) = F0(1:n_xcp) + cp_in%profiles_1d(i_time)%neutral(iT)%density(1:n_xcp)/denA2D
      end if
      AMAIN(1:NA1) = (F1(1:NA1) + 2.*F2(1:NA1) + 3.*F3(1:NA1))/ &
                   & (F1(1:NA1) + F2(1:NA1) + F3(1:NA1))
    !!!!!!!!!!!!!!!!!!!!!!!!!! separate ions, av mass, zeff, pei
      do j = 1, n_xcp
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
      end do

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
         XEQ(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(1:n_xeq)
         VOLe(1:n_xeq) = dmax1(0., eq_in%time_slice(j_time)%profiles_1d%VOLUME(1:n_xeq))
         XCP(1:n_xcp) = max(1.d-19,cp_in%profiles_1d(i_time)%grid%rho_tor_norm(1:n_xcp))
         ROC = eq_in%time_slice(j_time)%profiles_1d%phi(n_xeq)
         ROC = dsqrt(dabs(ROC/BTOR/M_PI))
         RHO(1:n_xcp) = ROC*XCP(1:n_xcp)

         ametre(1:n_xeq) = dmax1(0., (eq_in%time_slice(j_time)%profiles_1d%r_outboard(1:n_xeq) - &
                              & eq_in%time_slice(j_time)%profiles_1d%r_inboard(1:n_xeq))/2.)
         shife(1:n_xeq) = (eq_in%time_slice(j_time)%profiles_1d%r_outboard(1:n_xeq) + &
                      & eq_in%time_slice(j_time)%profiles_1d%r_inboard(1:n_xeq))/2.-RTOR
!         FPe(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%psi(1:n_xeq)
!
FPe(1:n_xeq) = -eq_in%time_slice(j_time)%profiles_1d%psi(1:n_xeq) + eq_in%time_slice(j_time)%profiles_1d%psi(1)
         IPOLe(1:n_xeq) = dabs(eq_in%time_slice(j_time)%profiles_1d%f(1:n_xeq)/RTOR/BTOR)
         SLATe(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%surface(1:n_xeq)
         G11e(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%gm3(1:n_xeq)
         G33e(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%gm1(1:n_xeq)*RTOR**2

!      BMINTe(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%b_field_min(1:n_xeq)
!      BMAXTe(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%b_field_max(1:n_xeq)
!      BDB0e(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%b_field_average(1:n_xeq)
         BDB02e(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%gm5(1:n_xeq)/BTOR**2
         B0DB2e(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%gm4(1:n_xeq)*BTOR**2
         G22e(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%gm2(1:n_xeq)
!   write(*,*) 'BDB02e, B0DB2e, G22e, G11e, G33e, SLATe, IPOLe'
!   write(*,*) BDB02e(1), B0DB2e(1), G22e(1), G11e(1), G33e(1), SLATe(1), IPOLe(1)
      !write (*, *) 'Ip =',
         IPL = dabs(eq_in%time_slice(j_time)%global_quantities%Ip/cuA2D)
         ABC = ametre(n_xeq)
         SHIFT = Shife(n_xeq) ! - RTOR
         HRO = (XCP(3) - XCP(2))*ROC
         ALFA = 0.001d0
      !SMTH(ALFA,NO,FO,XO,N,FN,XN) mapping from EQ to CP grids
      !        write(*,*) 'trace 184'
!      if (n_xeq .ne. n_xcp) then
         call SMTH(ALFA, n_xeq, VOLe, XEQ, n_xcp, VOL, XCP, NRD)
         call SMTH(ALFA, n_xeq, ametre, XEQ, n_xcp, ametr, XCP, NRD)
         call SMTH(ALFA, n_xeq, FPe, XEQ, n_xcp, FP, XCP, NRD)
         call SMTH(ALFA, n_xeq, shife, XEQ, n_xcp, shif, XCP, NRD)
         call SMTH(ALFA, n_xeq, IPOLe, XEQ, n_xcp, IPOL, XCP, NRD)
         call SMTH(ALFA, n_xeq, G11e, XEQ, n_xcp, G11, XCP, NRD)
         call SMTH(ALFA, n_xeq, G33e, XEQ, n_xcp, G33, XCP, NRD)
         call SMTH(ALFA, n_xeq, SLATe, XEQ, n_xcp, SLAT, XCP, NRD)
!         call SMTH(ALFA, n_xeq, BMINTe, XEQ, n_xcp,BMINT, XCP, NRD)
!         call SMTH(ALFA, n_xeq, BMAXTe, XEQ, n_xcp,BMAXT, XCP, NRD)
!         call SMTH(ALFA, n_xeq, BDB0e, XEQ, n_xcp, BDB0, XCP, NRD)
         call SMTH(ALFA, n_xeq, BDB02e, XEQ, n_xcp, BDB02, XCP, NRD)
         call SMTH(ALFA, n_xeq, B0DB2e, XEQ, n_xcp, B0DB2, XCP, NRD)
         call SMTH(ALFA, n_xeq, G22e, XEQ, n_xcp, G22, XCP, NRD)
         VR(1) = VOL(2)/HRO
         G11(1) = VR(1)*G11(1)
      do j = 2, n_xcp
         VR(j) = (VOL(j) - VOL(j - 1))/HRO
         G11(j) = VR(j)*G11(j)
      end do
		G22(1:n_xcp)= &
     	VR(1:n_xcp)*G22(1:n_xcp)/IPOL(1:n_xcp)*RTOR/4./M_PI**2

!		SHIF(1:n_xcp)=SHIF(1:n_xcp)-SHIFT
		SQEPS(1) = SQRT(AMETR(2)*0.5/(RTOR+SHIF(1)))
		SQEPS(2:n_xcp)  = SQRT(AMETR(2:n_xcp)/(RTOR+SHIF(2:n_xcp)))

      do j=1,n_xcp
         BMAXT(j) = BTOR*RTOR/(RTOR+SHIF(j)-AMETR(j)) ! temporary
         y=min(.99d0,BTOR*RTOR/BMAXT(j)/(RTOR+SHIF(j)))
		FOFB(j) = B0DB2(j)*(1.-sqrt(1.d0-y)*(1.+0.5*y))
      enddo
          write(*,*) 'BMAXT(1), BMAXT(NA1)',BMAXT(1), BMAXT(NA1)
      write(*,*) 'FOFB(1), FOFB(NA1)', FOFB(1), FOFB(NA1)
       write(*,*) 'G22(1), G22(NA1)', G22(1), G22(NA1)
      write(*,*) 'G33(1), G33(NA1)', G33(1), G33(NA1)
      write(*,*) 'G11(1), G11(NA1)', G11(1), G11(NA1)
         write(*,*) 'IPOL(1), IPOL(NA1)', IPOL(1), IPOL(NA1)
         write(*,*) 'RTOR, AMETR(NA1),SHIF(1),SHIF(NA1)', RTOR, AMETR(NA1),SHIF(1),SHIF(NA1)
         write(*,*) 'q(0), q(na1), BTOR, VOL ',1./MU(1), 1./MU(NA1), BTOR, VOL(NA1)
!         FP(1:n_xcp) = -(cp_in%profiles_1d(i_time)%grid%psi(1:n_xcp) &
!                  - cp_in%profiles_1d(i_time)%grid%psi(1))

!   write(*,*) 'BDB02, B0DB2, G22, G11, G33, SLAT, IPOL'
!   write(*,*) BDB02(1), B0DB2e(1), G22(1), G11(1), G33(1), SLAT(1), IPOL(1)
   write(*,*) 'FP(1), FP(n_xcp)', FP(1), FP(n_xcp), FP(1)- FP(n_xcp)
   write(*,*) 'FPe(1), FPe(n_xeq)', FPe(1), FPe(n_xeq), FPe(1)- FPe(n_xeq)
      TAU = smart_in%TAU
      dtau = smart_in%dtau
      write(*,*) 'NA1', NA1,n_xcp
         write(*,*) 'FP(1), FP(n_xcp)', FP(1), FP(n_xcp), FP(1)- FP(n_xcp)
         write(*,*) 'CU',(CU(j),j=1,NA1,30)
         write(*,*) 'CUbs',(CUbs(j),j=1,NA1,30)
         write(*,*) 'q',(1./MU(j),j=1,NA1,30)
         write(*,*) 'Fp',(FP(j),j=1,NA1,30)
         write(*,*) 'CC',(CC(j),j=1,NA1,30)
      !=========================================================== time loop
      !open (1, file='out_Peltran.dat')
      TIME = TIMBEG
      !TIMPEL = 0.d0
      open(1,file='out_VAR.dat')
      write(1,997) VARNAME
 997  format(30A14)
 write(*,*) 'TE(1), TI(1)=', TE(1), TI(1)
      do jtime=1,20
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
      if (smart_in%sw_ech2a .ne. 0) then
         RHOEC  = smart_in%ROCEC*ROC ! EC location
         RHODR  = smart_in%ROCDR*ROC ! EC width
         QECR   = 0.d0
         !smart_in%QECR      ! QEC= 10 MW
         YEFFec = smart_in%YEFFec    ! IEC/QEC MA/MW

      !  ECH2a(YR0,YDR,YQ,YEFF,YP,YC,NA1,RHO,VR)
         call ECH2a&
            (RHOEC, RHODR, QECR, YEFFec, PEECR, CUECR, NA1, RHO, VR, G33, IPOL, TE, NE)
!         PE(1:NA1) = PEECR(1:NA1)
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

      F0B = F01B + F02B + F03B
      !F0B is used only for normalization: puffing is controlled by QNB
      !fractions of neutral species: nH0B=  F01B/F0B, nD0B=  F02B/F0B, nT0B= F03B/F0B
      !        write(*,*) 'F0B=',F0B
      !========================================================transport coefficients
      !======================================== charged species
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
      end do

      ! for Pereverzev-Corrigan scheme
      do j=1,NA1
         DSE(J)=0.
         DSI(J)=0.
         DSN(J)=0.
      enddo

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
!write(*,*) 'SF0(J), SFF0(J)=', SF0(1), SFF0(1)
      !======================================= neutrals
      call STEPUPN0(&
         NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO,&
         VF0, DF0, F0o, F0, F0X, QNB, SF0, SFF0, SF0TOT, QF0, GF0, GF0X)
!write(*,*) 'VF0, DF0 TI1 TI2=',VF0(2), DF0(1), TI(1), TI(2)
      !=================================================
      do j = 1, NA1
! particle sources/sinks  due to fusion and neutrals
		Y	= f2(J)*f3(J)*SVDBHy(TI(j))		! n(14MeV)+4He(3.52MeV)
		Y1	= F2(j)*F2(j)*SVD1y(TI(j))/2.	! n2.45(MeV)+3He(0.817MeV)
		Y2	= F2(j)*F2(j)*SVD2y(TI(j))/2.	! t(1MeV)+p(3MeV)
		Y3	= F2(j)*stbrny(TE(j),NE(j),Z2NdA(j)) ! probability of burn out t(1MeV+dth=> 4He+n14MeV
		SF1(j)	= NEo(J)*(SVIEy(TE(J))*F0(J)*F01B/F0B-SVRCy(TE(J))*F1o(J)) + Y2
		SF2(j)	= NEo(J)*(SVIEy(TE(J))*F0(J)*F02B/F0B-SVRCy(TE(J))*F2o(J)) - Y - Y1 - Y2*Y3
		SF3(j)	= NEo(J)*(SVIEy(TE(J))*F0(J)*F03B/F0B-SVRCy(TE(J))*F3o(J)) -Y + Y2*(1.-Y3)
		SF4(j)	= Y + Y2*Y3
		SF5(j)	= Y1
		Sn14(j)	= Y + Y2*Y3
		Sn245(j)= Y1

! distributions of fusion products
			PDT 	= 5.632*Sn14(j)	! 3.52*1.6
			PDD1 	= 1.3072*Y1		! 0.817*1.6
			PDD2 	= 1.6128*Y2		!1.008*1.6
			PDD3	= 4.8*Y2		!3.*1.6

			Y4	=	(NE(j)/Z2NdA(j))**.66667
			Y	=	60.27*Y4			! 3520/14.6/4
			Y1	=	18.653*Y4			! 817./14.6/3.
			Y2	= 	23.01*Y4			! 1008/14.6/3
			Y3	= 	205.48*Y4			! 3000/14.6/1
!		TAUSp=2.d0*yABEAM/yNEJ*yTEJ*DSQRT(yTEJ)/YLE slowing down of proton
		TAUSp	=2.d0/NE(J)*TE(J)*DSQRT(TE(J))/(15.85d0+DLOG(TE(J)/DSQRT(NE(J))))
! nfast = S*tauSp*Afast/Zfast**2*ln(1+(Vb/Vc)**3)/3 fast ion desnity
		F1fast(j) = SF1(j)*TAUSp*DLOG(1.+(Y3/TE(J))**1.5)/3.	  !fast protons
		F3fast(j) = SF3(j)*TAUSp*3.*DLOG(1.+(Y2/TE(J))**1.5)/3.	  !fast t
		F4fast(j) = SF4(j)*TAUSp*DLOG(1.+(Y/TE(J))**1.5)/3.	      !fast 4He
		F5fast(j) = SF5(j)*TAUSp*0.75*DLOG(1.+(Y1/TE(J))**1.5)/3.  !fast 3He
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
	call CUBSy( &
      NA1, RTOR, BTOR, IPL, &
      FP, MU, ZEF, TE, TI, NE, NI, AMAIN, ZMAIN, &
      BMINT, BMAXT, BDB0, BDB02, FOFB, SQEPS, RHO, &
      CUBS, CC)	! output: bootsrap current density and curent conductivity by Sauter
!write(*,*) 'CUBS, CC ', CUBS(1), CC(1)
      cubs(1:na1) =0.
	call RHSEQy( &
      NA1, RTOR, BTOR, RHO, NE, NI, TE, TI, PBLON, PBPER, PFAST, &
      MU, CU, G22, G33, IPOL, AMETR, &
      CUTOR, EQFF, EQPF)	!out: toroidal current density, RHS for equilibrium equation

!write(*,*) 'CUTOR, EQFF, EQPF',  CUTOR(1), EQFF(1), EQPF(1)
!============================================ current diffusion (+equilibrium)
      call STEPUPF( &
         NA1, RHO, TAU, RTOR, BTOR, IPL, CUBS, CD, CC, G22, G33, IPOL, &
         FP, FPo, MU, CU, UPL, ULON, FV)
! ======================================= Ohmic heating
! PJOUL=CUTOR(J)*UPL(J)/(M_PI2*RTOR)
!		PJOUL(1:NA1)=CUTOR(1:NA1)*UPL(1:NA1)/(M_PI2*RTOR)
! POH [MW/m#3]:	Power of Ohmic Heating
!	P=sigma*Ez**2
!		(Pereverzev 12-FEB-90)
    do j=1,NA1
      PJOUL(j) =CC(j)*(ULON(j)/(2.*M_PI*RTOR*IPOL(j)))**2/G33(j)
	enddo
! 111 continue
!write(*,*) 'PJOUL, CC, ULON ',PJOUL(1), CC(1), ULON(1)
!=======================================================
		PEAUX(1:NA1) = PEECR(1:NA1)
!============================= total heat sources w/o equipartition
		PE(1:NA1) = PEAUX(1:NA1) +PJOUL(1:NA1) +PEFUS(1:NA1) +PEN(1:NA1) -PRAD(1:NA1)
		PI(1:NA1) = PIAUX(1:NA1) +PIFUS(1:NA1) +PIN(1:NA1)
		CD(1:NA1) = CUECR(1:NA1)

      !======================================= hydrogen species
      if (iH .ne. 0) call STEPUPN( &
         NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
         VF1, DF1, DSN, F1o, F1, F1X, F1B, SF1, SFF1, SF1TOT, QF1, GF1, GF1X)
      if (iD .ne. 0) call STEPUPN( &
         NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
         VF2, DF2, DSN, F2o, F2, F2X, F2B, SF2, SFF2, SF2TOT, QF2, GF2, GF2X)
      if (iT .ne. 0) call STEPUPN( &
         NA1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
         VF3, DF3, DSN, F3o, F3, F3X, F3B, SF3, SFF3, SF3TOT, QF3, GF3, GF3X)
      !============================ electron density from quasineutrality
         ne(1:na1) = f1(1:na1) + f2(1:na1) + f3(1:na1) + 2.*(f4(1:na1) + f5(1:na1)) &
               + f1fast(1:na1) + f2fast(1:na1) + f3fast(1:na1) &
               + 2.*(f4fast(1:na1) + f5fast(1:na1))
         ni(1:na1) = f1(1:na1) + f2(1:na1) + f3(1:na1) + f4(1:na1) + f5(1:na1) ! themal ion density

!======
      if (ispec(4) .gt. 0) then
         do j = 4, n_ion + 3 - nhydr
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

      call stepupt(&
          NA1, NB1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
      XI, HE, DSI, DSE, &
      PE, PET, PETOT, PI, PIT, PITOT, Z2NdA, &
      TEX, TE, TEo, TEB, TIX, TI, TIo, TIB, &
      NEX, NEo, NE, NIX, NIo, NI, Qe, Qi, GNX, GN2E, GN2I &
      )
      !================================================================
      if (smart_in%sw_stdout .ne.0) then
         write (*, 100) 'QE, QI, Ge     = ', QE(NA1), QI(NA1), QF1(NA1)+QF2(NA1)+QF3(NA1)
      endif
 
      !============================================= pelshot
      !        Write(*,*) 'before pellet'
      !        time=time + TAU
         TIMPEL = TIMPEL + TAU
         YDABL = 0.d0
         YDDEP = 0.d0
      if (TIMPEL .ge. (dtau - 1.d-7)) then
         TIMPEL = 0.d0
         !== Pellet Ablation Model: SMART
         if (smart_in%sw_smart .ne. 0) then
            call pelIMAS1(smart_in%YAM, smart_in%YVP, smart_in%YVOL, &
                          smart_in%YCOS0, smart_in%YEFF, smart_in%YDL, &
                          YDABL, YDDEP, YPELSRS, smart_in%yswitch, &
                          ne, ni, Te, Ti, F1, F2, F3, FP, &
                          ametr, shif, vr, mu, &
                          HRO, ROC, BTOR, RTOR, NA1, NRD)
         end if
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
         write(*, 100)'Te(1),  Ti(1)  = ',Te(1),Ti(1)
         write(*, 100)'ne(1),  ni(1)  = ',ne(1),ni(1)
         write(*, 100)'n0(1),  n0(a)  = ',F0(1),F0(NA1)
         write(*, 100)'<ne>,   QF0B   = ',VINTa(ne, ROC, RHO, VR, NA1)/VOLe(n_xeq),QNB
         write(*, 100)'<Shdt>, <Sn0>  = ',VINTa(SF3TOT, ROC, RHO, VR, NA1) + &
                                          VINTa(SF2TOT, ROC, RHO, VR, NA1) + &
                                          VINTa(SF1TOT, ROC, RHO, VR, NA1),  &
                                          VINTa(SF0TOT, ROC, RHO, VR, NA1)
      end if

 100  format(A17,5(1PE15.6))
 200  format(A17,10i5)

      TIME = TIME + TAU
               write (*, 100) 'time   = ', time
         write (*, 100) 'Pec,Pe,Pi,cc0,J0 = ', VINTa(PEECR, ROC, RHO, VR, NA1),&
                                             VINTa(PE, ROC, RHO, VR, NA1),&
                                             VINTa(PI, ROC, RHO, VR, NA1), CC(1), CU(1)
         write (*, 100) 'Pei, POH, FP(a) = ', VINTa(Pei, ROC, RHO, VR, NA1),&
                                             VINTa(PJOUL, ROC, RHO, VR, NA1),FP(NA1)
         write (*, 100) '<ne>, q(0), q(a) = ', VINTa(NE, ROC, RHO, VR, NA1)/VOL(NA1),1./mu(1),1./mu(NA1)
         write (*, 100) 'Ibs, Itot, U||(a),Icd =', IINTa(CUBS,ROC,RHO,G33,IPOL,NA1), &
                                             IINTa(CU,ROC,RHO,G33,IPOL,NA1),ULON(NA1), &
                                             IINTa(CD,ROC,RHO,G33,IPOL,NA1)
 write(1,998) time,Te(1),Ti(1),ne(1),ni(1),Te(NA1),Ti(NA1),ne(NA1),ni(NA1), &
      F0(1),F0(NA1),VINTa(NE, ROC, RHO, VR, NA1)/VOL(NA1),IINTa(CUBS,ROC,RHO,G33,IPOL,NA1),&
      IINTa(CU, ROC, RHO, G33, IPOL, NA1), IINTa(CD, ROC, RHO, G33, IPOL, NA1),ULON(1),ULON(NA1),FP(1),FP(NA1), &
      VINTa(Pe, ROC, RHO, VR, NA1), VINTa(Pi, ROC, RHO, VR, NA1), VINTa(Pefus, ROC, RHO, VR, NA1), &
      VINTa(Pifus, ROC, RHO, VR, NA1), VINTa(PJOUL, ROC, RHO, VR, NA1), VINTa(PEECR, ROC, RHO, VR, NA1), &
      VINTa(Pei, ROC, RHO, VR, NA1), 1./mu(1), 1./mu(na1), VINTa(SF3TOT, ROC, RHO, VR, NA1) + &
                                          VINTa(SF2TOT, ROC, RHO, VR, NA1) + &
                                          VINTa(SF1TOT, ROC, RHO, VR, NA1),  &
                                          VINTa(SF0TOT, ROC, RHO, VR, NA1)
 enddo
      ! end of time loop
 	close(1)
 998	format(30(1XPE13.6))
         write(*,*) 'FP(1), FP(n_xcp)', FP(1), FP(n_xcp), FP(1)- FP(n_xcp)
         write(*,*) 'CU',(CU(j),j=1,NA1,30)
         write(*,*) 'CUbs',(CUbs(j),j=1,NA1,30)
         write(*,*) 'CD',(CD(j),j=1,NA1,30)
         write(*,*) 'q',(1./MU(j),j=1,NA1,30)
         write(*,*) 'Fp',(FP(j),j=1,NA1,30)
         write(*,*) 'CC',(CC(j),j=1,NA1,30)
      !========================================================
      !== conversion to IMAS units
      cp_out%time(i_time) = TIME

      cp_out%profiles_1d(i_time)%electrons%density(1:n_xcp) = ne(1:n_xcp)*denA2D
      cp_out%profiles_1d(i_time)%electrons%temperature(1:n_xcp) = Te(1:n_xcp)*temA2D
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
      !        write(*,*) 'vr, n_xcp',n_xcp, vr(1:n_xcp)
      !        write(*,*) 'vole, n_xcp',n_xcp, vr(1:n_xcp)
      !        write(*,*) 'shif, n_xcp',n_xcp, shif(1:n_xcp)
      if (smart_in%sw_stdout .ne.0) then
         write (*, 100) 'YDABL, YDDEP   = ', YDABL, YDDEP
         write (*, 100) 'RHOEC, RHODR   = ', RHOEC, RHODR
         write (*, 100) 'ROC, QECR      = ', ROC, QECR
         write (*, 100) 'YEFFec         = ', YEFFec
         write (*, 100) 'Pecr, Pe, Pi   = ', VINTa(PEECR, ROC, RHO, VR, NA1),&
                                             VINTa(PEECR, ROC, RHO, VR, NA1),&
                                             VINTa(PI, ROC, RHO, VR, NA1)
      endif

!,&

      deallocate (ispec)
      !        include 'dealloc.corprf
      deallocate (ne, ni, Te, Ti, nex, nix, TEX, TIX, TN, NN,&
                  F0, F1, F2, F3, F4, F5, F0x, F1x, F2x, F3x, F4x, F5x, FP,&
                  neo, nio, Teo, Tio, F0o, F1o, F2o, F3o, F4o, F5o, FPo,&
                  cu, cutor, cd, cubs, UPL, ULON, EZ, ZEF, AMAIN, Z2NdA, ZMAIN &
                  )
      !        include 'dealloc.corsrs'
      deallocate (SN, SNN, SNTOT, QN, GN, GNX, PE, PET, PETOT,&
                  SF0, SFF0, SF0TOT, QF0, GF0, GF0X,&
                  SF1, SFF1, SF1TOT, QF1, GF1, GF1X,&
                  SF2, SFF2, SF2TOT, QF2, GF2, GF2X,&
                  SF3, SFF3, SF3TOT, QF3, GF3, GF3X,&
                  SF4, SFF4, SF4TOT, QF4, GF4, GF4X,&
                  SF5, SFF5, SF5TOT, QF5, GF5, GF5X,&
                  PI, PIT, PITOT, PEI, QE, QI, PEECR, CUECR, YPELSRS, &
                  PEFUS, PIFUS, PEAUX, PIAUX, PEN, PIN, PJOUL, PRAD, &
                  Sn14, Sn245 )
      !        include 'dealloc.cortran'
      deallocate (DF0, VF0, DF1, VF1, DF2, VF2, DF3, VF3, DF4, VF4, DF5, VF5, &
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
