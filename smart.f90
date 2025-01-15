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

      integer :: i, j, i_time, j_time, n_xcp, n_xeq, n_ion, nrd, NA1
      double precision &
         denA2D, temA2D, presA2D, cuA2D, &
         HRO, ROC, BTOR, SHIFT, ABC, RTOR, ALFA, &
         YDABL, YDDEP, ai, zi, Z2NdA, Dtpel, &
         NEB, TEB, TIB, F0B, F1B, F2B, F3B, F01B, F02B, F03B, &
         GN2E, GN2I, SVRC, SVCX, SVIE, RHOEC, RHODR, QECR, YEFFec, &
         y, y2, pdt, svdt, paion2, VINTa, QNB
      !
      double precision &
         TAU, dtau, TIME, TIMBEG
      double precision, save :: TIMPEL = 0.0

      integer, allocatable :: ispec(:)
      !
      integer:: iH, iD, iT, na, nb1, Nhydr, JABS, istep
      !
      double precision, allocatable :: &
         ne(:), ni(:), Te(:), Ti(:), TN(:), F0(:), F1(:), F2(:), F3(:), FP(:), &
         neo(:), nio(:), Teo(:), Tio(:), F0o(:), F1o(:), F2o(:), F3o(:), FPo(:), &
         neX(:), niX(:), TeX(:), TiX(:), F0x(:), F1x(:), F2x(:), F3x(:), NN(:), &
         cu(:), cutor(:), cd(:), cubs(:), UPL(:), UPOL(:), EZ(:), &
         ZEF(:), AMAIN(:)

      double precision, allocatable :: &
         SN(:), SNN(:), SNTOT(:), QN(:), GN(:), GNX(:), QE(:), QI(:), &
         SF0(:), SFF0(:), SF0TOT(:), QF0(:), GF0(:), GF0X(:), &
         SF1(:), SFF1(:), SF1TOT(:), QF1(:), GF1(:), GF1X(:), &
         SF2(:), SFF2(:), SF2TOT(:), QF2(:), GF2(:), GF2X(:), &
         SF3(:), SFF3(:), SF3TOT(:), QF3(:), GF3(:), GF3X(:), &
         PE(:), PET(:), PETOT(:), PI(:), PIT(:), PITOT(:), PEI(:), &
         PECR(:), CUECR(:), YPELSRS(:)

      double precision, allocatable :: &
         DF0(:), VF0(:), DF1(:), VF1(:), DF2(:), VF2(:), DF3(:), VF3(:), &
         DN(:), CN(:), HE(:), XI(:), cc(:)

      double precision, allocatable :: &
         ametr(:), shif(:), mu(:), VOL(:), vr(:), VOLo(:), vro(:), &
         IPOL(:), G11(:), G33(:), SLAT(:)
      double precision, allocatable :: &
         ametre(:), shife(:), VOLe(:), XCP(:), XEQ(:), RHO(:), FPe(:), &
         IPOLe(:), G11e(:), G33e(:), SLATe(:)

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
            write (*, *) 'Input pellets IDS NOT detected'
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
      write (*, *) 'i_time,time=', i_time, TIMBEG
      !              i_time=size(cp_in%profiles_1d(:)%time)
      j_time = size(eq_in%time_slice)
      n_xcp = size(cp_in%profiles_1d(i_time)%grid%rho_tor_norm)
      n_ion = size(cp_in%profiles_1d(i_time)%ion)
      n_xeq = size(eq_in%time_slice(j_time)%profiles_1d%psi)
      nrd = 2*n_xcp
      NA1 = n_xcp
      NB1 = NA1
      NA = NA1 - 1
      write (*, *) 'NA1,nrd,n_xeq', NA1, nrd, n_xeq
      !==========================
      allocate (ispec(n_ion + 3))
      !==========================
      !        include 'alloc.corprf'
      allocate ( &
         ne(n_xcp), ni(n_xcp), Te(n_xcp), Ti(n_xcp), &
         nex(n_xcp), nix(n_xcp), TEX(n_xcp), TIX(n_xcp), &
         F1x(n_xcp), F2x(n_xcp), F3x(n_xcp), &
         F1(n_xcp), F2(n_xcp), F3(n_xcp), FP(n_xcp), &
         neo(n_xcp), nio(n_xcp), Teo(n_xcp), Tio(n_xcp), &
         F1o(n_xcp), F2o(n_xcp), F3o(n_xcp), FPo(n_xcp), &
         F0(n_xcp), F0o(n_xcp), F0x(n_xcp), NN(n_xcp), TN(n_xcp), &
         cu(n_xcp), cutor(n_xcp), cd(n_xcp), cubs(n_xcp), &
         UPL(n_xcp), UPOL(n_xcp), EZ(n_xcp), ZEF(n_xcp), AMAIN(n_xcp) )
      !        include 'alloc.eq'
      ne = 0.
      ni = 0.
      Te = 0.
      Ti = 0.
      nex = 0.
      nix = 0.
      TEX = 0.
      TIX = 0.
      F1x = 0.
      F2x = 0.
      F3x = 0.
      F1 = 0.
      F2 = 0.
      F3 = 0.
      FP = 0.
      neo = 0.
      nio = 0.
      Teo = 0.
      Tio = 0.
      F1o = 0.
      F2o = 0.
      F3o = 0.
      FPo = 0.
      F0 = 0.
      F0o = 0.
      F0x = 0.
      NN = 0.
      TN = 0.
      cu = 0.
      cutor = 0.
      cd = 0.
      cubs = 0.
      UPL = 0.
      UPOL = 0.
      EZ = 0.
      ZEF = 0.
      AMAIN = 0.

      allocate ( &
         ametr(n_xcp), shif(n_xcp), vr(n_xcp), vro(n_xcp), &
         mu(n_xcp), VOL(n_xcp), VOLo(n_xcp), XCP(n_xcp), &
         IPOL(n_xcp), G11(n_xcp), G33(n_xcp), SLAT(n_xcp) )
      ametr = 0.
      shif = 0.
      vr = 0.
      vro = 0.
      mu = 0.
      VOL = 0.
      VOLo = 0.
      XCP = 0.
      IPOL = 0.
      G11 = 0.
      G33 = 0.
      SLAT = 0.

      allocate ( &
         ametre(n_xeq), shife(n_xeq), VOLe(n_xeq), &
         FPe(n_xeq), XEQ(n_xeq), RHO(n_xeq), &
         IPOLe(n_xeq), G11e(n_xeq), G33e(n_xeq), SLATe(n_xeq))
      !        include 'alloc.corsrs'
      ametre = 0.
      shife = 0.
      VOLe = 0.
      FPe = 0.
      XEQ = 0.
      RHO = 0.
      IPOLe = 0.
      G11e = 0.
      G33e = 0.
      SLATe = 0.

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
         DF0(n_xcp), VF0(n_xcp), DF1(n_xcp), VF1(n_xcp), &
         DF2(n_xcp), VF2(n_xcp), DF3(n_xcp), VF3(n_xcp), &
         PE(n_xcp), PET(n_xcp), PETOT(n_xcp), QE(n_xcp), QI(n_xcp), &
         PI(n_xcp), PIT(n_xcp), PITOT(n_xcp), PEI(n_xcp), &
         YPELSRS(n_xcp) )

      SN = 0.
      SNN = 0.
      SNTOT = 0.
      QN = 0.
      GN = 0.
      GNX = 0.
      SF0 = 0.
      SFF0 = 0.
      SF0TOT = 0.
      QF0 = 0.
      GF0 = 0.
      GF0X = 0.
      SF1 = 0.
      SFF1 = 0.
      SF1TOT = 0.
      QF1 = 0.
      GF1 = 0.
      GF1X = 0.
      SF2 = 0.
      SFF2 = 0.
      SF2TOT = 0.
      QF2 = 0.
      GF2 = 0.
      GF2X = 0.
      SF3 = 0.
      SFF3 = 0.
      SF3TOT = 0.
      QF3 = 0.
      GF3 = 0.
      GF3X = 0.
      DF0 = 0.
      VF0 = 0.
      DF1 = 0. 
      VF1 = 0.
      DF2 = 0.
      VF2 = 0.
      DF3 = 0.
      VF3 = 0.
      PE = 0.
      PET = 0.
      PETOT = 0.
      QE = 0.
      QI = 0.
      PI = 0.
      PIT = 0.
      PITOT = 0.
      PEI = 0.
      YPELSRS = 0.

      allocate ( PECR(n_xcp), CUECR(n_xcp) )
      !        include 'alloc.cortran'
      PECR = 0.
      CUECR = 0.

      allocate ( DN(n_xcp), CN(n_xcp), HE(n_xcp), XI(n_xcp), CC(n_xcp) )
      DN = 0.
      CN = 0.
      HE = 0.
      XI = 0.
      CC = 0.

      !====================================================== core profiles
      ne(1:n_xcp) = cp_in%profiles_1d(i_time)%electrons%density(1:n_xcp)/denA2D
      Te(1:n_xcp) = cp_in%profiles_1d(i_time)%electrons%temperature(1:n_xcp)/temA2D
      Ti(1:n_xcp) = cp_in%profiles_1d(i_time)%ion(1)%temperature(1:n_xcp)/temA2D
      MU(1:n_xcp) = 1./cp_in%profiles_1d(i_time)%q(1:n_xcp)

      !============================== detect hydrogen isotopes
      F0(1:n_xcp) = 0.
      F1(1:n_xcp) = 0.
      F2(1:n_xcp) = 0.
      F3(1:n_xcp) = 0.
      ispec(1:n_ion + 3) = 0
      Nhydr = 0
      istep = 0
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
      write (*, *) 'ispec', ispec
      if (ispec(1) .ne. 0) then
         iH = ispec(1)
         write (*, *) 'iH,iD,iT Nhydr', iH, iD, iT, Nhydr, n_xcp, i_time, denA2D
         write (*, *) cp_in%profiles_1d(i_time)%ion(iH)%density(1:n_xcp)/denA2D
         write (*, *) cp_in%profiles_1d(i_time)%neutral(iH)%density(1:n_xcp)/denA2D
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
         ni(j) = 0.
         Z2NdA = 0.
         do i = 1, n_ion
            ai = cp_in%profiles_1d(i_time)%ion(i)%element(1)%a
            zi = cp_in%profiles_1d(i_time)%ion(i)%z_ion_1D(j)
            ni(j) = ni(j) + cp_in%profiles_1d(i_time)%ion(i)%density(j)/denA2D
            if (ai .le. 0.) then
               write (*, *) 'wrong ion mass, i=', i
            else
               Z2NdA = Z2NdA + cp_in%profiles_1d(i_time)%ion(i)%density(j)/denA2D*zi**2/ai

            end if
            ZEF(j) = cp_in%profiles_1d(i_time)%ion(i)%density(j)/denA2D*&
                   & cp_in%profiles_1d(i_time)%ion(i)%z_ion_1D(j)**2/ne(j)
         end do                                        ! j ion
         if (ne(j) .ge. 0. .and. te(j) .gt. 0.) then
            pei(j) = 0.00246*(15.9 - .5*log(NE(j)) + log(TE(j)))*NE(j)*Z2NdA/TE(j)/sqrt(TE(j))
         else
            !                write(*,*) 'ne,te,j',ne(j),te(j),j
         end if
      end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! equilibrium
      RTOR = (eq_in%time_slice(j_time)%profiles_1d%r_outboard(n_xeq) + &
           &  eq_in%time_slice(j_time)%profiles_1d%r_inboard(n_xeq))/2.
      if (RTOR .gt. 0.) then
         BTOR = dabs(eq_in%vacuum_toroidal_field%b0(i_time))*eq_in%vacuum_toroidal_field%r0/RTOR
      else
         write (*, *) 'wrong RTOR =', RTOR
         stop
      end if
      XEQ(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(1:n_xeq)
      VOLe(1:n_xeq) = max(0., eq_in%time_slice(j_time)%profiles_1d%VOLUME(1:n_xeq))
      XCP(1:n_xcp) = cp_in%profiles_1d(i_time)%grid%rho_tor_norm(1:n_xcp)
      ROC = eq_in%time_slice(j_time)%profiles_1d%phi(n_xeq)
      ROC = sqrt(dabs(ROC/BTOR/3.141592))
      RHO(1:n_xcp) = ROC*XCP(1:n_xcp)
      ametre(1:n_xeq) = max(0., (eq_in%time_slice(j_time)%profiles_1d%r_outboard(1:n_xeq) - &
                              & eq_in%time_slice(j_time)%profiles_1d%r_inboard(1:n_xeq))/2.)
      shife(1:n_xeq) = (eq_in%time_slice(j_time)%profiles_1d%r_outboard(1:n_xeq) + &
                      & eq_in%time_slice(j_time)%profiles_1d%r_inboard(1:n_xeq))/2.-RTOR
      FPe(1:n_xeq) = -eq_in%time_slice(j_time)%profiles_1d%psi(1:n_xeq) + eq_in%time_slice(j_time)%profiles_1d%psi(1)
      IPOLe(1:n_xeq) = dabs(eq_in%time_slice(j_time)%profiles_1d%f(1:n_xeq)/RTOR/BTOR)
      SLATe(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%surface(1:n_xeq)
      G11e(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%gm3(1:n_xeq)
      G33e(1:n_xeq) = eq_in%time_slice(j_time)%profiles_1d%gm1(1:n_xeq)*RTOR**2
      write (*, *) 'Ip =', dabs(eq_in%time_slice(j_time)%global_quantities%Ip/cuA2D)
      ABC = ametre(n_xeq)
      SHIFT = Shife(n_xeq) - RTOR
      HRO = (XCP(3) - XCP(2))*ROC
      ALFA = 0.001d0
      !SMTH(ALFA,NO,FO,XO,N,FN,XN) mapping from EQ to CP grids
      !        write(*,*) 'trace 184'
      if (n_xeq .ne. n_xcp) then
         call SMTH(ALFA, n_xeq, VOLe, XEQ, n_xcp, VOL, XCP)
         call SMTH(ALFA, n_xeq, ametre, XEQ, n_xcp, ametr, XCP)
         call SMTH(ALFA, n_xeq, FPe, XEQ, n_xcp, FP, XCP)
         call SMTH(ALFA, n_xeq, shife, XEQ, n_xcp, shif, XCP)
         call SMTH(ALFA, n_xeq, IPOLe, XEQ, n_xcp, IPOL, XCP)
         call SMTH(ALFA, n_xeq, G11e, XEQ, n_xcp, G11, XCP)
         call SMTH(ALFA, n_xeq, G33e, XEQ, n_xcp, G33, XCP)
         call SMTH(ALFA, n_xeq, SLATe, XEQ, n_xcp, SLAT, XCP)
      else
         VOL(1:n_xcp) = VOLe(1:n_xeq)
         ametr(1:n_xcp) = ametre(1:n_xeq)
         FP(1:n_xcp) = FPe(1:n_xeq)
         shif(1:n_xcp) = shife(1:n_xeq)
         IPOL(1:n_xcp) = IPOLe(1:n_xeq)
         G11(1:n_xcp) = G11e(1:n_xeq)
         G33(1:n_xcp) = G33e(1:n_xeq)
         SLAT(1:n_xcp) = SLATe(1:n_xeq)
      end if
      VR(1) = VOL(2)/HRO
      G11(1) = VR(1)*G11(1)
      do j = 2, n_xcp
         VR(j) = (VOL(j) - VOL(j - 1))/HRO
         G11(j) = VR(j)*G11(j)
      end do


      TAU = smart_in%TAU
      dtau = smart_in%dtau
      !=========================================================== time loop
      open (1, file='out_Peltran.dat')
      TIME = TIMBEG
      !TIMPEL = 0.
      write (1, 997)   'TIME,s', '  Te(1)','   Ti(1)','   ne(1)',&
                     '   ni(1)', '   <ne>','   n0(1)','   n0(a)',&
                    '   <Shdt>','   <Sn0>',   '  QF0B'
997   format(11A14)
      !do 999 jtime=1,100
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
      VRo(1:NA1) = VR(1:NA1)

      !================================================== heat sources

      !=================================================== EC heating
      if (smart_in%sw_ech2a .ne. 0) then
         RHOEC  = smart_in%ROCEC*ROC ! EC location
         RHODR  = smart_in%ROCDR*ROC ! EC width
         QECR   = smart_in%QECR      ! QEC= 10 MW
         YEFFec = smart_in%YEFFec    ! IEC/QEC MA/MW

      !  ECH2a(YR0,YDR,YQ,YEFF,YP,YC,NA1,RHO,VR)
         call ECH2a(RHOEC, RHODR, QECR, YEFFec, PECR, CUECR, NA1, RHO, VR)
         PE(1:NA1) = PECR(1:NA1)
      endif
      !=================================================== alpha heating (simplified)
      do j = 1, na1
         ! SVDT [10#19m#3/s]:        The formula is a fit to D-T reaction rate
         !        according to Putvinskiy
         !        D+T=alpha(3.52MeV)+n(14.07MeV)
         !        Use: Palpha=Nd*Nt*SVDT*3520./625. [MW/m#3]
         !        (Yushmanov 11-JUN-87)
         SVDT = TI(J)**(-0.33333333)
         SVDT = 8.972*EXP(-19.9826*SVDT)*SVDT*SVDT*&
         &((TI(J) + 1.0134)/(1.+6.386E-3*(TI(J) + 1.0134)**2) +&
         &1.877*EXP(-.16176*TI(J)*SQRT(TI(J))))
         PDT = 5.632*f2(J)*f3(J)*SVDT
         ! PAION2 [MW/m#3]   D-T Fraction of fusion alpha power deposited to ions
         !     P.Pavlo  22.06.89/A.Polevoi 20-MAY-94
         y2 = 88./TE(J)
         y = sqrt(y2)
         PAION2 = 2.*(0.166666667*LOG((1.-y + y2)/(1.+2.*y + y2)) +&
         &0.57735026*(ATAN(0.57735026*(2.*y - 1.)) + 0.52359874))/y2

         pi(j) = PDT*paion2
         pe(j) = pe(j) + PDT*(1.-paion2)
      end do

      !        Write(*,*) 'before pellet'

      !=============================================== boundary conditions
      QNB = smart_in%QNB
      F01B = smart_in%F01B
      F02B = smart_in%F02B
      F03B = smart_in%F03B

      NEB = NE(NA1)
      TEB = TE(NA1)
      TIB = TI(NA1)
      F1B = F1(NA1)
      F2B = F2(NA1)
      F3B = F3(NA1)
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
         CN(J) = 0.
         DF1(j) = (HE(j) + XI(j))/10.
         VF1(J) = 0.
         DF2(j) = (HE(j) + XI(j))/10.
         VF2(J) = 0.
         DF3(j) = (HE(j) + XI(j))/10.
         VF3(J) = 0.
      end do

      !======================================== for neutral transport
      do j = 1, NA1
         svcx = 0.
         if (ti(j) .gt. 0 .and. amain(j) .ge. 1.) then
            SVCX = 10.**(5.9 + 0.3*LOG10(TI(J)/AMAIN(J)))
         end if
         DF0(J) = 9.584d10*(TI(J) + 1.d-9)/(SVCX + 1.d-10)/AMAIN(J)/NE&
         &(J)
         if (j .lt. NA) then
            VF0(J) = -7.6d-1*DF0(J)*(TI(J + 1) - TI(J))/HRO/(TI(J) + 1.d-9)
         else
            VF0(J) = -7.6d-1*DF0(J)*(TI(NA1) - TI(NA))/HRO/(TI(J) + 1.d-9)
         end if
         IF (TE(J) .LT. .0001) THEN
            SVRC = 0.
         ELSE
            SVRC = 13.6E-3/TE(j)
            SVRC = 1.27*SVRC*sqrt(SVRC)/(SVRC + .59)
         END IF
         SF0(J) = (F1(J) + F2(J) + F3(J))*NE(J)*SVRC
         SVIE = .0136/TE(J)
         IF (TE(J) .GT. .01) THEN
            SVIE = 9.7E5*EXP(-SVIE)*SQRT(SVIE/(1.+SVIE))/(SVIE + .73)
         ELSE
            SVIE = 2.958E5*EXP(-SVIE)*SQRT(SVIE)
         END IF
         SFF0(J) = -NEo(j)*SVIE
         !              SFF1(J)=-NEo(J)*SVRC
         !              SFF2(J)=-NEo(J)*SVRC
         !              SFF3(J)=-NEo(J)*SVRC

         SF1(j) = NEo(J)*(SVIE*F0(J)*F01B/F0B - SVRC*F1o(J))
         SF2(j) = NEo(J)*(SVIE*F0(J)*F02B/F0B - SVRC*F2o(J))
         SF3(j) = NEo(J)*(SVIE*F0(J)*F03B/F0B - SVRC*F3o(J))

         !        SF1(j)=NEo(J)*SVIE*F0(J)*F01B/F0B
         !        SF2(j)=NEo(J)*SVIE*F0(J)*F02B/F0B
         !        SF3(j)=NEo(J)*SVIE*F0(J)*F03B/F0B
         TN(j) = Ti(j) +&
         &(Te(j) - Ti(j))*195.*PEI(j)/SVCX/(ni(j) + f0(j))/ni(j)

      end do
      !        write(*,*) 'a1=',TEo(1:NA1)
      !==============================================density stepup
      !======================================= neutrals
      call STEPUPN0(&
         NA1, NB1, TAU, HRO, VRo, VR, G11, SLAT, RHO,&
         VF0, DF0, F0o, F0, F0X, QNB, SF0, SFF0, SF0TOT, QF0, GF0, GF0X)
      !=================================================
      do j = 1, NA1
         svcx = 0.
         if (ti(j) .gt. 0 .and. amain(j) .ge. 1.) then
            SVCX = 10.**(5.9 + 0.3*LOG10(TI(J)/AMAIN(J)))
         end if
         DF0(J) = 9.584d10*(TI(J) + 1.d-9)/(SVCX + 1.d-10)/AMAIN(J)/NE&
         &(J)
         if (j .lt. NA) then
            VF0(J) = -7.6d-1*DF0(J)*(TI(J + 1) - TI(J))/HRO/(TI(J) + 1.d-9)
         else
            VF0(J) = -7.6d-1*DF0(J)*(TI(NA1) - TI(NA))/HRO/(TI(J) + 1.d-9)
         end if
         IF (TE(J) .LT. .0001) THEN
            SVRC = 0.
         ELSE
            SVRC = 13.6E-3/TE(j)
            SVRC = 1.27*SVRC*sqrt(SVRC)/(SVRC + .59)
         END IF
         SF0(J) = (F1(J) + F2(J) + F3(J))*NE(J)*SVRC
         SVIE = .0136/TE(J)
         IF (TE(J) .GT. .01) THEN
            SVIE = 9.7E5*EXP(-SVIE)*SQRT(SVIE/(1.+SVIE))/(SVIE + .73)
         ELSE
            SVIE = 2.958E5*EXP(-SVIE)*SQRT(SVIE)
         END IF
         SFF0(J) = -NEo(j)*SVIE
         !              SFF1(J)=-NEo(J)*SVRC
         !              SFF2(J)=-NEo(J)*SVRC
         !              SFF3(J)=-NEo(J)*SVRC

         SF1(j) = NEo(J)*(SVIE*F0(J)*F01B/F0B - SVRC*F1o(J))
         SF2(j) = NEo(J)*(SVIE*F0(J)*F02B/F0B - SVRC*F2o(J))
         SF3(j) = NEo(J)*(SVIE*F0(J)*F03B/F0B - SVRC*F3o(J))

         !        SF1(j)=NEo(J)*SVIE*F0(J)*F01B/F0B
         !        SF2(j)=NEo(J)*SVIE*F0(J)*F02B/F0B
         !        SF3(j)=NEo(J)*SVIE*F0(J)*F03B/F0B
         TN(j) = Ti(j) + (Te(j) - Ti(j))*195.*PEI(j)/SVCX/(ni(j) + f0(j))/ni(j)
      end do

      !======================================= hydrogen species
      if (iH .ne. 0) call STEPUPN( &
         NA1, NB1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
         VF1, DF1, F1o, F1, F1X, F1B, SF1, SFF1, SF1TOT, QF1, GF1, GF1X)
      if (iD .ne. 0) call STEPUPN( &
         NA1, NB1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
         VF2, DF2, F2o, F2, F2X, F2B, SF2, SFF2, SF2TOT, QF2, GF2, GF2X)
      if (iT .ne. 0) call STEPUPN( &
         NA1, NB1, TAU, HRO, VRo, VR, G11, SLAT, RHO, &
         VF3, DF3, F3o, F3, F3X, F3B, SF3, SFF3, SF3TOT, QF3, GF3, GF3X)
      !============================ electron density from quasineutrality
      ne(1:na1) = f1(1:na1) + f2(1:na1) + f3(1:na1)
      ni(1:na1) = ne(1:na1)
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
      call STEPUPT( &
         NA1, NB1, TAU, HRO, VRo, VR, G11, SLAT, RHO, XI, HE, &
         PE, PET, PETOT, PI, PIT, PITOT, PEI, &
         TEX, TE, TEo, TEB, TIX, TI, TIo, TIB, &
         NEX, NEo, NE, NIX, NIo, NI, Qe, Qi, GNX, GN2E, GN2I)

      write (*, *) 'QE,QI,Ge', QE(NA1), QI(NA1), QF1(NA1) + QF2(NA1) + QF3(NA1)
 
      !============================================= pelshot
      !        Write(*,*) 'before pellet'
      !        time=time + TAU
      TIMPEL = TIMPEL + TAU
      if (TIMPEL .ge. (dtau - 1.d-7)) then
         TIMPEL = 0.
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
      jabs = (1.-YDABL)*NA1
      if (smart_in%YCOS0 .gt. 0.) then
         Dtpel = (smart_in%yglength + smart_in%YDL*(ametr(NA1) - ametr(jabs) + shif(jabs)))/smart_in%YVP/1000.
      else
         Dtpel = (smart_in%yglength + smart_in%YDL*(ametr(NA1) - ametr(jabs) - shif(jabs)))/smart_in%YVP/1000.
      end if
      !        write(*,*) 'Dtpel =',Dtpel
      !============================================ current diffusion (+equilibrium)
      !============================================ end of itterations
      write (1, 998) &
         TIME, Te(1), Ti(1), ne(1), ni(1), &
         VINTa(ne, ROC, RHO, VR, NA1)/VOLe(n_xeq), F0(1), F0(NA1), &
         VINTa(SF3TOT, ROC, RHO, VR, NA1) + &
         VINTa(SF2TOT, ROC, RHO, VR, NA1) + &
         VINTa(SF1TOT, ROC, RHO, VR, NA1),  &
         VINTa(SF0TOT, ROC, RHO, VR, NA1), QNB

      TIME = TIME + TAU

      !999 continue
      close (1)
998   format(11(1X0PE13.6))
      ! end of time loop
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
      write (*, *) 'YDABL,YDDEP', YDABL, YDDEP
      write (*, *) 'RHOEC,RHODR,ROC,QEC,YEFFec=',RHOEC, RHODR, ROC, QECR, YEFFec
      Write (*, *) 'Pecr', VINTa(PECR, ROC, RHO, VR, NA1),&
                   'Pe', VINTa(PECR, ROC, RHO, VR, NA1),&
                   'Pi', VINTa(PI, ROC, RHO, VR, NA1)
      deallocate (ispec)
      !        include 'dealloc.corprf
      deallocate (ne, ni, Te, Ti, nex, nix, TEX, TIX, TN, NN,&
                  F0, F1, F2, F3, F0x, F1x, F2x, F3x, FP,&
                  neo, nio, Teo, Tio, F0o, F1o, F2o, F3o, FPo,&
                  cu, cutor, cd, cubs, UPL, UPOL, EZ, ZEF)
      !        include 'dealloc.corsrs'
      deallocate (SN, SNN, SNTOT, QN, GN, GNX, PE, PET, PETOT,&
                  SF0, SFF0, SF0TOT, QF0, GF0, GF0X,&
                  SF1, SFF1, SF1TOT, QF1, GF1, GF1X,&
                  SF2, SFF2, SF2TOT, QF2, GF2, GF2X,&
                  SF3, SFF3, SF3TOT, QF3, GF3, GF3X,&
                  PI, PIT, PITOT, PEI, QE, QI, PECR, CUECR,YPELSRS)
      !        include 'dealloc.cortran'
      deallocate (DF0, VF0, DF1, VF1, DF2, VF2, DF3, VF3,DN, CN, HE, XI, CC)
      !        write(*,*) '308'
      !        include 'dealloc.eq'
      deallocate (IPOL, G11, G33, SLAT, ametr, shif, vr, vro,&
                  mu, VOL, XCP, ametre, shife, VOLe,&
                  FPe, XEQ, RHO, IPOLe, G11e, G33e, SLATe)
      !        write(*,*) '313'

   end subroutine smart

end module mod_smart

!======================================================================|
SUBROUTINE ECH2a(YR0, YDR, YQ, YEFF, YP, YC, NA1, RHO, VR)
   !                 RHOEC,RHODR,QEC,YEFFec,PEC,CUECR,NA1,RHO,VR
   !-----------------------------------------25-JUN-93
   ! Artificial parabolic Heating + CD                POLEVOY
   ! YR0   distance from plasma centre [m]
   ! YDR   half width [m]
   ! Q     ECH Power [MW]
   ! YP(r) ECR power density [MW/m3] ~ Q*exp(-((r-r0)/Dr)**2)
   ! YC(r) ECR current density [MW/m3]
   ! YEFF ECR current drive efficiency [A/W]
   implicit none
   integer j, NA1
   double precision YP(*), YC(*), RHO(*), VR(*)
   !        double precision IINT, VINT
   double precision VINTa, YS, YM, YR0, YDR, YQ, YEFF, ROC
   YS = 0.d0
   YM = 0.d0
   ROC = RHO(NA1)
   if (ydr .le. 1.d-7) then
      write (*, *) 'sbr/ech2a: too small iput power width'
      !                 YDR=ABC
      YDR = RHO(NA1)
   end if
   do J = 1, NA1
      YP(J) = -((RHO(J) - YR0)/YDR)**2
      if (YP(J) .gt. -18.d0) then
         YP(J) = dexp(YP(J))
      else
         YP(J) = 0.d0
      end if
   end do
   !                YM=IINT(YP,ROC)
   YS = VINTa(YP, ROC, RHO, VR, NA1)
   !        write(*,*) 'YS=',YS
   !        if(ym.le.1.d-7.or.ys.le.1.d-7) then

   if (ys .le. 0.) then
      YP(1:na1) = 0.
   else
      YP(1:na1) = YQ/(YS + 1.d-7)*YP(1:na1)
   end if
   if (ym .le. 0.) then
      YC(1:na1) = 0.
   else
      YC(1:na1) = YP(1:na1)*YEFF/(YM + 1.d-7)
   end if
   return
end subroutine ECH2a 

!=======================================================================
! VINT: Volume integral {0,R} of any array
! Only a radially dependent array may be the 1st parameter of the function
!                                                  (Yushmanov 26-DEC-90)
! Examples:
!    out\Vint(CAR3)            !Radial profile of CAR3 volume integral
!    out_Vint(CAR3,Ro); !Volume integral {0,Ro} of CAR3
!    out_Vint(CAR3B)    !Total volume integral of CAR3 (0,ROC)
double precision function VINTa(ARR, YR, RHO, VR, NA1)
   implicit none
   double precision ARR(*), VR(*), RHO(*), YR, YDR, YR1, HRO, HROA
   integer JK, J, NA1, NA
   VINTa = 0.
   if (YR .le. 0.) return
   NA = NA1 - 1
   HRO = RHO(3) - RHO(2)
   HROA = RHO(NA1) - RHO(NA)
   if (YR .le. RHO(NA)) then
      JK = YR/HRO + 1
      YR1 = YR
   else
      JK = NA
      YR1 = min(YR, RHO(NA) + .5d0*HROA)
   end if
   YDR = (JK - YR1/HRO)*VR(JK)
   do J = 1, JK
      VINTa = VINTa + ARR(J)*VR(J)
   end do
   VINTa = HRO*(VINTa - ARR(JK)*YDR)
end function VINTa
!======================================================================|
