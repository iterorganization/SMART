module mod_smart

contains

  subroutine smart(eq_in, cp_in, pellets_in, cp_out, codeparam, error_flag, error_message)

    ! ---------------------------------------
    ! PELLET ABLATION MODEL FROM ASTRA: SMART
    ! IDS INPUT: EQUILIBRIUM, CORE_PROFILES
    ! IDS OUTPUT: CORE_PROFILES
    ! ---------------------------------------

    use ids_schemas, only: ids_equilibrium, ids_core_profiles, ids_pellets
    use ids_schemas, only: ids_parameters_input, ids_is_valid
    use ids_routines, only: ids_copy
    use mod_codeparam_smart

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
    double precision ::                                         &
        denA2D, temA2D, presA2D, cuA2D,                         &
        HRO, ROC, BTOR, GP, SHIFT, ABC, RTOR, ALFA,             &
        YAM, YVP, YVOL, YCOS0, YEFF, YDL,YDABL, YDDEP, yswitch
    double precision, allocatable ::                            &
        ne(:), ni(:), Te(:), Ti(:), F1(:), F2(:), F3(:), FP(:), &
        ametr(:), shif(:), vr(:), mu(:), YPELSRS(:), VOL(:),    &
        ametre(:), shife(:), VOLe(:), XCP(:), XEQ(:), RHO(:), FPe(:)

    ! INITIALISATION OF ERROR FLAG
    error_flag = 0
    nullify(error_message) ! do this otherwise gfortran behaviour is undefined

    ! INITIAL DISPLAY
    write(*,*) ' '
    write(*,*) '======================================='
    write(*,*) 'START OF PELLET ABLATION MODEL: SMART'

    ! CHECK IF INPUT PELLETS IDS IS FILLED: IF SO, USE IT INSTEAD OF CODEPARAM
    ! FOR PELLET DESCRIPTION
    if(pellets_in%ids_properties%homogeneous_time.ge.0) from_pellets_ids = .True.
    
    ! CHECK IF INPUT IDS IS VALID
    if (ids_is_valid(eq_in%ids_properties%homogeneous_time)      .and. &
        size(eq_in%time)>0                                       .and. &
        ids_is_valid(cp_in%ids_properties%homogeneous_time)      .and. &
        size(cp_in%time)>0) then

       call assign_codeparam(codeparam%parameters_value,smart_in)

       if(from_pellets_ids.eq..true.) then ! REPLACE PELLET INFORMATION
          write(*,*) 'Input pellets IDS detected'
          smart_in%YAM  = pellets_in%time_slice(1)%pellet(1)%species(1)%a
          smart_in%YVP  = pellets_in%time_slice(1)%pellet(1)%velocity_initial*1.e-3
          smart_in%YVOL = pellets_in%time_slice(1)%pellet(1)%shape%size(1)*1e9
       else
          write(*,*) 'Input pellets IDS NOT detected'
       endif
       
       write(*,*) '------------------------------------'
       write(*,*) 'Parameters read from input xml file:'
       write(*,'(a25,f7.3)') ' YAM                   = ', smart_in%YAM
       write(*,'(a25,f7.3)') ' YVP                   = ', smart_in%YVP
       write(*,'(a25,f7.3)') ' YVOL                  = ', smart_in%YVOL
       write(*,'(a25,f7.3)') ' YCOS0                 = ', smart_in%YCOS0
       write(*,'(a25,f7.3)') ' YEFF                  = ', smart_in%YEFF
       write(*,'(a25,f7.3)') ' YDL                   = ', smart_in%YDL
       write(*,'(a25,f7.3)') ' yswitch               = ', smart_in%yswitch
       write(*,*) '------------------------------------'

       ! COPY THE INPUT IDS IN THE OUTPUT IDS
       call ids_copy(cp_in, cp_out)
    else
       error_flag = -1
       allocate(character(50):: error_message)
       error_message = 'Error in SMART: input IDS not valid'
       return
    endif


    !== pelTRY.f

    i_time = size(cp_in%profiles_1d(:)%time)
    j_time = size(eq_in%time_slice)

    n_xcp  = size(cp_in%profiles_1d(i_time)%grid%rho_tor_norm)
    n_ion  = size(cp_in%profiles_1d(i_time)%ion)
    n_xeq  = size(eq_in%time_slice(j_time)%profiles_1d%psi)

    nrd    = 2*n_xcp
    NA1    = n_xcp

    denA2D  = 1.d19
    temA2D  = 1.d3
    presA2D = 1.6d3 ! in ASTRA pressure is keV 10^19m-3, have to change Pa
    cuA2D   = 1.d6


    allocate(ne(n_xcp), ni(n_xcp), Te(n_xcp), Ti(n_xcp), &
             F1(n_xcp), F2(n_xcp), F3(n_xcp), FP(n_xcp), &
             ametr(n_xcp), shif(n_xcp), vr(n_xcp),       &
             mu(n_xcp), YPELSRS(n_xcp), VOL(n_xcp), XCP(n_xcp))

    allocate(ametre(n_xeq), shife(n_xeq), VOLe(n_xeq), &
             FPe(n_xeq), XEQ(n_xeq), RHO(n_xeq))

    BTOR = dabs(eq_in%vacuum_toroidal_field%b0(i_time))
    RTOR = eq_in%vacuum_toroidal_field%r0

    ne(:) = cp_in%profiles_1d(i_time)%electrons%density(:)/denA2D
    Te(:) = cp_in%profiles_1d(i_time)%electrons%temperature(:)/temA2D
    F1(:) = cp_in%profiles_1d(i_time)%ion(1)%density(:)/denA2D
    F2(:) = cp_in%profiles_1d(i_time)%ion(2)%density(:)/denA2D
    F3(:) = cp_in%profiles_1d(i_time)%ion(3)%density(:)/denA2D
    Ti(:) = cp_in%profiles_1d(i_time)%ion(1)%temperature(:)/temA2D
    MU(:) = 1./cp_in%profiles_1d(i_time)%q(:)
    XEQ(:)= eq_in%time_slice(j_time)%profiles_1d%rho_tor_norm(:)
    VOLe(:)= eq_in%time_slice(j_time)%profiles_1d%VOLUME(:)
    XCP(:)= cp_in%profiles_1d(i_time)%grid%rho_tor_norm(:)
    ROC = eq_in%time_slice(j_time)%profiles_1d%phi(n_xeq)
    ROC = sqrt(dabs(ROC/BTOR/3.141592))
    ametre(:) = (eq_in%time_slice(j_time)%profiles_1d%r_outboard(:) - &
                       eq_in%time_slice(j_time)%profiles_1d%r_inboard(:))/2.
    shife(:) = (eq_in%time_slice(j_time)%profiles_1d%r_outboard(:) + &
                      eq_in%time_slice(j_time)%profiles_1d%r_inboard(:))/2. - RTOR
    FPe(:) = - eq_in%time_slice(j_time)%profiles_1d%psi(:) &
                   + eq_in%time_slice(j_time)%profiles_1d%psi(1)

    ABC = ametre(n_xeq)
    SHIFT = Shife(n_xeq) - RTOR
    HRO = (XCP(3) - XCP(2)) * ROC
    ALFA = 0.001d0

    !SMTH(ALFA,NO,FO,XO,N,FN,XN) mapping from EQ to CP grids
    call SMTH(ALFA, n_xeq, VOLe, XEQ, n_xcp, VOL, XCP, NRD)
    call SMTH(ALFA, n_xeq, ametre, XEQ, n_xcp, ametr, XCP, NRD)
    call SMTH(ALFA, n_xeq, FPe, XEQ, n_xcp, FP, XCP, NRD)
    call SMTH(ALFA, n_xeq, shife, XEQ, n_xcp, shif, XCP, NRD)

    VR(1) = VOL(1) / HRO
    do j=2, n_xcp
        VR(j) = (VOL(j) - VOL(j-1)) / HRO
    enddo

    do i=1, n_xcp
        ni(i) = 0.
        do j=1, n_ion
            ni(i) = ni(i) + cp_in%profiles_1d(i_time)%ion(j)%density(i)/denA2D
        enddo
    enddo

    !== Pellet Ablation Model: SMART 
    call pelIMAS(smart_in%YAM,   smart_in%YVP,  smart_in%YVOL, &
                 smart_in%YCOS0, smart_in%YEFF, smart_in%YDL,  &
                 YDABL, YDDEP, YPELSRS, smart_in%yswitch,      &
                 ne, ni, Te, Ti, F1, F2, F3, FP,               &
                 ametr, shif, vr, mu,                          &
                 HRO, ROC, BTOR, RTOR, NA1, NRD)

    !== conversion to IMAS units
    cp_out%profiles_1d(i_time)%electrons%density(:) = ne(:)*denA2D
    cp_out%profiles_1d(i_time)%electrons%temperature(:) = Te(:)*temA2D
    cp_out%profiles_1d(i_time)%ion(1)%density(:) = F1(:)*denA2D
    cp_out%profiles_1d(i_time)%ion(2)%density(:) = F2(:)*denA2D
    cp_out%profiles_1d(i_time)%ion(3)%density(:) = F3(:)*denA2D
    do j=1, n_ion
        cp_out%profiles_1d(i_time)%ion(j)%temperature(:) = Ti(:)*temA2D
    enddo

    deallocate(ne, ni, Te, Ti, F1, F2, F3, FP, &
               ametr, shif, vr, mu, YPELSRS, VOL, XCP)
    deallocate(ametre, shife, VOLe, FPe, XEQ, RHO)


    ! FINAL DISPLAY
    write(*,*) 'END OF SMART'
    write(*,*) '======================================='
    write(*,*) ' '

  end subroutine smart

end module mod_smart
