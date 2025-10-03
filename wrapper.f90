program wrapper

  use nemo_module
  use mod_kind
  use ids_schemas
  use ids_routines
  use f90_file_reader, only: file2buffer
  use codeparam_input_wrapper ! ROUTINES FOR READING INPUT XML FILE FOR SPOT WRAPPER
#ifdef AMNS
  use itm_types
  use amns_module
  use amns_types
#endif

  implicit none

#ifdef AMNS
  type (amns_handle_type) :: amns                            ! AMNS global handle
  type (amns_handle_rx_type)  :: amns_rx                     ! AMNS table handle
  type (amns_reaction_type) :: xx_rx
  type (amns_reactants_type) :: species
  type (amns_version_type) :: version
  type (amns_error_type) :: error_status
  real (kind=R8) :: te=100.0_R8, ne=1e20_R8, rate, ne_tab(5)
  integer :: i
#endif

  character(len=:), pointer:: output_message
  integer:: output_flag

  character(len=200):: testmode

  integer(kind=inm):: shot_number,input_run,output_run,itime,ntime
  integer(kind=inm):: ierror,idx_in,idx_out,error_flag,myrank
  integer(kind=inm), parameter:: master_rank=0
  integer(kind=inm), parameter:: iounit=1

  type(ids_parameters_input):: codeparam_nemo,codeparam_wrapper
  type(type_wrapper_codeparam_data):: codeparam_wrapper_data

  type(ids_equilibrium):: equilibrium_input
  type(ids_core_profiles):: core_profiles_input
  type(ids_nbi):: nbi_input
  type(ids_wall):: wall_input
  type(ids_distribution_sources):: distribution_sources_output

  type(ids_equilibrium):: equilibrium_input_1t
  type(ids_core_profiles):: core_profiles_input_1t
  type(ids_nbi):: nbi_input_1t
  type(ids_wall):: wall_input_1t
  type(ids_distribution_sources):: distribution_sources_output_1t

  double precision:: timenow,time_to_simulate
  character(len=200):: user_in,user_local,machine
  integer :: return_status

  ! ---------------------------------------------------------------

#ifdef MPI
  call mpi_init(ierror)
#endif
  myrank = 0

  ! ARE WE IN REGRESSION/UNITARY TEST MODE OR NORMAL MODE?
#ifdef FRUIT
  call getenv('MS_TEST_MODE',testmode)
#else
  testmode = '0' ! Thomas: Note that the variable testmode is used also when FRUIT is turned off.
#endif

  ! ---------------------------------------
  ! USING THE WRAPPER FOR REGRESSION TESTS
  ! ---------------------------------------
  if(trim(testmode).eq.'1'.or.trim(testmode).eq.'2') then
     ! READING FROM XML DATA FOR NEMO
     call file2buffer('input/input_nemo_reg.xml',iounit,codeparam_nemo%parameters_value)

     ! READING FROM XML DATA FOR WRAPPER
     call file2buffer('input/input_wrapper_reg.xml',iounit,codeparam_wrapper%parameters_value)

     ! -----------------------------------
     ! USING THE WRAPPER FOR ALL THE REST
     ! -----------------------------------
  else
     ! READING FROM XML DATA FOR NEMO
     call file2buffer('input/input_nemo.xml',iounit,codeparam_nemo%parameters_value)

     ! READING FROM XML DATA FOR WRAPPER
     call file2buffer('input/input_wrapper.xml',iounit,codeparam_wrapper%parameters_value)

  endif

  ! READ INPUT XML FILE
  call parse_wrapper_codeparam(codeparam_wrapper%parameters_value,codeparam_wrapper_data)

  ! INPUT AND OUTPUT PARAMETERS
  shot_number      = codeparam_wrapper_data%shot_number
  input_run        = codeparam_wrapper_data%input_run
  output_run       = codeparam_wrapper_data%output_run
  user_in          = codeparam_wrapper_data%user_in
  machine          = codeparam_wrapper_data%machine
  time_to_simulate = codeparam_wrapper_data%time_to_simulate

  ! USERNAME DEFINED BY ENVIRONMENT VARIABLE USERNAME
  call getenv('USER',user_local)
  if(user_in.eq.'user_environment_variable') then
     user_in = user_local
  endif

  ! OPEN INPUT DATAFILE
  if(myrank.eq.master_rank) write(*,*) 'Open input datafile'
  call imas_open_env('ids',shot_number,input_run,idx_in,user_in,machine,'3', return_status)
  if(return_status.ne.0) then
     if(myrank.eq.master_rank) then
        write(*,*) '-------------------------------------------------------'
        write(*,*) ' Cannot read input datafile:'
        write(*,*) '    shot = ',shot_number,' run = ',input_run
        write(*,*) '    user = "', trim(adjustl(user_in)), '", machine = "', trim(adjustl(machine)),'"'
        write(*,*) '=> Program stopped.'
        write(*,*) '-------------------------------------------------------'
        call abort('')
     endif
  endif

  ! CREATE OUTPUT DATAFILE
  if(myrank.eq.master_rank) then
     write(*,*) 'Create output datafile'
     write(*,*) ' - shot=',shot_number,', run=',output_run
     write(*,*) ' - user="',trim(adjustl(user_local)),'", machine="',trim(adjustl(machine)),'"'
  end if
  call imas_create_env('ids',shot_number,output_run,0,0,idx_out,user_local,machine,'3',return_status)  
  if (return_status /= 0) call abort('Error received from imas_create_env')

  ! READ CORE_PROFILES, EQUILBRIUM, NBI, WALL
  write(*,*) 'Read input IDSs'
  write(*,*) '=> Read equilibrium'
  call ids_get(idx_in,'equilibrium',equilibrium_input,return_status)
  if (return_status /= 0) call abort('Error received from ids_get while reading equilibrium')

  write(*,*) '=> Read core_profiles'
  call ids_get(idx_in,'core_profiles',core_profiles_input,return_status)
  if (return_status /= 0) call abort('Error received from ids_get while reading core_profiles')

  write(*,*) '=> Read nbi'
  call ids_get(idx_in,'nbi',nbi_input,return_status)
  if (return_status /= 0) call abort('Error received from ids_get while reading nbi')

  write(*,*) '=> Read wall'
  call ids_get(idx_in,'wall',wall_input,return_status)
  if (return_status /= 0) call abort('Error received from ids_get while reading wall')

  write(*,*) 'Finished reading input IDSs.'

  ! ---------------------------------------------------------

#ifdef AMNS
  ne_tab = (/1e19_R8, 3e19_R8, 6e19_R8, 1e20_R8, 2e20_R8/)

  ! Define major version and localisation of IMAS data entry where the atomic data tables are stored
  version%string = '3'
  version%number = 3
  shot = 4
  run  = 1
  call imas_amns_setup(amns, error_status, version, shot, run, user_local ) ! set up the AMNS system
  allocate(species%components(1)) ! set up reactants (for which we want to access the tables for a particular reaction (defined below)
  species%components(1)%zn=4
  species%components(1)%za=1
  species%components(1)%mi=9
  xx_rx%string='EI'                                           ! set up reaction
  call imas_amns_setup_table(amns, xx_rx, species, amns_rx)   ! set up table
  do i = 1,5
     call imas_amns_rx(amns_rx, rate, te, ne_tab(i))          ! get results
     write(*,*) 'Rate = ', rate
  enddo
  call imas_amns_finish_table(amns_rx)                        ! finish with table
  call imas_amns_finish(amns)                                 ! finish with amns
#endif

  ! ---------------------------------------------------------

  ! CHOOSE TIME SLICE FROM INPUT, I.E. ONE TIME SLICE (IF TIME_TO_SIMULATE<0) OR FROM IDS (THEN ACTIVATE TIME LOOP)
  ! 1) IF TIME_TO_SIMULATE < 0: USE IDS FIRST TIME SLICE
  ! 2) IF TIME_TO_SIMULATE = 0: USE IDS WHOLE TIME VECTOR
  ! 3) IF TIME_TO_SIMULATE > 0: USE THE SPECIFIED TIME SLICE
  if(time_to_simulate.eq.0) then
     ntime = size(core_profiles_input%time)
  else
     ntime = 1
  endif

  call ids_copy(equilibrium_input,equilibrium_input_1t)
  call ids_copy(core_profiles_input,core_profiles_input_1t)
  call ids_copy(nbi_input,nbi_input_1t)
  call ids_copy(wall_input,wall_input_1t)

  do itime=1,ntime

     if(time_to_simulate.eq.0.) then
        timenow = core_profiles_input%time(itime)
     elseif(time_to_simulate.lt.0.) then
        timenow = 0.
     elseif(time_to_simulate.gt.0.) then
        timenow = time_to_simulate
     endif

     if(myrank.eq.master_rank) write(*,*) 'get equilibrium slice'
     call ids_get_slice(idx_in,'equilibrium',equilibrium_input_1t,timenow,1,return_status)
     if (return_status /= 0) call abort('Error received from ids_get_slice while reading equilibrium')

     if(myrank.eq.master_rank) write(*,*) 'get core_profiles slice'
     call ids_get_slice(idx_in,'core_profiles',core_profiles_input_1t,timenow,1,return_status)
     if (return_status /= 0) call abort('Error received from ids_get_slice while reading equilibrium')

     if(myrank.eq.master_rank) write(*,*) 'get nbi slice'
     call ids_get_slice(idx_in,'nbi',nbi_input_1t,timenow,1,return_status)
     if (return_status /= 0) call abort('Error received from ids_get_slice while reading nbi')

     if(myrank.eq.master_rank) write(*,*) 'get wall slice'
     call ids_get_slice(idx_in,'wall',wall_input_1t,timenow,1,return_status)
     if (return_status /= 0) call abort('Error received from ids_get_slice while reading wall')

     ! test provisoire
     !nbi_input_1t%unit(1)%energy%data(1) = 530e3
     !nbi_input_1t%unit(2)%energy%data(1) = 530e3
     !nbi_input_1t%unit(1)%power_launched%data(1) = 4.7e6
     !nbi_input_1t%unit(2)%power_launched%data(1) = 4.7e6
     !nbi_input_1t%unit(1)%species%a=1
     !nbi_input_1t%unit(1)%species%z_n=1
     !nbi_input_1t%unit(2)%species%a=1
     !nbi_input_1t%unit(2)%species%z_n=1

     ! RE-ADJUST THE TIME TO THE ACTUAL ONE AS WRITTEN IN THE IDS
     timenow = core_profiles_input_1t%time(1)

     if(myrank.eq.master_rank.and.itime.eq.1) write(*,'(a9,i5)')      'ntime = ',ntime
     if(myrank.eq.master_rank.and.ntime.gt.1) write(*,'(a9,i5)')      'itime = ',itime
     if(myrank.eq.master_rank)                write(*,'(a9,f8.2,a3)') 'Time  = ',timenow,' s'

     if(nbi_input_1t%unit(1)%power_launched%data(1).gt.0.) then

        ! EXECUTION OF THE NEMO ACTOR
        write(*,*) 'Execute nemo actor'
        !print*,'rho_tor',equilibrium_input_1t%time_slice(1)%profiles_1d%rho_tor
        call nemo(equilibrium_input_1t,core_profiles_input_1t,nbi_input_1t,wall_input_1t &
             ,distribution_sources_output_1t,codeparam_nemo,output_flag,output_message)

        if(output_flag.eq.0) then

           ! EXPORT RESULTS TO LOCAL DATABASE
           write(*,*) '--------------------------------'
           write(*,*) 'EXPORT RESULTS TO LOCAL DATABASE'
           write(*,*) '--------------------------------'

           if(itime.eq.1) then
              call ids_put(idx_out,'equilibrium',equilibrium_input_1t)
              call ids_put(idx_out,'core_profiles',core_profiles_input_1t)
              call ids_put(idx_out,'nbi',nbi_input_1t)
              if(wall_input_1t%ids_properties%homogeneous_time.ge.0) &
                   call ids_put(idx_out,'wall',wall_input_1t)
              call ids_put(idx_out,'distribution_sources',distribution_sources_output_1t)
           endif

           write(*,*) 'Write equilibrium IDS'
           call ids_put_slice(idx_out,'equilibrium',equilibrium_input_1t)
           write(*,*) 'Write core_profiles IDS'
           call ids_put_slice(idx_out,'core_profiles',core_profiles_input_1t)
           write(*,*) 'Write nbi IDS'
           nbi_input_1t%ids_properties%homogeneous_time = 1
           call ids_put_slice(idx_out,'nbi',nbi_input_1t)
           write(*,*) 'Write wall IDS'
           if(wall_input_1t%ids_properties%homogeneous_time.ge.0) &
                call ids_put_slice(idx_out,'wall',wall_input_1t)
           distribution_sources_output_1t%ids_properties%homogeneous_time = 1
           write(*,*) 'Write distribution_sources IDS'
           call ids_put(idx_out,'distribution_sources',distribution_sources_output_1t) ! homogenous_time=0 for ggd(2) => ids_put_slice does not work

        else

           write(*,*) output_message
           write(*,*) '=> Program stopped.'
           exit

        endif

     else
        write(*,*) 'No NBI power for this time slice'
     endif

  enddo

  if(myrank.eq.master_rank) write(*,*) 'Close input and output datafiles'
  call imas_close(idx_in)
  call imas_close(idx_out)

#ifdef MPI
  call mpi_finalize(ierror)
#endif

contains

  subroutine abort(message)
    character(*) :: message
    write(0,*)message
#ifdef MPI
    call mpi_abort()
#else
    stop
#endif
  end subroutine abort

end program wrapper
