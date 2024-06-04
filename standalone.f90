program standalone

  use ids_schemas
  use ids_routines
  use f90_file_reader, only: file2buffer
  use mod_codeparam_standalone
  use mod_smart

  implicit none

  integer:: idx,run_in,run_out,pulse,error_flag
  integer:: iounit = 1
  integer:: status
  type(ids_equilibrium):: equilibrium_in
  type(ids_core_profiles):: core_profiles_in, core_profiles_out
  type(ids_pellets) :: pellets_in
  type(ids_parameters_input):: codeparam_standalone,codeparam_smart
  type(type_standalone_data):: standalone_in
  character(len=200):: input_db,input_machine,local_db,local_machine
  character(len=:), pointer:: error_message
  character (STRMAXLEN) :: uri

  ! READ STANDALONE XML INPUT FILE, DEFINED BY SPECIFIC XSD FILE
  call file2buffer('input/standalone.xml',iounit, codeparam_standalone%parameters_value)
  call assign_codeparam(codeparam_standalone%parameters_value,standalone_in)
  write(*,*) '------------------------------------'
  write(*,*) 'Parameters read from input xml file:'
  write(*,'(a17,i7)')  ' pulse         = ',standalone_in%pulse
  write(*,'(a17,i7)')  ' run_in        = ',standalone_in%run_in
  write(*,'(a17,i7)')  ' run_out       = ',standalone_in%run_out
  write(*,'(a17,a30)') ' input_db      = ',standalone_in%input_db
  write(*,'(a17,a30)') ' input_machine = ',standalone_in%input_machine
  write(*,'(a17,a30)') ' local_db      = ',standalone_in%local_db
  write(*,'(a17,a30)') ' local_machine = ',standalone_in%local_machine
  write(*,*) '------------------------------------'

  ! READ SMART XML INPUT FILE (ASSIGN_CODEPARAM IS EXECUTED INSIDE THE SMART ROUTINE)
  call file2buffer('input/smart.xml',iounit, codeparam_smart%parameters_value)

  ! DEFINE LOCAL DATABASE (TO BE FIRST CREATED WITH "imasdb ITER")
  pulse         = standalone_in%pulse
  run_in        = standalone_in%run_in
  run_out       = standalone_in%run_out
  input_db      = standalone_in%input_db
  input_machine = standalone_in%input_machine
  local_db      = standalone_in%local_db
  local_machine = standalone_in%local_machine
  
  ! USERNAME DEFINED BY ENVIRONMENT VARIABLE USERNAME
  if(local_db.eq.'user_environment_variable') then
     call getenv('USER',local_db)
  endif

  ! OPEN INPUT DATAFILE FROM OFFICIAL IMAS SCENARIO DATABASE
  write(*,*) '=> Read input IDSs'
  !call imas_open_env('ids',pulse,run_in,idx,input_db,input_machine,'3')
  call al_build_uri_from_legacy_parameters(MDSPLUS_BACKEND, pulse, run_in, input_db, input_machine, "3", "", uri, status)
  call al_begin_dataentry_action(uri, OPEN_PULSE, idx, status);
  call ids_get(idx,'equilibrium',equilibrium_in)
  call ids_get(idx,'core_profiles',core_profiles_in)
  call imas_close(idx)
  write(*,*) 'Finished reading input IDSs'


  ! EXECUTE PELLET ABLATION MODEL: SMART
  call smart(equilibrium_in,core_profiles_in,pellets_in,core_profiles_out,codeparam_smart, &
       error_flag,error_message)

  if(error_flag.eq.0) then
     ! EXPORT RESULTS TO LOCAL DATABASE
     write(*,*) '=> Export output IDSs to local database'
     call al_build_uri_from_legacy_parameters(HDF5_BACKEND, pulse, run_out, local_db, local_machine, "3", "", uri, status)
     call al_begin_dataentry_action(uri, FORCE_CREATE_PULSE, idx, status);
     call ids_put(idx,'core_profiles',core_profiles_out)
     call imas_close(idx)
     write(*,*) 'Done exporting.'
     write(*,*) ' '
     write(*,*) 'End of standalone'
  else
     write(*,*) error_message
     write(*,*) '=> Program stopped.'
  endif

end program standalone

