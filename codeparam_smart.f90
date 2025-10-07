module mod_codeparam_smart

  implicit none
  type type_smart_data
     ! STEPUP
     double precision :: TAU, dtau, QNB, F01B, F02B, F03B, GN2E, GN2I
     double precision :: sw_smart, sw_ech2a, sw_stdout
     ! SMART
     double precision :: YAM, YVP, YVOL, YCOS0, YEFF, YDL, yswitch, yglength
     ! ECH2a
     double precision :: ROCEC, ROCDR, QECR, YEFFec
     ! Control Keys (expected length = 15)
     integer, allocatable :: key4control(:)
  end type type_smart_data

contains

  subroutine assign_codeparam(string, smart)
    use xml2eg_mdl, only: xml2eg_parse_memory, xml2eg_get, xml2eg_getAllocatable, &
                          type_xml2eg_document, xml2eg_free_doc, set_verbose
    character(len=132), pointer :: string(:)
    type(type_smart_data), intent(out) :: smart
    type(type_xml2eg_document) :: doc
    logical :: errflag
    integer, allocatable :: tmp(:)
    integer :: nread, expected_keys

    expected_keys = 15

    call xml2eg_parse_memory(string, doc)
    call set_verbose(.True.)   ! suppress library verbosity by default

    ! Scalars
    call xml2eg_get(doc,'TAU',      smart%TAU)
    call xml2eg_get(doc,'dtau',     smart%dtau)
    call xml2eg_get(doc,'QNB',      smart%QNB)
    call xml2eg_get(doc,'F01B',     smart%F01B)
    call xml2eg_get(doc,'F02B',     smart%F02B)
    call xml2eg_get(doc,'F03B',     smart%F03B)
    call xml2eg_get(doc,'GN2E',     smart%GN2E)
    call xml2eg_get(doc,'GN2I',     smart%GN2I)
    call xml2eg_get(doc,'sw_smart', smart%sw_smart)
    call xml2eg_get(doc,'sw_ech2a', smart%sw_ech2a)
    call xml2eg_get(doc,'sw_stdout',smart%sw_stdout)

    call xml2eg_get(doc,'YAM',      smart%YAM)
    call xml2eg_get(doc,'YVP',      smart%YVP)
    call xml2eg_get(doc,'YVOL',     smart%YVOL)
    call xml2eg_get(doc,'YCOS0',    smart%YCOS0)
    call xml2eg_get(doc,'YEFF',     smart%YEFF)
    call xml2eg_get(doc,'YDL',      smart%YDL)
    call xml2eg_get(doc,'yswitch',  smart%yswitch)
    call xml2eg_get(doc,'yglength', smart%yglength)

    call xml2eg_get(doc,'ROCEC',    smart%ROCEC)
    call xml2eg_get(doc,'ROCDR',    smart%ROCDR)
    call xml2eg_get(doc,'QECR',     smart%QECR)
    call xml2eg_get(doc,'YEFFec',   smart%YEFFec)

    ! key4control list (space separated)
    if (allocated(smart%key4control)) deallocate(smart%key4control)
    call xml2eg_getAllocatable(doc,'key4control', smart%key4control, errflag)
    if (errflag) then
       write(6,*) 'Warning: key4control element missing or parse error – defaults applied'
       allocate(smart%key4control(expected_keys))
       smart%key4control = 0
    end if

    if (.not. allocated(smart%key4control)) then
       allocate(smart%key4control(expected_keys))
       smart%key4control = 0
    end if

    nread = size(smart%key4control)
    if (nread /= expected_keys) then
       allocate(tmp(expected_keys))
       tmp = 0
       if (nread > 0) tmp(1:min(nread,expected_keys)) = smart%key4control(1:min(nread,expected_keys))
       if (allocated(smart%key4control)) deallocate(smart%key4control)
       allocate(smart%key4control(expected_keys))
       smart%key4control = tmp
       deallocate(tmp)
       write(6,*) 'Notice: key4control adjusted to length', expected_keys, '(input length=', nread, ')'
    end if

    call xml2eg_free_doc(doc)
  end subroutine assign_codeparam

  subroutine deallocate_smart_data(smart)
    type(type_smart_data), intent(inout) :: smart
    if (allocated(smart%key4control)) deallocate(smart%key4control)
  end subroutine deallocate_smart_data

end module mod_codeparam_smart
