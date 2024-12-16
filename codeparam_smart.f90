module mod_codeparam_smart

  implicit none
  type type_smart_data
     ! STEPUP
     double precision:: TAU 
     double precision:: dtau
     double precision:: QNB 
     double precision:: F01B
     double precision:: F02B
     double precision:: F03B
     double precision:: GN2E
     double precision:: GN2I
     double precision:: sw_smart
     double precision:: sw_ech2a
     ! SMART
     double precision:: YAM
     double precision:: YVP
     double precision:: YVOL
     double precision:: YCOS0
     double precision:: YEFF
     double precision:: YDL
     double precision:: yswitch
     double precision:: yglength
     ! ECH2a
     double precision:: ROCEC
     double precision:: ROCDR
     double precision:: QECR
     double precision:: YEFFec
  end type type_smart_data

contains

  subroutine assign_codeparam(string, smart)

    use xml2eg_mdl, only: xml2eg_parse_memory, xml2eg_get, type_xml2eg_document, xml2eg_free_doc, set_verbose
    
    ! Input/Output
    character(len=132), pointer :: string(:)
    type(type_smart_data), intent(out) :: smart

    ! Internal
    type(type_xml2eg_document) :: doc

    ! Parse the "string". This means that the data is put into a document "doc"
    call xml2eg_parse_memory(string, doc)
    call set_verbose(.TRUE.) ! Only needed if you want to see what's going on in the parsing

    call xml2eg_get(doc,'TAU',smart%TAU)
    call xml2eg_get(doc,'dtau',smart%dtau)
    call xml2eg_get(doc,'QNB',smart%QNB)
    call xml2eg_get(doc,'F01B',smart%F01B)
    call xml2eg_get(doc,'F02B',smart%F02B)
    call xml2eg_get(doc,'F03B',smart%F03B)
    call xml2eg_get(doc,'GN2E',smart%GN2E)
    call xml2eg_get(doc,'GN2I',smart%GN2I)
    call xml2eg_get(doc,'sw_smart',smart%sw_smart)
    call xml2eg_get(doc,'sw_ech2a',smart%sw_ech2a)

    call xml2eg_get(doc,'YAM',smart%YAM)
    call xml2eg_get(doc,'YVP',smart%YVP)
    call xml2eg_get(doc,'YVOL',smart%YVOL)
    call xml2eg_get(doc,'YCOS0',smart%YCOS0)
    call xml2eg_get(doc,'YEFF',smart%YEFF)
    call xml2eg_get(doc,'YDL',smart%YDL)
    call xml2eg_get(doc,'yswitch',smart%yswitch)
    call xml2eg_get(doc,'yglength',smart%yglength)

    call xml2eg_get(doc,'ROCEC',smart%ROCEC)
    call xml2eg_get(doc,'ROCDR',smart%ROCDR)
    call xml2eg_get(doc,'QECR',smart%QECR)
    call xml2eg_get(doc,'YEFFec',smart%YEFFec)

    ! Make sure to clean up after you!!
    ! When calling "xml2eg_parse_memory" memory was allocated in the "doc" object.
    ! This memory is freed by "xml2eg_free_doc(doc)"
    call xml2eg_free_doc(doc)

  end subroutine assign_codeparam
end module mod_codeparam_smart
