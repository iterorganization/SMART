module mod_codeparam_smart

  implicit none
  type type_smart_data
     double precision:: YAM
     double precision:: YVP
     double precision:: YVOL
     double precision:: YCOS0
     double precision:: YEFF
     double precision:: YDL
     double precision:: yswitch
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
    
    call xml2eg_get(doc,'YAM',smart%YAM)
    call xml2eg_get(doc,'YVP',smart%YVP)
    call xml2eg_get(doc,'YVOL',smart%YVOL)
    call xml2eg_get(doc,'YCOS0',smart%YCOS0)
    call xml2eg_get(doc,'YEFF',smart%YEFF)
    call xml2eg_get(doc,'YDL',smart%YDL)
    call xml2eg_get(doc,'yswitch',smart%yswitch)

    ! Make sure to clean up after you!!
    ! When calling "xml2eg_parse_memory" memory was allocated in the "doc" object.
    ! This memory is freed by "xml2eg_free_doc(doc)"
    call xml2eg_free_doc(doc)

  end subroutine assign_codeparam
end module mod_codeparam_smart
