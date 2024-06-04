module mod_codeparam_standalone

  implicit none
  type type_standalone_data
     integer :: pulse,run_in,run_out
     character(len=200):: input_db,input_machine,local_db,local_machine
  end type type_standalone_data

contains

  subroutine assign_codeparam(string, standalone)

    use xml2eg_mdl, only: xml2eg_parse_memory, xml2eg_get, type_xml2eg_document, xml2eg_free_doc, set_verbose
    
    ! Input/Output
    character(len=132), pointer :: string(:)
    type(type_standalone_data), intent(out) :: standalone

    ! Internal
    type(type_xml2eg_document) :: doc

    ! Parse the "string". This means that the data is put into a document "doc"
    call xml2eg_parse_memory(string, doc)
    call set_verbose(.TRUE.) ! Only needed if you want to see what's going on in the parsing
    
    call xml2eg_get(doc,'pulse',standalone%pulse)
    call xml2eg_get(doc,'run_in',standalone%run_in)
    call xml2eg_get(doc,'run_out',standalone%run_out)
    call xml2eg_get(doc,'input_db',standalone%input_db)
    call xml2eg_get(doc,'input_machine',standalone%input_machine)
    call xml2eg_get(doc,'local_db',standalone%local_db)
    call xml2eg_get(doc,'local_machine',standalone%local_machine)

    ! Make sure to clean up after you!!
    ! When calling "xml2eg_parse_memory" memory was allocated in the "doc" object.
    ! This memory is freed by "xml2eg_free_doc(doc)"
    call xml2eg_free_doc(doc)

  end subroutine assign_codeparam
end module mod_codeparam_standalone
