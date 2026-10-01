!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Read ini-style parameter file
!!
!! @copyright
!!   Copyright 2013-2023 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! modified: 2026 Masashi Ogiso mogiso@mri-jma.go.jp
!! ----
module m_readini

  !! -- Dependency
  use nrtype, only : fp

  !! -- Declarations
  implicit none
  private
  save

  !! -- Public Procedures
  public :: readini
  public :: readini__strict_mode
  logical :: strict_mode = .false. 

  !--

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! find keyword 'key' from ini-style file, and returns var
  !!
  !!
  !! @par Usage
  !!   call readini( io, key, var, def )
  !!     - io    (integer)   ! file io number
  !!     - key   (character) ! keyword
  !!     - var (any type)    ! output parameter
  !!     - def (any type)    ! default value
  !<
  !! --
  interface readini

    module procedure readini_f, readini_i, readini_c, readini_l

  end interface readini
  !! --------------------------------------------------------------------------------------------------------------------------- !!

contains


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine readini_c( io, key, var, def )

    !! -- Arguments
    integer,       intent(in)           :: io
    character(*),  intent(in)           :: key
    character(*),  intent(in), optional :: def
    character(*),  intent(out)          :: var

    character(256) :: keyword
    integer        :: ierr
    character(256) :: cline
    integer        :: keylen
    logical        :: isopen
    !! ----

    !! file status
    inquire( io, OPENED=isopen )

    if( .not. isopen ) then

      if(present(def)) then
        write(0,'(A)') 'ERROR [readini]: file not open.'
        var = trim(def)
        return
      else
        write(0, '(a)') "Error, [readini]: file not open, default not given."
        error stop
      endif

    end if



    !! initialize file I/O location
    rewind(io)

    keyword = trim(adjustl(key))
    keylen  = len_trim(keyword)


    do

      !! get one line
      read(io,'(A)', iostat=ierr) cline

      !! reach to the last line
      if( ierr /= 0 ) then
        write(0, '(3a)') "key " // trim(keyword) // " is not found." 
        if( .not. strict_mode) then
          if(present(def)) then
            write(0, '(3a)') "Use default value " // trim(def) // " instead."
            var = trim(def)
          else
            write(0, '(a)') 'Program terminate ... '
            error stop
          end if
        else
          write(0, '(a)') 'Program terminate ... '
          error stop
        end if
        return
      end if

      cline = adjustl(cline)


      !! comment line
      if( cline(1:1) == '#' .or. cline(1:1) == '!' ) cycle


      !! find keyword
      if( cline(1:keylen) == trim(keyword) ) then

        cline = adjustl( cline( keylen+1: ) )

        if( cline(1:1) == '=' ) then
          cline = adjustl(cline(2:))
          read(cline,*) var
          exit
        end if
      end if
    end do

    rewind(io)

    !! expand environmental variable
    call system__expenv( var )

  end subroutine readini_c
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine readini_f( io, key, var, def )

    integer,      intent(in)           :: io
    character(*), intent(in)           :: key
    real(fp),     intent(in), optional :: def
    real(fp),     intent(out)          :: var
    !! -
    character(256) :: avar, adef
    !! ----
    if(present(def)) then
      write(adef,*) def
      call readini_c( io, key, avar, def = adef )
    else
      call readini_c( io, key, avar )
    endif
    read(avar,*) var

  end subroutine readini_f
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine readini_i( io, key, var, def )
    integer,      intent(in)           :: io
    character(*), intent(in)           :: key
    integer,      intent(in), optional :: def
    integer,      intent(out)          :: var
    !! --
    character(256) :: avar, adef
    !! ----
    if(present(def)) then
      write(adef,*) def
      call readini_c( io, key, avar, def = adef )
    else
      call readini_c( io, key, avar )
    endif
    read(avar,*) var

  end subroutine readini_i
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine readini_l( io, key, var, def )

    integer,      intent(in)           :: io
    character(*), intent(in)           :: key
    logical,      intent(in), optional :: def
    logical,      intent(out)          :: var
    !! --
    character(256) :: avar, adef
    !! ----

    if(present(def)) then
      write(adef,*) def
      call readini_c( io, key, avar, def = adef )
    else
      call readini_c( io, key, avar )
    endif
    read(avar,*) var

  end subroutine readini_l
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  subroutine readini__strict_mode( mode )
    logical, intent(in) :: mode

    strict_mode = mode
    
  end subroutine readini__strict_mode


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Expand shell environmental variables wrapped in ${...}
  !<
  !! --
  subroutine system__expenv( str )
    character(*), intent(inout) :: str
    character(256) :: str2
    integer :: ikey1, ikey2
    integer :: iptr
    character(256) :: str3

    iptr = 1
    str2 = ''
    do
      ikey1 = scan( str(iptr:), "${" ) + iptr - 1
      if( ikey1==iptr-1 ) exit

      ikey2 = scan( str(iptr:), "}" ) + iptr -1
      str2=trim(str2) // str(iptr:ikey1-1)
      call system__getenv( str(ikey1+2:ikey2-1), str3 )
      str2 = trim(str2) // trim(str3)
      iptr = ikey2+1

    end do
    str2 = trim(str2) // trim(str(iptr:))

    str = trim(str2)

  end subroutine system__expenv
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Obtain environmental variable "name".
  !<
  !! --
  subroutine system__getenv( name, value )

    !! -- Arguments
    character(*), intent(in)  :: name
    character(*), intent(out) :: value

    !! ----

    call get_environment_variable( name, value )

  end subroutine system__getenv


end module m_readini
!! ----------------------------------------------------------------------------------------------------------------------------- !!
