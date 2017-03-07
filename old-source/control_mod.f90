module Control
  use Types
!  use f90_unix_proc
  private

  character(len=40),public :: basename
  character(len=40),private :: runfile,paramfile
  character(len=40),public :: paramarchive
  integer,public :: nrun

  integer, dimension(:),allocatable,private :: ipvals
  real(kind=r8), dimension(:),allocatable,private :: dpvals
  character(len=64), dimension(:),allocatable,private :: ipnames,dpnames

  logical, dimension(:),allocatable,private :: lpvals
  character(len=64), dimension(:),allocatable,private :: cpvals
  character(len=64), dimension(:),allocatable,private :: cpnames,lpnames

!  public :: basename,nrun
!  private :: runfile, paramfile, paramarchive

  interface param
     module procedure iparam,dparam,lparam,cparam
  end interface

  interface addparam
     module procedure addiparam,adddparam,addlparam,addcparam
  end interface

  interface modparam
     module procedure modiparam,moddparam,modlparam,modcparam
  end interface

  public :: param,addparam,modparam
  private :: iparam,dparam,addiparam,adddparam,modiparam,moddparam
  private :: lparam,cparam,addcparam,addlparam,modlparam,modcparam
!  private :: ipnames,dpnames,dpvals,ipvals,lpnames,lpvals,cpnames,cpvals
  public :: get_run_info

  public :: writeparams,purgeparams

  public :: spawn_run, get_meta_run_info

contains

  subroutine get_run_info()
    character(len=1) :: c2, typei, typed, typel, typec
    character(len=160) :: line,filler
    character(len=40) :: name
    integer,dimension(3) :: blanks
    integer :: ios,ivalue,namelen
    real(kind=r8) :: dvalue
    logical :: lvalue
    character(len=64) :: cvalue

    ! open the run file and get the base name for the run.
    runfile = "run"
    open(unit=29,file=runfile,form="formatted",action="read",&
                    status="old",position="rewind")
    read(unit=29,fmt="(a)") basename
    close(unit=29)

    ! get the name of the param file
    paramfile = trim(adjustl(basename)) //".param"

    ! open the param file and read the params
    open(unit=29,file=paramfile,form="formatted",action="read",&
             status="old",position="rewind")

    typei = "i"
    typed = "d"
    typel = "l"
    typec = "c"

    ! read a line of the form "run #"
    read(unit=29, fmt="(a)") line
    line = trim(adjustl(line))
    read(unit=line, fmt="(a4,i5)") filler,nrun

    write(unit=paramarchive,fmt=*) nrun
    paramarchive = trim(adjustl(basename))//"."//trim(adjustl(paramarchive))&
         // ".param"
    open(unit=31,file=paramarchive,form="formatted",action="write",&
         status="replace")
    write(unit=31,fmt="(a)") trim(adjustl(line))

    ! now read the rest of the file.
    ! use iostat to bail out of the loop when we get eof
    readfile: do
       read(unit=29,fmt="(a)",iostat=ios) line
       if(ios < 0) exit readfile
       line = trim(adjustl(line))
       write(unit=31,fmt="(a)") trim(adjustl(line))
       ! find the lengths of the words
       blanks(1) = scan(line," ")
       blanks(2) = blanks(1)+scan(line(blanks(1)+1:)," ")
       blanks(3) = blanks(2)+scan(line(blanks(2)+1:)," ")
       ! length of the name is blanks(2) - blanks(1)
       namelen = blanks(2)-blanks(1)

       read(unit=line,fmt="(a1)") c2
       name = "                                        "
       read(unit=line(blanks(1):blanks(2)),fmt="(a)") name
       select case(c2)
       case("i")
          read(unit=line(blanks(2):),fmt="(i10)") ivalue
          call addparam(name,ivalue)
       case("d")
          read(unit=line(blanks(2):),fmt="(f15.5)") dvalue
          call addparam(name,dvalue)
       case("l")
          read(unit=line(blanks(2):),fmt=*) lvalue
          call addparam(name,lvalue)
       case("c")
          read(unit=line(blanks(2):),fmt="(a)") cvalue
!          print *, "name is ", name, " and value is ",cvalue
          call addparam(name,cvalue)
       end select
    enddo readfile
    close(unit=29)
    close(unit=31)

    return
  end subroutine get_run_info
  !-------------------------------------------------------
  subroutine addiparam(name, ival)
    character(len=*),intent(in) :: name
    integer,intent(in) :: ival
    integer :: niparams
    integer,dimension(:),allocatable :: itemp
    character(len=64),dimension(:),allocatable :: ntemp

    if (allocated(ipvals)) then
       niparams = size(ipvals)+1
       allocate(itemp(niparams),ntemp(niparams))
       itemp(1:niparams-1) = ipvals
       ntemp(1:niparams-1) = ipnames
       deallocate(ipnames,ipvals)
       itemp(niparams) = ival
       ntemp(niparams) = trim(name)
       allocate(ipvals(niparams),ipnames(niparams))
       ipvals = itemp
       ipnames = ntemp
       deallocate(itemp,ntemp)
    else
       allocate(ipvals(1),ipnames(1))
       ipvals(1) = ival
       ipnames(1) = trim(name)
    end if

    return
  end subroutine addiparam
  !-------------------------------------------------------
  subroutine adddparam(name, dval)
    character(len=*),intent(in) ::  name
    real(kind=r8),intent(in) :: dval
    integer :: ndparams
    real(kind=r8),dimension(:),allocatable :: dtemp
    character(len=64),dimension(:),allocatable :: ntemp

    if (allocated(dpvals)) then
       ndparams = size(dpvals)+1
       allocate(dtemp(ndparams),ntemp(ndparams))
       dtemp(1:ndparams-1) = dpvals
       ntemp(1:ndparams-1) = dpnames
       deallocate(dpnames,dpvals)
       dtemp(ndparams) = dval
       ntemp(ndparams) = trim(name)
       allocate(dpvals(ndparams),dpnames(ndparams))
       dpvals = dtemp
       dpnames = ntemp
       deallocate(dtemp,ntemp)
    else
       allocate(dpvals(1),dpnames(1))
       dpvals(1) = dval
       dpnames(1) = trim(name)
    end if

    return
  end subroutine adddparam
  !-------------------------------------------------------
  subroutine addlparam(name, lval)
    character(len=*),intent(in) :: name
    logical,intent(in) :: lval
    integer :: nlparams
    logical,dimension(:),allocatable :: ltemp
    character(len=64),dimension(:),allocatable :: ntemp

    if (allocated(lpvals)) then
       nlparams = size(lpvals)+1
       allocate(ltemp(nlparams),ntemp(nlparams))
       ltemp(1:nlparams-1) = lpvals
       ntemp(1:nlparams-1) = lpnames
       deallocate(lpnames,lpvals)
       ltemp(nlparams) = lval
       ntemp(nlparams) = trim(name)
       allocate(lpvals(nlparams),lpnames(nlparams))
       lpvals = ltemp
       lpnames = ntemp
       deallocate(ltemp,ntemp)
    else
       allocate(lpvals(1),lpnames(1))
       lpvals(1) = lval
       lpnames(1) = trim(name)
    end if

    return
  end subroutine addlparam
  !-------------------------------------------------------
  subroutine addcparam(name, cval)
    character(len=*),intent(in) :: name
    character(len=*),intent(in) :: cval
    integer :: ncparams
    character(len=64),dimension(:),allocatable :: ltemp
    character(len=64),dimension(:),allocatable :: ntemp

    if (allocated(cpvals)) then
       ncparams = size(cpvals)+1
       allocate(ltemp(ncparams),ntemp(ncparams))
       ltemp(1:ncparams-1) = cpvals
       ntemp(1:ncparams-1) = cpnames
       deallocate(cpnames,cpvals)
       ltemp(ncparams) = trim(adjustl(cval))
       ntemp(ncparams) = trim(adjustl(name))
       allocate(cpvals(ncparams),cpnames(ncparams))
       cpvals = ltemp
       cpnames = ntemp
       deallocate(ltemp,ntemp)
    else
       allocate(cpvals(1),cpnames(1))
       cpvals(1) = trim(adjustl(cval))
       cpnames(1) = trim(adjustl(name))
    end if

    return
  end subroutine addcparam
  !----------------------------------------------------
  function dparam(name, defval, iuse) result(dval)
    real(kind=r8),intent(in) :: defval
    integer,intent(in) :: iuse
    character(len=*),intent(in):: name
    real(kind=r8) :: dval
    integer :: i

    do i=1,size(dpnames)
       if(trim(adjustl(name))==trim(adjustl(dpnames(i)))) then
          dval = dpvals(i)
          return
       end if
    end do

    ! not in the list yet.
    ! add to the list unless we shouldn't (iuse /= 0)
    if(iuse == 0) then
!       call adddparam(name,defval)
       dval = defval
       return
    endif

    ! if iuse /= 0, flag an error, tell user to fix param file.
    if (iuse /= 0) then
!       write(unit=*,fmt=*)"unable to find double parameter:"
!       write(unit=*,fmt=*) trim(adjustl(name))
!       stop
    endif
    return
  end function dparam
  !---------------------------------------------------------
  function iparam(name, idefval, iuse) result(iprm)
    integer :: iprm
    integer,intent(in):: idefval,iuse
    character(len=*),intent(in) :: name
    integer:: i

    do i=1,size(ipnames)
       if(trim(adjustl(name))==trim(adjustl(ipnames(i)))) then
          iprm = ipvals(i)
          return
       end if
    end do

    ! not in the list yet.
    ! add to the list unless we shouldn't (iuse /= 0)
    if(iuse == 0) then
!       call addiparam(name,idefval)
       iprm = idefval
       return
    endif

    ! if iuse /= 0, flag an error, tell user to fix param file.
    if (iuse /= 0) then
!       write(unit=*,fmt=*)"unable to find integer parameter:"
!       write(unit=*,fmt=*) trim(adjustl(name))
!       stop
    endif
    return
  end function iparam
  !----------------------------------------------------
  function lparam(name, ldefval, iuse) result(lval)
    logical,intent(in) :: ldefval
    integer,intent(in) :: iuse
    character(len=*),intent(in):: name
    logical :: lval
    integer :: i

    do i=1,size(lpnames)
       if(trim(adjustl(name))==trim(adjustl(lpnames(i)))) then
          lval = lpvals(i)
          return
       end if
    end do

    ! not in the list yet.
    ! add to the list unless we shouldn't (iuse /= 0)
    if(iuse == 0) then
!       call addlparam(name,defval)
       lval = ldefval
       return
    endif

    ! if iuse /= 0, flag an error, tell user to fix param file.
    if (iuse /= 0) then
!       write(unit=*,fmt=*)"unable to find logical parameter:"
!       write(unit=*,fmt=*) trim(adjustl(name))
!       stop
    endif
    return
  end function lparam
  !----------------------------------------------------
  function cparam(name, defval, iuse) result(dval)
    character(len=*),intent(in) :: defval
    integer,intent(in) :: iuse
    character(len=*),intent(in):: name
    character(len=64) :: dval
    integer :: i

    do i=1,size(cpnames)
       if(trim(adjustl(name))==trim(adjustl(cpnames(i)))) then
          dval = cpvals(i)
          return
       end if
    end do

    ! not in the list yet.
    ! add to the list unless we shouldn't (iuse /= 0)
    if(iuse == 0) then
!       call addcparam(name,defval)
       dval = defval
       return
    endif

    ! if iuse /= 0, flag an error, tell user to fix param file.
    if (iuse /= 0) then
!       write(unit=*,fmt=*)"unable to find character parameter:"
!       write(unit=*,fmt=*) trim(adjustl(name))
!       stop
    endif
    return
  end function cparam
  !-------------------------------------------------------------------
  subroutine modiparam(name, ival)
    character(len=*),intent(in) :: name
    integer,intent(in) :: ival
    integer :: i

    do i=1,size(ipnames)
       if(trim(adjustl(name))==trim(adjustl(ipnames(i)))) then
          ipvals(i) = ival
          return
       end if
    end do

    call addiparam(name, ival)

    return
  end subroutine modiparam
  !-------------------------------------------------------------------
  subroutine modcparam(name, ival)
    character(len=*),intent(in) :: name
    character(len=*),intent(in) :: ival
    integer :: i

    do i=1,size(cpnames)
       if(trim(adjustl(name))==trim(adjustl(cpnames(i)))) then
          cpvals(i) = ival
          return
       end if
    end do

    call addcparam(name, ival)

    return
  end subroutine modcparam
  !------------------------------------------------------------------------
  subroutine moddparam(name, dval)
    character(len=*),intent(in) :: name
    real(kind=r8),intent(in) :: dval
    integer :: i

    do i=1,size(dpnames)
       if(trim(adjustl(name))==trim(adjustl(dpnames(i)))) then
          dpvals(i) = dval
          return
       end if
    end do

    call adddparam(name, dval)

    return
  end subroutine moddparam
  !------------------------------------------------------------------------
  subroutine modlparam(name, dval)
    character(len=*),intent(in) :: name
    logical,intent(in) :: dval
    integer :: i

    do i=1,size(lpnames)
       if(trim(adjustl(name))==trim(adjustl(lpnames(i)))) then
          lpvals(i) = dval
          return
       end if
    end do

    call addlparam(name, dval)

    return
  end subroutine modlparam
!----------------------------------------------------
  subroutine writeparams()
    integer :: i
    character(len=200) :: line

    ! open the file
    print *, "in writeparams, opening paramfile ",trim(adjustl(paramfile))
    open(unit=59, file=paramfile, action="write",status="old",&
         form="formatted",position="rewind")

    ! give the new run information
    write(unit=line,fmt="(a4,i5)") "run",nrun+1
    write(unit=59,fmt="(a)") trim(adjustl(line))

    ! do the integer params
    if(allocated(ipvals)) then
       do i=1,size(ipvals)
          write(unit=line,fmt="(i10)") ipvals(i)
          write(unit=line,fmt=*) "i "//trim(adjustl(ipnames(i)))//" "&
               //trim(adjustl(line))
          write(unit=59,fmt="(a)") trim(adjustl(line))
       enddo
    endif

    ! do the logical params
    if(allocated(lpvals)) then
       do i=1,size(lpvals)
          write(unit=line,fmt=*) lpvals(i)
          write(unit=line,fmt=*) "l "//trim(adjustl(lpnames(i)))//" "&
               //trim(adjustl(line))
          write(unit=59,fmt="(a)") trim(adjustl(line))
       enddo
    endif

    ! do the character params
    if(allocated(cpvals)) then
       do i=1,size(cpvals)
          write(unit=line,fmt="(a)") trim(adjustl(cpvals(i)))
          write(unit=line,fmt="(a)") "c "//trim(adjustl(cpnames(i)))//" "&
               //trim(adjustl(line))
          write(unit=59,fmt="(a)") trim(adjustl(line))
    enddo
    endif

    ! do the real params
    if(allocated(dpvals)) then
       do i=1,size(dpvals)
          write(unit=line,fmt="(es15.8)") dpvals(i)
          write(unit=line,fmt=*) "d "//trim(adjustl(dpnames(i)))//" "&
               //trim(adjustl(line))
          write(unit=59,fmt="(a)") trim(adjustl(line))
       enddo
    endif

    close(unit=59)

    return
  end subroutine writeparams
!-----------------------------------------------
  subroutine purgeparams()

    if(allocated(ipvals)) deallocate(ipvals,ipnames)
    if(allocated(dpvals)) deallocate(dpvals,dpnames)
    if(allocated(cpvals)) deallocate(cpvals,cpnames)
    if(allocated(lpvals)) deallocate(lpvals,lpnames)

    return
  end subroutine purgeparams
!-----------------------------------------------
  subroutine get_meta_run_info()
    character(len=1) :: c2, typei, typed
    character(len=80) :: line
    character(len=40) :: name
    integer,dimension(3) :: blanks
    integer :: ios,ivalue,namelen
    real(kind=r8) :: dvalue

    ! open the run file and get the base name for the run.
    runfile = "metarun"
    open(unit=29,file=runfile,form="formatted",action="read",&
          status="old",position="rewind")
    read(unit=29,fmt="(a)") basename
    close(unit=29)

    ! get the name of the param file
    paramfile = trim(adjustl(basename)) //".meta"

    ! open the param file and read the params
    open(unit=29,file=paramfile,form="formatted",action="read",status="old")

    typei = "i"
    typed = "d"

    ! now read the rest of the file.
    ! use iostat to bail out of the loop when we get eof
    readfile: do
       read(unit=29,fmt="(a)",iostat=ios) line
       if(ios < 0) exit readfile
       line = trim(adjustl(line))
       ! find the lengths of the words
       blanks(1) = scan(line," ")
       blanks(2) = blanks(1)+scan(line(blanks(1)+1:)," ")
       blanks(3) = blanks(2)+scan(line(blanks(2)+1:)," ")
       ! length of the name is blanks(2) - blanks(1)
       namelen = blanks(2)-blanks(1)

       read(unit=line,fmt="(a1)") c2
       read(unit=line(blanks(1):blanks(2)),fmt="(a)") name
       select case(c2)
       case("i")
          read(unit=line(blanks(2):),fmt="(i10)") ivalue
          call addparam(name,ivalue)
       case("d")
          read(unit=line(blanks(2):),fmt="(f15.5)") dvalue
          call addparam(name,dvalue)
       end select
    enddo readfile
    close(unit=29)

    return
  end subroutine get_meta_run_info
!----------------------------------------------
  subroutine spawn_run()
    character(len=200) :: command,rundir

    ! assume basename, parameters have been set

    ! write the run file
    open(unit=11,file="run",action="write",status="old")
    write(unit=11,fmt="(a)") trim(adjustl(basename))
    close(unit=11)

    ! set the value of paramfile
    paramfile = trim(adjustl(basename)) // ".param"

    ! write the params
    nrun = -1
    call modparam("start",1)

    call writeparams()

    ! create the run directory
    rundir = "run-" // trim(adjustl(basename))
    command = "mkdir "// trim(adjustl(rundir))
    call system(command)

    ! move the files
    command = "mv run " // trim(adjustl(paramfile)) &
         // " " // trim(adjustl(rundir))
    call system(command)

    ! build spawn command
    command = "(cd " // trim(adjustl(rundir)) //"; " &
         // "/atmos01/cmckay/ns2dsrc/f90/inittest >" // &
         trim(adjustl(basename)) //".out)"

    ! spawn
    call system(command)

    return
  end subroutine spawn_run

end module control
