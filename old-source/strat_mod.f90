module Strat
  use Types
  use Control
  use NewTrans
  use Initialization
  use Geometry
  use Diagnostics
  private

  type,public :: layer
    integer :: nx,ny,nz
    real(kind=r8) :: s,phi0,phi1,phi2
    real(kind=r8),dimension(:,:),pointer :: psi,q,k1,k2,k3,dx,dy,ddx,ddy,del2,res,rot
    real(kind=r8) :: tote,ave,enstrophy,vmax,umax,ske,deltapsi,ximax
    logical :: get_energy,need_vmax,get_enstrophy
    real(kind=r8),dimension(:,:),pointer :: jac_dx,jac_dy
    real(kind=r8),dimension(:,:),pointer :: use_q,jac_result,jac_temp
    real(kind=r8),dimension(:),pointer :: evec
  end type layer
  
  real(kind=r8),public :: gprime,rd,gamma,omega,bft,mindelta
  real(kind=r8),dimension(:,:),pointer,public :: psibc,psibt
  
  integer,public :: nmodes
  integer,public :: efnum,skenum,ensnum,force_mode
  
  logical,public :: reduced_gravity,mode_force,layer_force
  logical,public :: do_force,do_nl,do_rot,do_visc
  
  character(len=200),public :: energyfile,plotname
  type(layer),dimension(:),pointer,public :: slab,meanslab
  
  ! mean state
  real(kind=r8),public :: meant,dtold
  
  ! run constants
  real(kind=r8),public :: rhozero,ndt
  
  ! strat model constants
  real(kind=r8),public :: sdelta, nbsq, beta_st, rh, tau,basemax
  
  ! triple interaction coefficients
  real(kind=r8),dimension(:,:,:),allocatable,public :: xi
  
  
  public :: init_strat_run
  public :: layer_step_strat_split
  public :: strat_jac,do_strat_constants
  public :: save_slab,restore_slab,read_slab,checkstable_slab,save_mean
  public :: makemode
  
contains
  subroutine init_strat_run()
    integer :: istart,i
    real(kind=r8) :: time
  
    call get_run_info()
    call init_nondim()
  
    nx = param("nx",128,0)
    ny = param("ny",128,0)
    nz = param("nz",5000,0)
    nmodes = param("nmodes",2,0)
    reduced_gravity = param("reduced_gravity",.false.,0)
    beta_st = param("beta_st",0.9_r8,0)
  
  ! set up storage
    allocate(slab(nmodes),meanslab(nmodes))
    do i=1,nmodes
       call makemode(slab(i),nx,ny,nz)
       call makemode(meanslab(i),nx,ny,nz)
    enddo
  
  ! set up constants
    call do_strat_constants(slab)
  
  ! set up for mean state
    meant = 0.0_r8
    do i=1,nmodes
       meanslab(i)%psi = 0.0_r8
       meanslab(i)%q = 0.0_r8
       meanslab(i)%evec = slab(i)%evec
       meanslab(i)%s = slab(i)%s
    enddo
  
  
    time = param("t",0.0_r8,0)
    ngyres = param("ngyres",2,0)
  !  tau = param("tau",0.04_r8,0)
  
  ! grids are used to initialize the plots.
  
    call init_geometry(nnx=nx,nny=ny,chebx=.true.,cheby=.true.)
  
  ! set up functional pieces
    call init_trans(nx,ny)
  ! this sets the shape of the force, but not the amplitude
    call init_force(ngyres)
  ! rescale forcing here if necessary
    forcing = forcing
    call tcxcy(forcing) ! I'd like to put this in equations, but ho ho
    call init_integrator()
  
    mindelta = minval(deltas(1,:))
  
  
  ! now initialize the slab.  either start from scratch or read a file.
  ! this is set up to make other possibilities (other initial conditions)
  ! available.  I am not implementing anything fancy at the moment.
       istart = param("start",1,0)
       select case(istart)
       case(0) ! restore from a saved slab
          print *,"reading old field"
          call restore_slab(slab)
          print *, "restored"
       case default   ! start from scratch
          print *,"starting from scratch"
          ! we can get away with zero on the second layer since it
          ! doesn't affect the timestep.
          do i=1,nmodes
             slab(i)%q = forcing * 1.0e-8_r8
          enddo
       end select
  
  ! write out initial q
    if (istart >= 0) then
  
  ! assume the next run will continue
       call modparam("start", 0)
  
  ! write param file for next run
       call writeparams()
    endif
  
  ! open a file for the energy
    write(unit=energyfile,fmt=*) nrun
    energyfile = trim(adjustl(basename))//"."//trim(adjustl(energyfile))&
                    //".energy"
    efnum = 97
    open(file=energyfile,unit=efnum,action="write",status="replace")
  
    skenum = 98
    open(file=trim(adjustl(energyfile))//"-ske",unit=skenum,&
         action="write",status="replace")
  
    ensnum = 99
    open(file=trim(adjustl(energyfile))//"-ens",unit=ensnum,&
         action="write",status="replace")
  
  
    return
  end subroutine init_strat_run
  subroutine layer_step_strat_split(n,t)
    integer,intent(in) :: n
    real(kind=r8),intent(inout) :: t
    integer :: i
  
  
    do i=1,nmodes
      slab(i)%use_q => slab(i)%q
      slab(i)%jac_result => slab(i)%k1
      slab(i)%jac_temp => slab(i)%res
    enddo
    
    call strat_jac(slab,xi)
    
    do i=1,nmodes
      ! generate del2psi
      slab(i)%ddx = slab(i)%psi
      slab(i)%ddy = slab(i)%psi
      call d_dx_cheb(slab(i)%ddx)
      call d_dx_cheb(slab(i)%ddx)
      call d_dy_cheb(slab(i)%ddy)
      call d_dy_cheb(slab(i)%ddy)
      slab(i)%del2 = alpha**2*slab(i)%ddx+slab(i)%ddy
      
    ! no force
      slab(i)%k1 = slab(i)%k1 - (a_2ds-a_4dm3*slab(i)%s)*slab(i)%del2
    enddo
    
    if(do_force) then
    if (mode_force) then
      slab(force_mode)%k1 = slab(force_mode)%k1 + forcing
    else if(layer_force) then
      do i=1,nmodes
        slab(i)%k1 = slab(i)%k1 + forcing * slab(i)%phi2
      enddo
    else
    ! bc force
      do i=1,nmodes
         slab(i)%k1 = slab(i)%k1 + forcing * slab(i)%phi0
      enddo
    endif
    endif
    
    do i=1,nmodes
       slab(i)%k3 = slab(i)%q + dt3*slab(i)%k1
    end do
    
    do i=1,nmodes
      call solve_p_nondim(slab(i)%k3,slab(i)%psi,slab(i)%s,alpha)
    enddo
    
    
    do i=1,nmodes
      slab(i)%use_q => slab(i)%k3
      slab(i)%jac_result => slab(i)%k2
      slab(i)%jac_temp => slab(i)%res
    enddo
    
    call strat_jac(slab,xi)
    
    do i=1,nmodes
      ! generate del2psi
      slab(i)%ddx = slab(i)%psi
      slab(i)%ddy = slab(i)%psi
      call d_dx_cheb(slab(i)%ddx)
      call d_dx_cheb(slab(i)%ddx)
      call d_dy_cheb(slab(i)%ddy)
      call d_dy_cheb(slab(i)%ddy)
      slab(i)%del2 = alpha**2*slab(i)%ddx+slab(i)%ddy
    
    ! no force
      slab(i)%k2 = slab(i)%k2 - (a_2ds-a_4dm3*slab(i)%s)*slab(i)%del2
    enddo
    
    
    if(do_force) then
    if(mode_force) then
      slab(force_mode)%k2 = slab(force_mode)%k2 + forcing
    else if (layer_force) then
      do i=1,nmodes
        slab(i)%k2 = slab(i)%k2 + forcing * slab(i)%phi2
      enddo
    else
    ! bc force
      do i=1,nmodes
         slab(i)%k2 = slab(i)%k2 + forcing * slab(i)%phi0
      enddo
    endif ! mode_force
    endif ! do_force
    
    do i=1,nmodes
       slab(i)%k2 = slab(i)%q + dt23*slab(i)%k2
    end do
    
    do i=1,nmodes
      call solve_p_nondim(slab(i)%k2,slab(i)%psi,slab(i)%s,alpha)
    enddo
    
    
    do i=1,nmodes
      slab(i)%use_q => slab(i)%k2
      slab(i)%jac_result => slab(i)%k3
      slab(i)%jac_temp => slab(i)%res
    enddo
    
    call strat_jac(slab,xi)
    
    
    do i=1,nmodes
       ! generate del2psi
       slab(i)%ddx = slab(i)%psi
       slab(i)%ddy = slab(i)%psi
       call d_dx_cheb(slab(i)%ddx)
       call d_dx_cheb(slab(i)%ddx)
       call d_dy_cheb(slab(i)%ddy)
       call d_dy_cheb(slab(i)%ddy)
       slab(i)%del2 = alpha**2*slab(i)%ddx+slab(i)%ddy
       
    
    ! no force
      slab(i)%k3 = slab(i)%k3 - (a_2ds-a_4dm3*slab(i)%s)*slab(i)%del2
    enddo
    
    if(do_force) then
    if(mode_force) then
      slab(force_mode)%k3 = slab(force_mode)%k3+ forcing
    else if(layer_force) then
      do i=1,nmodes
        slab(i)%k3 = slab(i)%k3 + forcing * slab(i)%phi2
      enddo
    else
    ! bc force
      do i=1,nmodes
         slab(i)%k3 = slab(i)%k3 + forcing * slab(i)%phi0
      enddo
    endif ! mode_force
    endif ! do_force
    
    ! commenting out the last part of this line removes the explicit
    ! terms from the equation.  We need the rest of the computation
    ! to avoid division by zero problems.
    
    do i=1,nmodes
       slab(i)%q = slab(i)%q + dt4*(slab(i)%k1+3*slab(i)%k3)
    end do
    
    do i=1,nmodes
      call solve_p_nondim(slab(i)%q,slab(i)%psi,slab(i)%s,alpha)
    enddo
    
    
  
  if(do_visc) then
    do i=1,nmodes
       slab(i)%ddx = slab(i)%q
       slab(i)%ddy = slab(i)%q
       call d_dx_cheb(slab(i)%ddx)
       call d_dx_cheb(slab(i)%ddx)
       call d_dy_cheb(slab(i)%ddy)
       call d_dy_cheb(slab(i)%ddy)
       slab(i)%del2 = -(alpha**2*slab(i)%ddx+slab(i)%ddy+ndt*slab(i)%q)
    
       call solve_p_nondim(slab(i)%del2,slab(i)%q,ndt,alpha)
    enddo
    
  endif ! do_visc
  
    do i=1,nmodes
      call solve_p_nondim(slab(i)%q,slab(i)%psi,slab(i)%s,alpha)
    enddo
    
  
    ! update time, timestep
    t = t+dt
    call checkstable_slab()
    
  
    return
  end subroutine layer_step_strat_split
  
  ! this follows the nondimensionalization of Ghil et. al. 2002
  subroutine strat_jac(sl,xi)
  
  ! can we trim any of this?
  
      type(layer),dimension(:),intent(inout) :: sl
      real(kind=r8),dimension(:,:,:),intent(in) :: xi
  
      integer :: i,j,k,nmodes,index
  
      nmodes = size(sl)
  
      do k = 1, nmodes
         sl(k)%jac_result = 0.0_r8
         sl(k)%jac_dx = 0.0_r8
         sl(k)%jac_dy = 0.0_r8
         sl(k)%dx = sl(k)%psi
         sl(k)%dy = sl(k)%psi
  
         call d_dx_cheb(sl(k)%dx)
         call d_dy_cheb(sl(k)%dy)
  
         sl(k)%rot = -sl(k)%dx
  
      ! to real space for the nonlinear part
         call fcxcy(sl(k)%use_q)
         call fcxcy(sl(k)%dy)
         call fcxcy(sl(k)%dx)
      enddo
  
      do i=1,nmodes
            ! for stability
            if (sl(i)%need_vmax) then
        ! here for now; move for efficiency when debugged.
               basemax = sl(i)%phi0*6.3/minval(deltas(1,:))
               sl(:)%vmax = basemax*sl(:)%ximax
               do index=0,ny
                  sl(i)%umax = sl(i)%ximax*(maxval(abs(sl(i)%dy(:,index))/deltas(1,:))+basemax)
                  if(sl(i)%umax > sl(i)%vmax) sl(i)%vmax = sl(i)%umax
               enddo
               do index=0,nx
                  sl(i)%umax = sl(i)%ximax*maxval(abs(sl(i)%dx(index,:)/deltas(2,:)))
                  if(sl(i)%umax > sl(i)%vmax) sl(i)%vmax = sl(i)%umax
               enddo
        
               if (sl(i)%vmax > 1.0e16_r8) then
        !          print *, "blowing up.  stopping at step",curstep," time ",t
        !          stop ! head off blowups
               endif
        
        ! this is handled by basemax above.
            ! for initial forcing/wave resolution:
        !       if (sl(i)%vmax < beta_st/mindelta) sl(i)%vmax = beta_st/mindelta
        !       if (sl(i)%vmax < 1.0_r8) sl(i)%vmax = 1.0_r8
        
               sl(i)%need_vmax = .false.
            endif
        
      enddo
  
      do k = 1, nmodes
        do i=1,nmodes
           do j=1,nmodes
  
              ! xi contains the nondimensional factors
              if(abs(xi(i,j,k)) > 1.0e-8_r8) then
  
  !              print *, "doing nonlinear for k,i,j:",k,i,j,xi(i,j,k)
                !    if(sl(i)%get_energy) then
                !      call compute_energy(sl(i)%dx,sl(i)%dy,sl(i)%tote,sl(i)%ave,sl(i)%ske,sl(i)%deltapsi)
                !      sl(i)%get_energy = .false.
                !    endif
                
                !    if(sl(j)%get_enstrophy) then
                ! I am not using this at the moment.
                !      call compute_enstrophy(sl(j)%use_q,sl(j)%enstrophy)
                !      sl(j)%get_enstrophy = .false.
                !    endif
                
  
                sl(k)%jac_dx = sl(k)%jac_dx+(sl(i)%dx*sl(j)%use_q)*xi(i,j,k)
                sl(k)%jac_dy = sl(k)%jac_dy+(sl(i)%dy*sl(j)%use_q)*xi(i,j,k)
              endif
           enddo
        enddo
      enddo
  
  
  
  
  
  
  
      do k=1,nmodes
        call tcxcy(sl(k)%jac_dx)
        call tcxcy(sl(k)%jac_dy)
        call tcxcy(sl(k)%use_q)
  
        call d_dx_cheb(sl(k)%jac_dy)
        call d_dy_cheb(sl(k)%jac_dx)
  
  
        if (do_nl) then
           if(do_rot) then ! Full problem
              sl(k)%jac_result = sl(k)%jac_dy - sl(k)%jac_dx + sl(k)%rot
           else ! nonlinearity, but without rotation:
              sl(k)%jac_result = sl(k)%jac_dy - sl(k)%jac_dx
           endif ! do_rot
         else  ! linear problem only:
           if(do_rot) then
              sl(k)%jac_result = sl(k)%rot
           endif
  ! with no nonlinearity and no rotation, we do nothing.
         endif
  
  
  
        ! the better dealiasing
        ! this is handled in the transforms
        where(truncate)
           sl(k)%jac_result = 0.0_r8
        end where
  
      end do ! k
  
     return
  end subroutine strat_jac
  subroutine do_strat_constants(slab)
    use shoot
    type(layer),dimension(:),intent(inout) :: slab
  
    real(kind=r8),dimension(:),allocatable :: evals
    real(kind=r8),dimension(:,:),allocatable :: evecs,prev,cur,next
    real(kind=r8) :: dlambda,sfac,forcingdelta,hfrac
    integer :: nz,i,n,nzstar
  
    nz = slab(1)%nz
    dlambda = 0.1_r8
  
    allocate(evals(nmodes),evecs(nmodes,nz),prev(4,nz),cur(4,nz),next(4,nz))
    allocate(xi(nmodes,nmodes,nmodes))
  
    rhozero = param("rhozero",1021.0_r8,0)
    dt = param("dt",0.00001_r8,0)
    nbsq = param("nbsq",7.0e-5_r8,0)
    fzero = param("fzero", 1.0e-4_r8,0)
    ly = param("ly", 1.0e6_r8,0)
    sdelta = param("sdelta", 0.5_r8,0)
    h = param("h",4000.0_r8,0)
    rh = param("rh", 0.4_r8,0)
    mode_force = param("mode_force", .false., 0)
    layer_force=param("layer_force",.false., 0)
    force_mode = param("force_mode", 1, 0)
    forcingdelta = param("forcingdelta",0.125_r8,0)
    do_force = param("do_force", .true., 0)
    do_rot = param("do_rot", .true., 0)
    do_nl = param("do_nl", .true., 0)
    do_visc = param("do_visc", .true., 0)
  
    dtold = dt
  
  
    sfac = nbsq*h**2/(fzero**2*ly**2)
    print *, "sfac is ",sfac
  
    cnconst = 2*alpha**4/(dt*deltam**3)
    ndt = cnconst
  
  
  ! setup z
  ! this is offset by half a step since phi' is offset
    do i=1,nz
       prev(4,i) = (i-0.5_r8)/(nz-1.0_r8)
    enddo
  
  ! set S
    do i=1,nz
       cur(4,i) = sfac/(prev(4,i)**(2.0_r8*sdelta))
       next(4,i) = cur(4,i)
    enddo
  
  ! ladder up the eigenmodes and eigenvalues
  
    !initial phi, phi'
    next(1,:) = 1.0_r8
    cur(1,:) = 1.0_r8
    next(2,:) = 0.0_r8
    cur(2,:) = 0.0_r8
  
  ! lambda = 0 is an eigenvalue
    cur(3,:) = 0.0_r8
    next(3,:) = 0.0_r8
  
    call fire(cur,next)
  
    evals(1)=0.0_r8
    evecs(1,:) = next(1,:)
  
  ! go on to the next one.
    prev = cur
    cur = next
    next(3,:) = cur(3,:) + 0.01_r8
    call fire(cur,next)
  
    prev = cur
    cur = next
  
  
  
    evals(1) = 0.0_r8
     n = 1
  !print *, "before ladder"
    do i=2,nmodes
  !print *, i
       call ladder(prev,cur,next,evals,evecs,n,dlambda)
    end do
  !print *, "after ladder"
  
    call normalize(evecs)
  !print *, "after normalize"
    call triples(evecs,xi)
  !print *, "after triples"
  
    ! put nondimensional factors in xi
    xi(:,:,:) = xi(:,:,:) * a_3di2
  
    do i=1,nmodes
       slab(i)%ximax = maxval(xi(i,:,:))
    enddo
  
  ! evals are lambda NOT squared.
    print *, "evals:"
    print *, evals(:)
  
  ! to get effective deformation radii, take ly/sqrt(evals(:))
    print *, "effective L_d (km):"
    print *, ly/evals(2:)
  
  
    ! slab(i)%s should have lambda_i^2; eval(i) is lambda_i
    slab(:)%s = evals(:)**2
    slab(:)%phi0 = evecs(:,1)
    slab(:)%phi1 = evecs(:,nz)
  
    hfrac = forcingdelta/(1+forcingdelta)
    nzstar = hfrac*nz
    print *, "forcingdelta, hfrac, nzstar:",forcingdelta,hfrac,nzstar
    do i=1,nmodes
      slab(i)%phi2 = integrate(evecs(i,1:nzstar),0.0_r8,hfrac)
    enddo
  
  
  !  do i=1,nmodes
  !    print *, "phi2 for mode",i,slab(i)%phi2,slab(i)%phi0
  !     slab(i)%evec(:) = evecs(i,:)
  !  enddo
  
  ! this is probably the wrong thing to do without rescaling something else.
  ! or maybe not doing it without rescaling deltai is wrong...
  ! rescale phi0, phi2
  !  slab(:)%phi0 = slab(:)%phi0/(sum(slab(:)%phi0))
  !  slab(:)%phi2 = slab(:)%phi2/(sum(slab(:)%phi2))
  
  
    do i=1,nmodes
      print *, "phi2 for mode",i,slab(i)%phi2,slab(i)%phi1,slab(i)%phi0
       slab(i)%evec(:) = evecs(i,:)
    enddo
  
    deallocate(evals,evecs)
  
    return
  end subroutine do_strat_constants
  subroutine save_slab(slab,n)
    type(layer),dimension(:),intent(in) :: slab
    integer,intent(in) :: n
    character(len=200) :: filename,nstr
    integer :: i
  
    write(unit=filename,fmt=*) nrun
    write(unit=nstr,fmt=*) n
    call modparam("nslab",n)
    filename = trim(adjustl(basename))//"."//trim(adjustl(filename))&
               //"."//trim(adjustl(nstr))//".slab"
  !  print *, "file is called ",filename
  
    open(unit=17,file=filename,status="replace",action="write",&
          form="unformatted")
  
    do i=1,nmodes
       write(unit=17) slab(i)%psi
       write(unit=17) slab(i)%q
       write(unit=17) slab(i)%evec
       write(unit=17) slab(i)%s
    enddo
  
    close(unit=17)
  
    return
  end subroutine save_slab
  subroutine restore_slab(slab)
    type(layer),dimension(:),intent(inout) :: slab
    character(len=200) :: filename,nstr
    integer :: i,nslab
    nslab = param("nslab",0,0)
    write(unit=filename,fmt=*) nrun-1
    write(unit=nstr,fmt=*) nslab
    filename = trim(adjustl(basename))//"."//trim(adjustl(filename))&
               //"."//trim(adjustl(nstr))//".slab"
  
    print *, "loading from file ",filename
  
  
    open(unit=13,file=filename,status="old",action="read",&
          form="unformatted",position="rewind")
    do i=1,nmodes
       read(unit=13) slab(i)%psi
       read(unit=13) slab(i)%q
  ! this would need to change for pmc
       read(unit=13) slab(i)%evec
       read(unit=13) slab(i)%s
    enddo
  
    close(unit=13)
  
  ! get mean also
    write(unit=filename,fmt=*) nrun-1
    write(unit=nstr,fmt=*) nslab
  
    filename = trim(adjustl(basename))//"-mean."//trim(adjustl(filename))&
               //"."//trim(adjustl(nstr))//".slab"
  
    print *, "loading mean from file ",filename
    open(unit=13,file=filename,status="old",action="read",&
          form="unformatted",position="rewind")
    do i=1,nmodes
       read(unit=13) smeanslab(i)%psi
       read(unit=13) smeanslab(i)%q
  ! this would need to change for pmc
       read(unit=13) smeanslab(i)%evec
       read(unit=13) smeanslab(i)%s
  !  set up mean u, v
       smeanslab(i)%dx = smeanslab(i)%psi
       smeanslab(i)%dy = smeanslab(i)%psi
       call d_dx_cheb(smeanslab(i)%dx)
       call d_dy_cheb(smeanslab(i)%dy)
       call fcxcy(smeanslab(i)%dx)
       call fcxcy(smeanslab(i)%dy)
    enddo
  
    close(unit=13)
  
    return
  end subroutine restore_slab
  subroutine read_slab(slab,filename)
    type(layer),dimension(:),intent(inout) :: slab
    character(len=*),intent(in) :: filename
    integer :: i
  
    open(unit=13,file=filename,status="old",action="read",&
          form="unformatted",position="rewind")
    do i=1,nmodes
       read(unit=13) slab(i)%psi
       read(unit=13) slab(i)%q
       read(unit=13) slab(i)%evec
       read(unit=13) slab(i)%s
    enddo
  
    close(unit=13)
  
    return
  end subroutine read_slab
  subroutine save_mean(slab,n,time)
    type(layer),dimension(:),intent(inout) :: slab
    integer,intent(in) :: n
    real(kind=r8),intent(in) :: time
    character(len=200) :: filename,nstr
    integer :: i
  
    write(unit=filename,fmt=*) nrun
    write(unit=nstr,fmt=*) n
  
    filename = trim(adjustl(basename))//"-mean."//trim(adjustl(filename))&
               //"."//trim(adjustl(nstr))//".slab"
  !  print *, "file is called ",filename
  
    open(unit=17,file=filename,status="replace",action="write",&
          form="unformatted")
    do i=1,nmodes
       slab(i)%dx = slab(i)%psi/time
       slab(i)%dy = slab(i)%q/time
  
       write(unit=17) slab(i)%dx
       write(unit=17) slab(i)%dy
       write(unit=17) slab(i)%evec
       write(unit=17) slab(i)%s
    enddo
  
    close(unit=17)
  
    return
  end subroutine save_mean
  subroutine checkstable_slab()
     real(kind=r8) :: cflnew,vm
  
     !get it again, next time step.
     slab(:)%need_vmax = .true.
  
     vm = maxval(slab(:)%vmax)
  !   print *, "vmax:", slab(1)%vmax,slab(2)%vmax,vm,cfl
  
  
     dtold = dt
  !return without changing dt--removes the adaptive timestep(removal removed)
     ! if we're within 1%, return without changing dt.
     cflnew = dt*vm
     if(abs(cfl-cflnew)/cfl < 0.01) then
  !      print *, "returning without changing dt"
      return
     endif
  
        print *, "setting dt from ",dt,"to",cfl/vm
  !   print *, "POSSIBLE DANGER: ",vm
     dt = cfl/vm
  
  ! if dt is too small, we won't get anywhere. bail out
     if(dt <= 1.0e-7_r8) then
        print *, "dt too small:",dt,vm
        print *, "bailing out"
        stop
     endif
  
  
     dt3 = dt/3.0_r8
     dt4 = dt/4.0_r8
     dt23= 2.0_r8*dt/3.0_r8
     cnconst = 2*alpha**4/(dt*deltam**3)
     ndt = cnconst
  
     return
  end subroutine checkstable_slab
  
  subroutine makemode(l,nx,ny,nz)
    type(layer),intent(out) :: l
    integer,intent(in) :: nx,ny,nz
  
    l%nx = nx
    l%ny = ny
    l%nz = nz
  
    allocate(l%psi(0:nx+1,0:ny+1),l%q(0:nx+1,0:ny+1))
    allocate(l%k1(0:nx+1,0:ny+1),l%k2(0:nx+1,0:ny+1))
    allocate(l%k3(0:nx+1,0:ny+1))
    allocate(l%dy(0:nx+1,0:ny+1),l%ddy(0:nx+1,0:ny+1))
    allocate(l%dx(0:nx+1,0:ny+1),l%ddx(0:nx+1,0:ny+1))
    allocate(l%del2(0:nx+1,0:ny+1),l%res(0:nx+1,0:ny+1))
    allocate(l%evec(nz))
    allocate(l%jac_dy(0:nx+1,0:ny+1),l%jac_dx(0:nx+1,0:ny+1))
    allocate(l%rot(0:nx+1,0:ny+1))
  
    l%need_vmax = .true.
    l%get_energy = .true.
    l%get_enstrophy = .true.
  
    return
  end subroutine makemode
  
  
end module Strat
