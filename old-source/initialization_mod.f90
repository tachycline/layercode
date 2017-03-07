module Initialization
  use Types
  use Geometry
  use Control
  implicit none
  private

   real(kind=r8),dimension(:,:),allocatable,save,public :: forcing
   integer,public :: ngyres
   real(kind=r8),public :: asymmetry
   logical, public :: biharmonic
   
      ! nondimensional parameters
      real(kind=r8),public :: alpha,deltai,deltam,a_3di2
      real(kind=r8),public :: delta_s,a_2ds,a_4dm3
      ! dimensionl parameters
      real(kind=r8),public :: kappa,beta,fzero,w,thick1,famp,r
      ! scaling constants
      real(kind=r8),public :: lx,ly
   
      logical,public :: halflayer
      real(kind=r8),public :: ld,ldsq
   

   public :: init_force

  public :: init_nondim,init_integrator
! time/timestep variables
real(kind=r8),public :: t,dt,dt3,dt23,dt4,cfl,cnconst,tdat,st,tmax,t_dimen,&
     tfac,tframe
integer,public :: nsteps,totalsteps,nframe,curstep,startsteps,nrg,nst

logical,dimension(:,:),allocatable,save,public:: truncate

contains
  subroutine init_nondim()
    real(kind=r8) :: third

    lx = param("lx",1.0e6_r8,0)
    ly = param("ly",2.0e6_r8,0)
    beta = param("beta",2.0e-11_r8,0)
    kappa = param("kappa",1.0e3_r8,0)
    fzero = param("fzero",1.0e-4_r8,0)
    w = param("w",1.0e-6_r8,0)
    thick1 = param("thick1",1.0e3_r8,0)
    r = param("r",1.0e-10_r8,0)
    asymmetry = param("asymmetry",0.0_r8,0)
    biharmonic = param("biharmonic", .true., 0)

    famp = fzero*w/thick1

    third = 1.0_r8/3.0_r8

    alpha = ly/lx
    deltai = (famp/(beta*beta*lx*lx))**(1.0_r8/2.0_r8)

    ! for some reason, this breaks in F (linux)
    !   deltam = ((kappa/beta)**third)/lx
    ! do it this way instead.
    deltam = exp(third*log(kappa/beta))/lx

    delta_s =r/(beta*lx)


    call modparam("deltai",deltai)
    call modparam("deltam",deltam)
    call modparam("delta_s",delta_s)
    call modparam("alpha",alpha)

    a_3di2 = deltai**2/(alpha**3)
    a_2ds = delta_s/(alpha**2)
    a_4dm3 = deltam**3/(alpha**4)

    return
  end subroutine init_nondim
  subroutine init_integrator()

    ! parameters
    !   nsteps = param("nsteps",1000,0)
    tframe = param("tframe",10.0_r8,0)  ! tframe in years
    startsteps = param("totalsteps",0,0)
    cfl = param("cfl",0.5_r8,0)
    dt = param("dt",1.0e-4_r8,0)
    t = param("t",0.0_r8,0)
    nrg = param("nrg",100,0)
    tdat = param("tdat",1.0_r8,0) ! tdat in years
    tmax = param("tmax",20.0_r8,0)! tmax in years

    lx = param("lx",5.0e5_r8,0)
    ly = param("ly",1.0e6_r8,0)
    beta = param("beta",2.0e-11_r8,0)

    ! convert tmax and tdat to nondimensional time units.

    tfac = 3600*24*365*beta*ly**2/lx

    tmax = tmax*tfac
    tframe = tframe*tfac
    tdat = tdat*tfac

    nsteps = tmax/tdat
    nframe = tframe/tdat

    dt3 = dt/3.0_r8
    dt4 = dt/4.0_r8
    dt23= 2.0_r8*dt/3.0_r8

    if(allocated(truncate)) then
       deallocate(truncate)
    end if
    allocate(truncate(0:nx+1,0:ny+1))

    ! cheb modes are numbered 0 .. nx
    ! fourier modes are numbered 0 .. ny/2,-ny/2-1 .. -1
    ! have to account for this.

    !do j=0,ny+1
    !   do i=0,nx+1
    !      ! this makes anything but a square domain a waste
    !      if ((i*i + ky2(j)/(PI**2)) >= (2*min(nx,ny)/3)**2) then
    !         truncate(i,j) = .true.
    !      else
    !         truncate(i,j) = .false.
    !      endif
    !   enddo
    !enddo

    truncate(:,:) = .false.
    truncate(:,2*ny/3-2:) = .true.
    truncate(2*nx/3-2:,:) = .true.

    cnconst = 2*alpha**4/(dt*deltam**3)

    return
  end subroutine init_integrator

  subroutine init_force(gyres)
    integer,intent(in) :: gyres
    ! this needs to be called AFTER geometry is initialized.
    integer ::i,j
    logical :: btforcefixed
    real(kind=r8) :: thick2

    thick1 = param("thick1",8.0e2_r8,0)
    thick2 = param("thick2",3.2e3_r8,0)
    btforcefixed = param("btforcefixed",.false.,0)

    if(allocated(forcing)) then
       deallocate(forcing)
    end if
    allocate(forcing(0:nx+1,0:ny+1))
    forcing = 0.0_r8

    if(abs(asymmetry) > 1.0_r8) asymmetry = 0.0_r8
    print *, "asymmetry: ",asymmetry
    do i=0,nx
       do j=0,ny

          forcing(i,j) = -sin((gyres/2.0_r8)*PI*(grid(2,j)+1.0_r8)) &
               + asymmetry*cos((gyres/4.0_r8)*PI*(grid(2,j)+1.0_r8))
          !             forcing(i,j) = exp(-((grid(1,i)/0.95_r8)**50))*forcing(i,j)
       enddo
    enddo
    if (btforcefixed) then
       forcing = forcing * (thick1+thick2)/thick2
    endif
    return
  end subroutine init_force
end module Initialization
