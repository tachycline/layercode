module Layers
  use Types
  use Control
  use NewTrans
  use Initialization
  use Geometry
  use Diagnostics
  implicit none
  private


  type,public :: layer
    integer :: nx,ny
    real(kind=r8) :: thickness,rho,deltarho,s,btfac,phi0,phi1
    real(kind=r8),dimension(:,:),pointer :: psi,q,k1,k2,k3,dx,dy,ddx,ddy,del2,res
    real(kind=r8) :: tote,ave,enstrophy,vmax,umax,ske,deltapsi
    logical :: get_energy,need_vmax,get_enstrophy
    real(kind=r8),dimension(:,:),pointer :: use_q,jac_result,jac_temp
  end type layer
  
  real(kind=r8),public :: gprime,rd,gamma,omega,bft
  real(kind=r8),dimension(:,:),pointer,public :: psibc,psibt
  
  integer,public :: nlayers,nmodes
  integer,public :: efnum,ekenum
  
  logical,public :: reduced_gravity
  
  ! for stochastic force
  logical,public :: stochastic
  real(kind=r8),public :: st_amp
  
  character(len=200),public :: energyfile,ekefile,plotname
  type(layer),dimension(:),pointer,public :: slab,meanslab,smeanslab
  
  ! mean state
  real(kind=r8),public :: meant,dtold
  
  ! run constants
  real(kind=r8),public :: rhozero,ndt
  
  ! barotropic/baroclinic energies
  real(kind=r8),public :: bte,bce,btave,bcave,btske,btdeltapsi,bcske,bcdeltapsi
  
  ! potential energies
  real(kind=r8),public :: ape,avape,sape,avsape
  
  real(kind=r8),public :: sdelta, nbsq
  
  real(kind=r8),dimension(:,:,:),allocatable,public :: xi
  
  ! energy output
  integer,public :: which_nrg
  logical,public :: saved_mean
  
  public :: setdeltas,makelayer
  public :: do_layer_constants,layer_step_2
  public :: save_slab,restore_slab,read_slab,save_mean
  public :: init_layer_run,solve_multi
  public :: checkstable_slab,layer_jac
  public :: test_multi
  public :: layer_step_rg
  public :: layer_step_bt
  public :: layer_step_bf
  public :: init_mode
  public :: layer_step_rg_split
  public :: layer_step_bf_split
  
contains
  subroutine makelayer(l,nx,ny,thickness,rho)
    type(layer),intent(out) :: l
    integer,intent(in) :: nx,ny
    real(kind=r8),intent(in) :: thickness,rho
  
    l%nx = nx
    l%ny = ny
    l%thickness = thickness
    l%rho = rho
  
    allocate(l%psi(0:nx+1,0:ny+1),l%q(0:nx+1,0:ny+1))
    allocate(l%k1(0:nx+1,0:ny+1),l%k2(0:nx+1,0:ny+1))
    allocate(l%k3(0:nx+1,0:ny+1))
    allocate(l%dy(0:nx+1,0:ny+1),l%ddy(0:nx+1,0:ny+1))
    allocate(l%dx(0:nx+1,0:ny+1),l%ddx(0:nx+1,0:ny+1))
    allocate(l%del2(0:nx+1,0:ny+1),l%res(0:nx+1,0:ny+1))
  
    l%need_vmax = .true.
    l%get_energy = .true.
    l%get_enstrophy = .true.
  
    return
  end subroutine makelayer
  
  subroutine setdeltas(lower,upper)
    type(layer),intent(inout) :: lower,upper
  
    lower%deltarho = lower%rho - upper%rho
    upper%deltarho = lower%deltarho
  
     return
  end subroutine setdeltas
  
  subroutine do_layer_constants(slab)
    type(layer),dimension(:),intent(inout) :: slab
  
    real(kind=r8) :: g
  
    rhozero = param("rhozero",1021.0_r8,0)
    dt = param("dt",0.00001_r8,0)
  
  ! for mean state
    dtold = dt
    meant = 0.0_r8
  
    g = 9.8_r8
  
    ndt = 2*alpha**4/(dt*deltam**3)
  
    if(nlayers == 2) then
      gamma = slab(1)%thickness/slab(2)%thickness
      gprime = abs(slab(2)%deltarho)*g/rhozero
      !!$rd = sqrt(gprime)/fzero*sqrt(slab(1)%thickness*slab(2)%thickness/&
      !!$                  (slab(1)%thickness+slab(2)%thickness))
      rd = sqrt(gprime*slab(1)%thickness)/fzero
      slab(1)%s = (ly/rd)**2
      print *, "lambdasq is ",slab(1)%s, rd, ly, ly/rd
      
      slab(1)%btfac = slab(1)%s*gamma/(1+gamma)
      slab(2)%btfac = slab(1)%s/(1+gamma)
      print *, "btfac:",slab(1)%btfac,slab(2)%btfac
      print *, "lambda^2:",slab(1)%s
      
    endif
    if(nlayers == 1) then
       if(reduced_gravity) then
          ! the baseline value is from mccalpin and haidvogel (1996)
          slab(1)%deltarho = param("deltarho",2.6_r8,0)
          slab(1)%thickness = param("thick1",1000.0_r8,0)
          slab(1)%s=fzero**2*rhozero/(g*slab(1)%thickness*slab(1)%deltarho)
          
          ! nondimensionalize now...
          slab(1)%s = slab(1)%s*ly**2
          print *, "lambda^2:",slab(1)%s
          
          
       endif
    endif
  
    return
  end subroutine do_layer_constants
  subroutine layer_jac(l)
      type(layer),intent(inout) :: l
      integer :: index
      l%jac_result(:,:) = 0.0_r8
      l%dx = l%psi
      l%dy = l%psi
      call d_dx_cheb(l%dx)
      call d_dy_cheb(l%dy)
  
      l%jac_result = -l%dx
  
      ! to real space for the nonlinear part
      call fcxcy(l%use_q)
      call fcxcy(l%dy)
      call fcxcy(l%dx)
  
      l%dy = l%dy * a_3di2
      l%dx = l%dx * a_3di2
  
      ! for stability
      if (l%need_vmax) then
  ! I am not completely satisfied with this.  It seems too restrictive.
         l%vmax = 0.0_r8
  !       l%vmax = maxval((abs(l%dx)+abs(l%dy))/(grid(1,1)-grid(1,0)))
  
         do index=0,ny
            l%umax = maxval(abs(l%dy(:,index)/deltas(1,:)))
            if(l%umax > l%vmax) l%vmax = l%umax
         enddo
         do index=0,nx
            l%umax = maxval(abs(l%dx(index,:)/deltas(2,:)))
            if(l%umax > l%vmax) l%vmax = l%umax
         enddo
  
         if (l%vmax > 1.0e16_r8) then
  !          print *, "blowing up.  stopping at step",curstep," time ",t
  !          stop ! head off blowups
         endif
  
      ! for initial forcing/wave resolution:
         if (l%vmax < 1.0_r8) l%vmax = 1.0_r8
  
         l%need_vmax = .false.
      endif
  
      l%dx = l%dx*l%use_q
      l%dy = l%dy*l%use_q
  
      call tcxcy(l%dx)
      call tcxcy(l%dy)
  ! see if this fixes stuff
      call tcxcy(l%use_q)
  
      call d_dx_cheb(l%dy)
      call d_dy_cheb(l%dx)
  
      l%jac_result = l%jac_result + l%dy - l%dx
  
      ! the better dealiasing
      ! this is handled in the transforms
      where(truncate)
         l%jac_result = 0.0_r8
      end where
  
     return
  end subroutine layer_jac
  subroutine layer_step_2(t)
  !  integer,intent(in) :: n
    real(kind=r8),intent(inout) :: t
    integer :: i
  
  !  i = n
  
    if (stochastic) then
      call random_number(st_amp)
      st_amp = st_amp*2.0_r8
      ! otherwise, st_amp is initialized to 1.
    end if
    
    do i=1,nlayers
      slab(i)%use_q => slab(i)%q
      slab(i)%jac_result => slab(i)%k1
    enddo
    
    do i=1,nlayers
         call layer_jac(slab(i))
    end do
    
    slab(1)%k1 = slab(1)%k1 + forcing*st_amp
    
    if (biharmonic) then
      ! if we're doing biharmonic psi dissipation
      ! get del2 psi1 - del2 psi2
      do i=1,nlayers
        ! generate del2psi
        slab(i)%ddx = slab(i)%psi
        slab(i)%ddy = slab(i)%psi
        call d_dx_cheb(slab(i)%ddx)
        call d_dx_cheb(slab(i)%ddx)
        call d_dy_cheb(slab(i)%ddy)
        call d_dy_cheb(slab(i)%ddy)
        slab(i)%del2 = alpha**2*slab(i)%ddx+slab(i)%ddy
        
      enddo
    
      slab(2)%del2 = slab(1)%del2 - slab(2)%del2
      slab(1)%del2 = -slab(2)%btfac*slab(2)%del2
      slab(2)%del2 = slab(1)%btfac*slab(2)%del2
    
      do i=1,nlayers
         slab(i)%k1 = slab(i)%k1 - a_4dm3*slab(i)%del2
         slab(i)%k3 = slab(i)%q + dt3*slab(i)%k1
      end do
    else
      do i=1,nlayers
         slab(i)%k3 = slab(i)%q + dt3*slab(i)%k1
      end do
    end if
    
    
    
    slab(1)%ddx = slab(1)%btfac*slab(1)%k3 + slab(2)%btfac*slab(2)%k3
    call solve_p_nondim(slab(1)%ddx,psibt,0.0_r8,alpha)
    
    slab(2)%ddx = slab(1)%k3-slab(2)%k3
    call solve_p_nondim(slab(2)%ddx,psibc,slab(1)%s,alpha)
    
    slab(1)%psi = psibt/slab(1)%s + psibc*slab(2)%btfac/slab(1)%s
    slab(2)%psi = psibt/slab(1)%s - psibc*slab(1)%btfac/slab(1)%s
    
    
    do i=1,nlayers
      slab(i)%use_q => slab(i)%k3
      slab(i)%jac_result => slab(i)%k2
    enddo
    
    do i=1,nlayers
         call layer_jac(slab(i))
    end do
    
    slab(1)%k2 = slab(1)%k2+forcing*st_amp
    if (biharmonic) then
      do i=1,nlayers
         ! generate del2psi
         slab(i)%ddx = slab(i)%psi
         slab(i)%ddy = slab(i)%psi
         call d_dx_cheb(slab(i)%ddx)
         call d_dx_cheb(slab(i)%ddx)
         call d_dy_cheb(slab(i)%ddy)
         call d_dy_cheb(slab(i)%ddy)
         slab(i)%del2 = alpha**2*slab(i)%ddx+slab(i)%ddy
         
      enddo
    
      slab(2)%del2 = slab(1)%del2 - slab(2)%del2
      slab(1)%del2 = -slab(2)%btfac*slab(2)%del2
      slab(2)%del2 = slab(1)%btfac*slab(2)%del2
    
      do i=1,nlayers
         slab(i)%k2 = slab(i)%k2 - a_4dm3*slab(i)%del2
         slab(i)%k2 = slab(i)%q + dt23*slab(i)%k2
      end do
    else
      do i=1,nlayers
         slab(i)%k2 = slab(i)%q + dt23*slab(i)%k2
      end do
    end if
    
    slab(1)%ddx = slab(1)%btfac*slab(1)%k2 + slab(2)%btfac*slab(2)%k2
    call solve_p_nondim(slab(1)%ddx,psibt,0.0_r8,alpha)
    
    slab(2)%ddx = slab(1)%k2-slab(2)%k2
    call solve_p_nondim(slab(2)%ddx,psibc,slab(1)%s,alpha)
    
    slab(1)%psi = psibt/slab(1)%s + psibc*slab(2)%btfac/slab(1)%s
    slab(2)%psi = psibt/slab(1)%s - psibc*slab(1)%btfac/slab(1)%s
    
    
    do i=1,nlayers
      slab(i)%use_q => slab(i)%k2
      slab(i)%jac_result => slab(i)%k3
    enddo
    
    do i=1,nlayers
         call layer_jac(slab(i))
    end do
    
    slab(1)%k3 = slab(1)%k3+forcing*st_amp
    if(biharmonic) then
      do i=1,nlayers
         ! generate del2psi
         slab(i)%ddx = slab(i)%psi
         slab(i)%ddy = slab(i)%psi
         call d_dx_cheb(slab(i)%ddx)
         call d_dx_cheb(slab(i)%ddx)
         call d_dy_cheb(slab(i)%ddy)
         call d_dy_cheb(slab(i)%ddy)
         slab(i)%del2 = alpha**2*slab(i)%ddx+slab(i)%ddy
         
      enddo
    
      slab(2)%del2 = slab(1)%del2 - slab(2)%del2
      slab(1)%del2 = -slab(2)%btfac*slab(2)%del2
      slab(2)%del2 = slab(1)%btfac*slab(2)%del2
    
    
    ! commenting out the last part of this line removes the explicit
    ! terms from the equation.  We need the rest of the computation
    ! to avoid division by zero problems.
    
      do i=1,nlayers
         slab(i)%k3 = slab(i)%k3 - a_4dm3*slab(i)%del2
         slab(i)%q = slab(i)%q + dt4*(slab(i)%k1+3*slab(i)%k3)
      end do
    else
      do i=1,nlayers
         slab(i)%q = slab(i)%q + dt4*(slab(i)%k1+3*slab(i)%k3)
      end do
    end if
    
    
    slab(1)%ddx = slab(1)%btfac*slab(1)%q + slab(2)%btfac*slab(2)%q
    call solve_p_nondim(slab(1)%ddx,psibt,0.0_r8,alpha)
    
    slab(2)%ddx = slab(1)%q-slab(2)%q
    call solve_p_nondim(slab(2)%ddx,psibc,slab(1)%s,alpha)
    
    slab(1)%psi = psibt/slab(1)%s + psibc*slab(2)%btfac/slab(1)%s
    slab(2)%psi = psibt/slab(1)%s - psibc*slab(1)%btfac/slab(1)%s
    
    
  
    do i=1,nlayers
       slab(i)%ddx = slab(i)%q
       slab(i)%ddy = slab(i)%q
       call d_dx_cheb(slab(i)%ddx)
       call d_dx_cheb(slab(i)%ddx)
       call d_dy_cheb(slab(i)%ddy)
       call d_dy_cheb(slab(i)%ddy)
       slab(i)%del2 = -(alpha**2*slab(i)%ddx+slab(i)%ddy+ndt*slab(i)%q)
    
       call solve_p_nondim(slab(i)%del2,slab(i)%q,ndt,alpha)
    enddo
    
  
    slab(1)%ddx = slab(1)%btfac*slab(1)%q + slab(2)%btfac*slab(2)%q
    call solve_p_nondim(slab(1)%ddx,psibt,0.0_r8,alpha)
    
    slab(2)%ddx = slab(1)%q-slab(2)%q
    call solve_p_nondim(slab(2)%ddx,psibc,slab(1)%s,alpha)
    
    slab(1)%psi = psibt/slab(1)%s + psibc*slab(2)%btfac/slab(1)%s
    slab(2)%psi = psibt/slab(1)%s - psibc*slab(1)%btfac/slab(1)%s
    
  
    ! update time, timestep
    t = t+dt
    call checkstable_slab()
    
  
    return
  end subroutine layer_step_2
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
  
    do i=1,nlayers
       write(unit=17) slab(i)%psi
       write(unit=17) slab(i)%q
       write(unit=17) slab(i)%thickness,slab(i)%rho,slab(i)%s
       write(unit=17) slab(i)%deltarho
       write(unit=17) slab(i)%btfac
    enddo
  
    close(unit=17)
  
    return
  end subroutine save_slab
  subroutine restore_slab(slab)
    type(layer),dimension(:),intent(inout) :: slab
    character(len=200) :: filename,nstr
    integer :: i,nslab
    real(kind=r8) :: thickness,rho,s,deltarho,btfac
    nslab = param("nslab",0,0)
    write(unit=filename,fmt=*) nrun-1
    write(unit=nstr,fmt=*) nslab
    filename = trim(adjustl(basename))//"."//trim(adjustl(filename))&
               //"."//trim(adjustl(nstr))//".slab"
  
    print *, "loading from file ",filename
  
  
    open(unit=13,file=filename,status="old",action="read",&
          form="unformatted",position="rewind")
    do i=1,nlayers
       read(unit=13) slab(i)%psi
       read(unit=13) slab(i)%q
  ! the rest of these should go with the values from the param file.
       read(unit=13) thickness,rho,s
       read(unit=13) deltarho
       read(unit=13) btfac
    enddo
  
    close(unit=13)
  ! get mean also
    write(unit=filename,fmt=*) nrun-1
    write(unit=nstr,fmt=*) nslab
  
    filename = trim(adjustl(basename))//"-mean."//trim(adjustl(filename))&
               //".slab"
  ! this is the old version
  !  filename = trim(adjustl(basename))//"-mean."//trim(adjustl(filename))&
  !             //"."//trim(adjustl(nstr))//".slab"
  
    print *, "loading mean from file ",filename
    open(unit=13,file=filename,status="old",action="read",&
          form="unformatted",position="rewind")
    do i=1,nlayers
       read(unit=13) smeanslab(i)%psi
       read(unit=13) smeanslab(i)%q
  ! this would need to change for pmc
       read(unit=13) thickness,rho,s
       read(unit=13) deltarho
       read(unit=13) btfac
  
       smeanslab(i)%ddx = smeanslab(i)%psi
       smeanslab(i)%ddy = smeanslab(i)%psi
       call d_dx_cheb(smeanslab(i)%ddx)
       call d_dx_cheb(smeanslab(i)%ddx)
       call d_dy_cheb(smeanslab(i)%ddy)
       call d_dy_cheb(smeanslab(i)%ddy)
       smeanslab(i)%del2 = alpha**2*smeanslab(i)%ddx + smeanslab(i)%ddy
       print *, "alpha is ",alpha
       call fcxcy(smeanslab(i)%del2)
       call fcxcy(smeanslab(i)%psi)
  
       call compute_energy(smeanslab(i)%del2,smeanslab(i)%psi,smeanslab(i)%tote,&
            smeanslab(i)%ave,smeanslab(i)%ske,smeanslab(i)%deltapsi)
  
       print *, "layer ",i,"mean energy:",smeanslab(i)%tote
    enddo
  
    if(nlayers == 2) then
    ! do mean bt/bc also. put bt in (1) and bc in (2)
    smeanslab(1)%dx = (slab(1)%btfac*smeanslab(1)%psi +&
              slab(2)%btfac*smeanslab(2)%psi)/slab(1)%s
       call tcxcy(smeanslab(1)%dx)
       smeanslab(1)%ddx = smeanslab(1)%dx
       smeanslab(1)%ddy = smeanslab(1)%dx
       call d_dx_cheb(smeanslab(1)%ddx)
       call d_dx_cheb(smeanslab(1)%ddx)
       call d_dy_cheb(smeanslab(1)%ddy)
       call d_dy_cheb(smeanslab(1)%ddy)
       smeanslab(1)%dy = alpha**2*smeanslab(1)%ddx + smeanslab(1)%ddy
       call fcxcy(smeanslab(1)%dy)
       call fcxcy(smeanslab(1)%dx)
  
    smeanslab(2)%dx = smeanslab(1)%psi - smeanslab(2)%psi
       call tcxcy(smeanslab(2)%dx)
       smeanslab(2)%ddx = smeanslab(2)%dx
       smeanslab(2)%ddy = smeanslab(2)%dx
       call d_dx_cheb(smeanslab(2)%ddx)
       call d_dx_cheb(smeanslab(2)%ddx)
       call d_dy_cheb(smeanslab(2)%ddy)
       call d_dy_cheb(smeanslab(2)%ddy)
       smeanslab(2)%dy = alpha**2*smeanslab(2)%ddx + smeanslab(2)%ddy
       call fcxcy(smeanslab(2)%dy)
       call fcxcy(smeanslab(2)%dx)
    endif
    close(unit=13)
  
    return
  end subroutine restore_slab
  subroutine read_slab(slab,filename)
    type(layer),dimension(:),intent(inout) :: slab
    character(len=*),intent(in) :: filename
    integer :: i
  
    open(unit=13,file=filename,status="old",action="read",&
          form="unformatted",position="rewind")
    do i=1,nlayers
       read(unit=13) slab(i)%psi
       read(unit=13) slab(i)%q
       read(unit=13) slab(i)%thickness,slab(i)%rho,slab(i)%s
       read(unit=13) slab(i)%deltarho
       read(unit=13) slab(i)%btfac
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
  
  !  filename = trim(adjustl(basename))//"-mean."//trim(adjustl(filename))&
  !             //"."//trim(adjustl(nstr))//".slab"
    filename = trim(adjustl(basename))//"-mean."//trim(adjustl(filename))&
               //".slab"
  !  print *, "file is called ",filename
  
    open(unit=17,file=filename,status="replace",action="write",&
          form="unformatted")
  
    do i=1,nlayers
       slab(i)%ddx = slab(i)%psi/time
       slab(i)%ddy = slab(i)%q/time
  
       write(unit=17) slab(i)%ddx
       write(unit=17) slab(i)%ddy
       write(unit=17) slab(i)%thickness,slab(i)%rho,slab(i)%s
       write(unit=17) slab(i)%deltarho
       write(unit=17) slab(i)%btfac
    enddo
  
    close(unit=17)
  
    write(unit=filename,fmt=*) nrun
  ! save variance
  !  filename = trim(adjustl(basename))//"-var."//trim(adjustl(filename))&
  !             //"."//trim(adjustl(nstr))//".slab"
    filename = trim(adjustl(basename))//"-var."//trim(adjustl(filename))&
               //".slab"
  !  print *, "file is called ",filename
  
    open(unit=17,file=filename,status="replace",action="write",&
          form="unformatted")
  
    do i=1,nlayers
       slab(i)%k1 = slab(i)%ddx
       slab(i)%k2 = slab(i)%ddy
       call fcxcy(slab(i)%k1)
       call fcxcy(slab(i)%k2)
  
       slab(i)%k1 = slab(i)%dx/time - slab(i)%k1**2
       slab(i)%k2 = slab(i)%dy/time - slab(i)%k2**2
  
       call tcxcy(slab(i)%k1)
       call tcxcy(slab(i)%k2)
  
       write(unit=17) slab(i)%k1
       write(unit=17) slab(i)%k2
       write(unit=17) slab(i)%thickness,slab(i)%rho,slab(i)%s
       write(unit=17) slab(i)%deltarho
       write(unit=17) slab(i)%btfac
    enddo
  
    close(unit=17)
  
  
    return
  end subroutine save_mean
  subroutine init_layer_run()
    integer :: istart,i
    real(kind=r8) :: time,thick,dens
    character(len=80) :: tstr,dstr
  
    call get_run_info()
    call init_nondim()
  
    nx = param("nx",128,0)
    ny = param("ny",128,0)
    nlayers = param("nlayers",2,0)
    reduced_gravity = param("reduced_gravity",.false.,0)
    biharmonic = param("biharmonic",.true.,0)
  
    stochastic = param("stochastic",.false.,0)
    st_amp = 1.0_r8  ! in case stochastic is false
  
    which_nrg = param("which_nrg",5,0)
    if (nlayers /= 2) then
       which_nrg = 1
       call modparam("which_nrg", 1)
    endif
  ! set up storage
    allocate(slab(nlayers),meanslab(nlayers),smeanslab(nlayers))
    do i=1,nlayers
       write(unit=dstr,fmt=*) i
       tstr = "thick"//trim(adjustl(dstr))
       dstr = "dens"//trim(adjustl(dstr))
       dens = param(dstr,1021.0_r8,0)
       thick = param(tstr,1000.0_r8,0)
       call makelayer(slab(i),nx,ny,thick,dens)
       call makelayer(meanslab(i),nx,ny,thick,dens)
       call makelayer(smeanslab(i),nx,ny,thick,dens)
    enddo
  
    do i=2,nlayers
       call setdeltas(slab(i),slab(i-1))
    enddo
  
  ! barotropic and baroclinic streamfunctions
  
    allocate(psibt(0:nx+1,0:ny+1),psibc(0:nx+1,0:ny+1))
  
  ! set up constants
    call do_layer_constants(slab)
  
    if (.not. biharmonic) then
     r = slab(1)%s*kappa/ly**2
     delta_s =r/(beta*lx)
     call modparam("r",r)
     call modparam("delta_s",delta_s)
     a_2ds = delta_s/(alpha**2)
     print *, "delta_s:",delta_s
    endif
  
  ! mean state
    meant = 0.0_r8
    do i=1,nlayers
       meanslab(i)%psi = 0.0_r8
       meanslab(i)%q = 0.0_r8
       meanslab(i)%dx = 0.0_r8
       meanslab(i)%dy = 0.0_r8
    enddo
  
  
    time = param("t",0.0_r8,0)
    ngyres = param("ngyres",2,0)
  ! grids are used to initialize the plots.
  
    call init_geometry(nnx=nx,nny=ny,chebx=.true.,cheby=.true.)
  
  ! set up functional pieces
    call init_trans(nx,ny)
    call init_force(ngyres)
    call init_integrator()
  
  ! now initialize the slab.  either start from scratch or read a file.
  ! this is set up to make other possibilities (other initial conditions)
  ! available.  I am not implementing anything fancy at the moment.
       istart = param("start",1,0)
       select case(istart)
       case(0) ! restore from a saved slab
          print *,"reading old field"
          call restore_slab(slab)
          print *, "restored"
          saved_mean = .true.
  
       case default   ! start from scratch
          print *,"starting from scratch"
          ! we can get away with zero on the second layer since it
          ! doesn't affect the timestep.
          do i=1,nlayers
             slab(i)%q = forcing * 1.0e-8_r8
          enddo
          saved_mean = .false.
  
       end select
  
    call tcxcy(forcing) ! I'd like to put this in equations, but ho ho
  
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
  
    ekenum = 99
    ekefile=trim(adjustl(energyfile))//".eke"
    open(file=ekefile,unit=ekenum,action="write",status="replace")
  
    energyfile = trim(adjustl(energyfile))//"-template"
    open(file=energyfile,unit=98,action="write",status="replace")
  
    select case(which_nrg)
      case(1)  ! single layer
         write(unit=98,fmt="(a)") "t, slab(1)%tote, slab(1)%enstrophy"
      case(2)  ! 2 layer, with everything
         write(unit=98,fmt="(a)") "t, slab(1)%tote, slab(1)%enstrophy, "&
              //"slab(2)%tote, slab(2)%enstrophy, bte, bce, slab(1)%ske, "&
              //"slab(2)%ske, btske, bcske, ape, sape"
      case(3)  ! 2 layer, without APE
         write(unit=98,fmt="(a)") "t, slab(1)%tote, slab(1)%enstrophy, "&
              //"slab(2)%tote, slab(2)%enstrophy, bte, bce, slab(1)%ske, "&
              //"slab(2)%ske, btske, bcske"
      case(4)  ! 2 layer, without bt/bc
         write(unit=98,fmt="(a)") "t, slab(1)%tote, slab(1)%enstrophy, "&
              //"slab(2)%tote, slab(2)%enstrophy, slab(1)%ske, "&
              //"slab(2)%ske, ape, sape"
      case(5)  ! 2 layer, without ape, bt/bc
         write(unit=98,fmt="(a)") "t, slab(1)%tote, slab(1)%enstrophy, "&
              //"slab(2)%tote, slab(2)%enstrophy, slab(1)%ske, "&
              //"slab(2)%ske"
      case(6)  ! 2 layer, with everything, including minmax
         write(unit=98,fmt="(a)") "t, slab(1)%tote, slab(1)%enstrophy, "&
              //"slab(2)%tote, slab(2)%enstrophy, bte, bce, slab(1)%ske, "&
              //"slab(2)%ske, btske, bcske, ape, sape, "&
              //"minmax1, minmax2, btminmax, bcminmax"
    end select
  
    close(unit=98)
  
  
  
  
  
    return
  end subroutine init_layer_run
  subroutine solve_multi(f,u,lambda,alpha)
    real(kind=r8),dimension(0:,0:),intent(in) :: f
    real(kind=r8),dimension(0:,0:),intent(out) :: u
    complex(kind=r8),dimension(:),intent(in) :: lambda
    real(kind=r8),intent(in) :: alpha
  
    complex(kind=r8),dimension(0:nxtrunc+2,&
                       0:nytrunc+2) :: uprime
  !  complex(kind=r8),dimension(0:nxtrunc+2,&
  !                     0:nytrunc+2) :: fprime
    real(kind=r8) :: a2
    integer :: i,j,k
  
  !  print *, "solve multi",size(lambda),alpha
  
    a2 = alpha**2
  !print *, "matmul"
    uprime = cmplx(f(0:nxtrunc+2,0:nytrunc+2),0.0_r8,r8)
  !print *, "matmul"
    do k=1,size(lambda)
  !print *, "matmul"
       uprime = matmul(pinv,matmul(uprime,trans_q))
  !     print *, "pass",k,"lambda:",lambda(k)
       do i=0,nxtrunc+2
          do j=0,nytrunc+2
             uprime(i,j) = uprime(i,j)/(a2*wra(i+1) + wrb(j+1) - lambda(k))
          enddo
       enddo
       uprime = matmul(p,matmul(uprime,trans_qinv))
  !     print *, "maxvals:",maxval(real(uprime)),maxval(aimag(uprime))
    enddo
  
  ! double check that im(u') is zero?
  ! what about instability/errors?
  
    u(0:nxtrunc+2,0:nytrunc+2) = real(uprime)
    u(nxtrunc-1:,:) = 0.0_r8
    u(:,nytrunc-1:) = 0.0_r8
  
  !  print *,"end solve multi"
  
    return
  end subroutine solve_multi
  subroutine checkstable_slab()
     real(kind=r8) :: cflnew,vm
  
     !get it again, next time step.
     slab(:)%need_vmax = .true.
  
     vm = maxval(slab(:)%vmax)
  
     dtold = dt
  
     ! if we're within 1%, return without changing dt.
     cflnew = dt*vm
  !   print *, "vmax:", vm,dt, cflnew
     if(abs(cfl-cflnew)/cfl < 0.01) then
  !      print *, "returning without changing dt"
      return
     endif
  
  !      print *, "setting dt from ",dt,"to",cfl/vm
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
  
  subroutine test_multi()
    use Geometry
  
    real(kind=r8) :: a,b
  !  real(kind=r8) :: rlambda
  !  complex(kind=r8),dimension(2) :: clambda
    integer :: i,j
  
    print *, "ndt: ",ndt
    a = 8.0_r8*PI
    b = 10.0_r8*PI
  
    do i=0,nx
       do j=0,ny
          slab(1)%psi(i,j) = sin(a*grid(1,i))*sin(b*grid(2,j))
          slab(2)%psi(i,j) = sin(b*grid(1,i))*sin(a*grid(2,j))
    !      slab(2)%psi(i,j) = 0.0_r8
       enddo
    enddo
    
    !call random_number(slab(2)%psi(1:4,3:5))
    
    call tcxcy(slab(1)%psi)
    call tcxcy(slab(2)%psi)
    
    
  
    return
  end subroutine test_multi
  subroutine layer_step_rg(t)
  !  integer,intent(in) :: n
    real(kind=r8),intent(inout) :: t
    integer :: i
  
    do i=1,nlayers
      slab(i)%use_q => slab(i)%q
      slab(i)%jac_result => slab(i)%k1
    enddo
    
    do i=1,nlayers
         call layer_jac(slab(i))
    end do
    
    slab(1)%k1 = slab(1)%k1 + forcing
    
    do i=1,nlayers
       slab(i)%k3 = slab(i)%q + dt3*slab(i)%k1
    end do
    call solve_p_nondim(slab(1)%k3,slab(1)%psi,slab(1)%s,alpha)
    
    
    do i=1,nlayers
      slab(i)%use_q => slab(i)%k3
      slab(i)%jac_result => slab(i)%k2
    enddo
    
    do i=1,nlayers
         call layer_jac(slab(i))
    end do
    
    slab(1)%k2 = slab(1)%k2+forcing
    do i=1,nlayers
       slab(i)%k2 = slab(i)%q + dt23*slab(i)%k2
    end do
    call solve_p_nondim(slab(1)%k2,slab(1)%psi,slab(1)%s,alpha)
    
    
    
    
    
    do i=1,nlayers
      slab(i)%use_q => slab(i)%k2
      slab(i)%jac_result => slab(i)%k3
    enddo
    
    do i=1,nlayers
         call layer_jac(slab(i))
    end do
    
    slab(1)%k3 = slab(1)%k3+forcing
    
    ! commenting out the last part of this line removes the explicit
    ! terms from the equation.  We need the rest of the computation
    ! to avoid division by zero problems.
    
    do i=1,nlayers
       slab(i)%q = slab(i)%q + dt4*(slab(i)%k1+3*slab(i)%k3)
    end do
    
    call solve_p_nondim(slab(1)%q,slab(1)%psi,slab(1)%s,alpha)
    
    
  
    do i=1,nlayers
       slab(i)%ddx = slab(i)%q
       slab(i)%ddy = slab(i)%q
       call d_dx_cheb(slab(i)%ddx)
       call d_dx_cheb(slab(i)%ddx)
       call d_dy_cheb(slab(i)%ddy)
       call d_dy_cheb(slab(i)%ddy)
       slab(i)%del2 = -(alpha**2*slab(i)%ddx+slab(i)%ddy+ndt*slab(i)%q)
    
       call solve_p_nondim(slab(i)%del2,slab(i)%q,ndt,alpha)
    enddo
    
  
    call solve_p_nondim(slab(1)%q,slab(1)%psi,slab(1)%s,alpha)
    
  
    ! update time, timestep
    t = t+dt
    call checkstable_slab()
    
  
    return
  end subroutine layer_step_rg
  subroutine layer_step_bt(t)
  !  integer,intent(in) :: n
    real(kind=r8),intent(inout) :: t
    integer :: i
  
    slab(1)%use_q => slab(1)%q
    slab(1)%jac_result => slab(1)%k1
    
    call layer_jac(slab(1))
    
    slab(1)%k1 = slab(1)%k1 + forcing
    slab(1)%k3 = slab(1)%q + dt3*slab(1)%k1
    
    call solve_p_nondim(slab(1)%k3,slab(1)%psi,0.0_r8,alpha)
    
    
    slab(1)%use_q => slab(1)%k3
    slab(1)%jac_result => slab(1)%k2
    
    call layer_jac(slab(1))
    
    slab(1)%k2 = slab(1)%k2+forcing
    slab(1)%k2 = slab(1)%q + dt23*slab(1)%k2
    
    call solve_p_nondim(slab(1)%k2,slab(1)%psi,0.0_r8,alpha)
    
    
    slab(1)%use_q => slab(1)%k2
    slab(1)%jac_result => slab(1)%k3
    
    call layer_jac(slab(1))
    
    slab(1)%k3 = slab(1)%k3+forcing
    
    ! commenting out the last part of this line removes the explicit
    ! terms from the equation.  We need the rest of the computation
    ! to avoid division by zero problems.
    
    slab(1)%q = slab(1)%q + dt4*(slab(1)%k1+3*slab(1)%k3)
    
    call solve_p_nondim(slab(1)%q,slab(1)%psi,0.0_r8,alpha)
    
    
  
    do i=1,nlayers
       slab(i)%ddx = slab(i)%q
       slab(i)%ddy = slab(i)%q
       call d_dx_cheb(slab(i)%ddx)
       call d_dx_cheb(slab(i)%ddx)
       call d_dy_cheb(slab(i)%ddy)
       call d_dy_cheb(slab(i)%ddy)
       slab(i)%del2 = -(alpha**2*slab(i)%ddx+slab(i)%ddy+ndt*slab(i)%q)
    
       call solve_p_nondim(slab(i)%del2,slab(i)%q,ndt,alpha)
    enddo
    
  
    call solve_p_nondim(slab(1)%q,slab(1)%psi,0.0_r8,alpha)
    
  
    ! update time, timestep
    t = t+dt
    call checkstable_slab()
    
  
    return
  end subroutine layer_step_bt
  subroutine layer_step_bf(t)
  !  integer,intent(in) :: n
    real(kind=r8),intent(inout) :: t
    integer :: i
  
    slab(1)%use_q => slab(1)%q
    slab(1)%jac_result => slab(1)%k1
    
    call layer_jac(slab(1))
    
    slab(1)%k1 = slab(1)%k1 + forcing
    slab(1)%k3 = slab(1)%q + dt3*slab(1)%k1
    
    call solve_p_nondim(slab(1)%k3,slab(1)%psi,0.0_r8,alpha)
    
    
    slab(1)%use_q => slab(1)%k3
    slab(1)%jac_result => slab(1)%k2
    
    call layer_jac(slab(1))
    
    slab(1)%k2 = slab(1)%k2+forcing
    slab(1)%k2 = slab(1)%q + dt23*slab(1)%k2
    
    call solve_p_nondim(slab(1)%k2,slab(1)%psi,0.0_r8,alpha)
    
    
    slab(1)%use_q => slab(1)%k2
    slab(1)%jac_result => slab(1)%k3
    
    call layer_jac(slab(1))
    
    slab(1)%k3 = slab(1)%k3+forcing
    
    ! commenting out the last part of this line removes the explicit
    ! terms from the equation.  We need the rest of the computation
    ! to avoid division by zero problems.
    
    slab(1)%q = slab(1)%q + dt4*(slab(1)%k1+3*slab(1)%k3)
    
    call solve_p_nondim(slab(1)%q,slab(1)%psi,0.0_r8,alpha)
    
    
  
    do i=1,nlayers
       slab(i)%ddx = slab(i)%q
       slab(i)%ddy = slab(i)%q
       call d_dx_cheb(slab(i)%ddx)
       call d_dx_cheb(slab(i)%ddx)
       call d_dy_cheb(slab(i)%ddy)
       call d_dy_cheb(slab(i)%ddy)
    
       bft = ndt-slab(1)%s
    
       slab(i)%del2 = -(alpha**2*slab(i)%ddx+slab(i)%ddy+bft*slab(i)%q)
    
       bft = ndt+slab(1)%s
    
       call solve_p_nondim(slab(i)%del2,slab(i)%q,bft,alpha)
    enddo
    
  
    call solve_p_nondim(slab(1)%q,slab(1)%psi,0.0_r8,alpha)
    
  
    ! update time, timestep
    t = t+dt
    call checkstable_slab()
    
  
    return
  end subroutine layer_step_bf
  subroutine init_mode(l,ain,bin,amp)
    type(layer),intent(inout) :: l
    real(kind=r8),intent(in) :: ain,bin,amp
    real(kind=r8) :: a,b,pi
    integer :: i,j
  
    pi = 4.0_r8*atan(1.0_r8)
    a = ain*pi
    b = bin*pi
  
    do i=0,nx
       do j=0,ny
          l%psi(i,j) = amp*sin(a*grid(1,i))*sin(b*grid(2,j))
          l%dx(i,j) = amp*a*cos(a*grid(1,i))*sin(b*grid(2,j))
          l%dy(i,j) = amp*b*sin(a*grid(1,i))*cos(b*grid(2,j))
  ! del**2 psi and its x and y derivatives
          l%q(i,j) = -(alpha**2*a**2+b**2)*l%psi(i,j)
          l%ddx(i,j) = -(alpha**2*a**2+b**2)*l%dx(i,j)
          l%ddy(i,j) = -(alpha**2*a**2+b**2)*l%dy(i,j)
       enddo
    enddo
  
    return
  end subroutine init_mode
  subroutine layer_step_rg_split(t)
  !  integer,intent(in) :: n
    real(kind=r8),intent(inout) :: t
    integer :: i
  
    do i=1,nlayers
      slab(i)%use_q => slab(i)%q
      slab(i)%jac_result => slab(i)%k1
    enddo
    
    do i=1,nlayers
         call layer_jac(slab(1))
    end do
    
    ! generate del2psi
    slab(1)%ddx = slab(1)%psi
    slab(1)%ddy = slab(1)%psi
    call d_dx_cheb(slab(1)%ddx)
    call d_dx_cheb(slab(1)%ddx)
    call d_dy_cheb(slab(1)%ddy)
    call d_dy_cheb(slab(1)%ddy)
    slab(1)%del2 = alpha**2*slab(1)%ddx+slab(1)%ddy
    
    slab(1)%k1 = slab(1)%k1+forcing - (a_2ds-a_4dm3*slab(1)%s)*slab(1)%del2
    
    do i=1,nlayers
       slab(i)%k3 = slab(i)%q + dt3*slab(i)%k1
    end do
    call solve_p_nondim(slab(1)%k3,slab(1)%psi,slab(1)%s,alpha)
    
    
    do i=1,nlayers
      slab(i)%use_q => slab(i)%k3
      slab(i)%jac_result => slab(i)%k2
    enddo
    
    do i=1,nlayers
         call layer_jac(slab(1))
    end do
    
    ! generate del2psi
    slab(1)%ddx = slab(1)%psi
    slab(1)%ddy = slab(1)%psi
    call d_dx_cheb(slab(1)%ddx)
    call d_dx_cheb(slab(1)%ddx)
    call d_dy_cheb(slab(1)%ddy)
    call d_dy_cheb(slab(1)%ddy)
    slab(1)%del2 = alpha**2*slab(1)%ddx+slab(1)%ddy
    
    slab(1)%k2 = slab(1)%k2+forcing - (a_2ds-a_4dm3*slab(1)%s)*slab(1)%del2
    
    do i=1,nlayers
       slab(i)%k2 = slab(i)%q + dt23*slab(i)%k2
    end do
    call solve_p_nondim(slab(1)%k2,slab(1)%psi,slab(1)%s,alpha)
    
    
    do i=1,nlayers
      slab(i)%use_q => slab(i)%k2
      slab(i)%jac_result => slab(i)%k3
    enddo
    
    do i=1,nlayers
         call layer_jac(slab(1))
    end do
    ! generate del2psi
    slab(1)%ddx = slab(1)%psi
    slab(1)%ddy = slab(1)%psi
    call d_dx_cheb(slab(1)%ddx)
    call d_dx_cheb(slab(1)%ddx)
    call d_dy_cheb(slab(1)%ddy)
    call d_dy_cheb(slab(1)%ddy)
    slab(1)%del2 = alpha**2*slab(1)%ddx+slab(1)%ddy
    
    
    slab(1)%k3 = slab(1)%k3+forcing - (a_2ds-a_4dm3*slab(1)%s)*slab(1)%del2
    
    ! commenting out the last part of this line removes the explicit
    ! terms from the equation.  We need the rest of the computation
    ! to avoid division by zero problems.
    
    do i=1,nlayers
       slab(i)%q = slab(i)%q + dt4*(slab(i)%k1+3*slab(i)%k3)
    end do
    
    call solve_p_nondim(slab(1)%q,slab(1)%psi,slab(1)%s,alpha)
    
    
  
    do i=1,nlayers
       slab(i)%ddx = slab(i)%q
       slab(i)%ddy = slab(i)%q
       call d_dx_cheb(slab(i)%ddx)
       call d_dx_cheb(slab(i)%ddx)
       call d_dy_cheb(slab(i)%ddy)
       call d_dy_cheb(slab(i)%ddy)
       slab(i)%del2 = -(alpha**2*slab(i)%ddx+slab(i)%ddy+ndt*slab(i)%q)
    
       call solve_p_nondim(slab(i)%del2,slab(i)%q,ndt,alpha)
    enddo
    
  
    call solve_p_nondim(slab(1)%q,slab(1)%psi,slab(1)%s,alpha)
    
  
    ! update time, timestep
    t = t+dt
    call checkstable_slab()
    
  
    return
  end subroutine layer_step_rg_split
  subroutine layer_step_bf_split(t)
  !  integer,intent(in) :: n
    real(kind=r8),intent(inout) :: t
    integer :: i
  
    slab(1)%use_q => slab(1)%q
    slab(1)%jac_result => slab(1)%k1
    
    call layer_jac(slab(1))
    
    slab(1)%k1 = slab(1)%k1 + forcing - (a_2ds-a_4dm3*slab(1)%s)*slab(1)%use_q
    slab(1)%k3 = slab(1)%q + dt3*slab(1)%k1
    
    call solve_p_nondim(slab(1)%k3,slab(1)%psi,0.0_r8,alpha)
    
    
    slab(1)%use_q => slab(1)%k3
    slab(1)%jac_result => slab(1)%k2
    
    call layer_jac(slab(1))
    
    slab(1)%k2 = slab(1)%k2+forcing - (a_2ds-a_4dm3*slab(1)%s)*slab(1)%use_q
    slab(1)%k2 = slab(1)%q + dt23*slab(1)%k2
    
    call solve_p_nondim(slab(1)%k2,slab(1)%psi,0.0_r8,alpha)
    
    
    slab(1)%use_q => slab(1)%k2
    slab(1)%jac_result => slab(1)%k3
    
    call layer_jac(slab(1))
    
    slab(1)%k3 = slab(1)%k3+forcing - (a_2ds-a_4dm3*slab(1)%s)*slab(1)%use_q
    
    ! commenting out the last part of this line removes the explicit
    ! terms from the equation.  We need the rest of the computation
    ! to avoid division by zero problems.
    
    slab(1)%q = slab(1)%q + dt4*(slab(1)%k1+3*slab(1)%k3)
    
    call solve_p_nondim(slab(1)%q,slab(1)%psi,0.0_r8,alpha)
    
    
  
    do i=1,nlayers
       slab(i)%ddx = slab(i)%q
       slab(i)%ddy = slab(i)%q
       call d_dx_cheb(slab(i)%ddx)
       call d_dx_cheb(slab(i)%ddx)
       call d_dy_cheb(slab(i)%ddy)
       call d_dy_cheb(slab(i)%ddy)
    
       bft = ndt-slab(1)%s
    
       slab(i)%del2 = -(alpha**2*slab(i)%ddx+slab(i)%ddy+bft*slab(i)%q)
    
       bft = ndt+slab(1)%s
    
       call solve_p_nondim(slab(i)%del2,slab(i)%q,bft,alpha)
    enddo
    
  
    call solve_p_nondim(slab(1)%q,slab(1)%psi,0.0_r8,alpha)
    
  
    ! update time, timestep
    t = t+dt
    call checkstable_slab()
    
  
    return
  end subroutine layer_step_bf_split
  
end module Layers
