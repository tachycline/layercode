
program test_bvp
  use Types
  use Geometry
  use NewTrans
  use OpenDX

  real(kind=r8),dimension(:,:),allocatable,save :: u,uprime,f,fprime
  real(kind=r8),dimension(:,:),allocatable,save :: uxx,uyy,forig,uorig
  real(kind=r8) :: lambda,xpt,ypt
  integer :: i,j

  call init_geometry(nnx=64,nny=64,chebx=.true.,cheby=.true.)
  call init_trans(nx,ny)
  call init_plots(basename="bvp",xvec=grid(1,:),&
                  yvec=grid(2,:),nx=nx,ny=ny,nrun=0)

  allocate(u(0:nx,0:ny),f(0:nx,0:ny),forig(0:nx,0:ny))
  allocate(uxx(0:nx,0:ny),uyy(0:nx,0:ny),uorig(0:nx,0:ny))
  allocate(uprime(0:nx-2,0:ny-2),fprime(0:nx-2,0:ny-2))
  

  lambda = 0.0_r8


  ! generate U, F in real space
  PI = 4.0_r8*atan(1.0_r8)
  do i=0,nx
     xpt = -cos(i*PI/nx)
     do j=0,ny
        ypt = -cos(j*PI/ny)
        u(i,j) = (xpt**6 - 15*xpt**2 + 14.0_r8) &
                 *(ypt**6 - 15*ypt**2 + 14.0_r8)
        f(i,j) = 120*(xpt**4-1.0_r8)*(ypt**6 - 15*ypt**2 + 14.0_r8) &
               + 30*(ypt**4-1.0_r8)*(xpt**6 - 15*xpt**2 + 14.0_r8) &
               - lambda * u(i,j)

     enddo
  enddo

  call plot_field(fieldname="forig",field=f,n=0,t=0.0_r8)


  call tcxcy(u)
  call tcxcy(f)

  uorig = u
  forig = f
  call plot_field(fieldname="uorig-cheb",field=u,n=0,t=0.0_r8)
  call plot_field(fieldname="forig-cheb",field=f,n=0,t=0.0_r8)

  call solve_p_nondim(f,u,lambda,2.0_r8)

  call plot_field(fieldname="u-cheb",field=u,n=0,t=0.0_r8)

  print *, "cheb space norms"
  print *,"linf norm, u,uorig:",maxval(abs(u-uorig)),&
          "at",maxloc(abs(u-uorig))
  
  print *,"linf norm, u,0:",maxval(abs(u-0)),&
          "at",maxloc(abs(u-0))
  
  print *,"linf norm, uorig,0:",maxval(abs(uorig-0)),&
          "at",maxloc(abs(uorig-0))
  

  uxx = u-uorig
  print *,"linf norm, uxx,0:",maxval(abs(uxx-0)),&
          "at",maxloc(abs(uxx-0))
  
  call plot_field(fieldname="udiff-cheb",field=u,n=0,t=0.0_r8)


  call fcxcy(u)
  call fcxcy(uorig)

  call plot_field(fieldname="uorig",field=uorig,n=0,t=0.0_r8)
  call plot_field(fieldname="u",field=u,n=0,t=0.0_r8)

  print *, "real space norms"
  print *,"linf norm, u,uorig:",maxval(abs(u-uorig)),&
          "at",maxloc(abs(u-uorig))
  
  print *,"l2 norm, u,uorig:",sqrt(sum((u-uorig)**2))
  

  u = u/2.5_r8-uorig
  call plot_field(fieldname="udiff",field=u,n=0,t=0.0_r8)

  stop
end program test_bvp
