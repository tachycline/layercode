
program eigtest
  use Types
  use Geometry
  use NewTrans
  use OpenDX

! for solving the problem 5.1.40 in Canuto
!  real(kind=r8),dimension(:,:),allocatable,save :: u,uprime,f,fprime
  real(kind=r8),dimension(:,:),allocatable,save :: uprime,f,fprime
  real(kind=r8),dimension(:,:),allocatable,save :: uxx,uyy,forig,uorig
  real(kind=r8) :: lambda,xpt,ypt

! matricies for the similarity transform diagonalizing a,b
  real(kind=r8),dimension(:,:),allocatable,save :: q,qinv,p,pinv

! eigenvalues of a, b
  real(kind=r8),dimension(:),allocatable,save:: wra,wia,wrb,wib

! original derivative matricies (and smaller forms)
  real(kind=r8),dimension(:,:),allocatable,save :: b,a,bs,as

! indicies and limits
  integer :: i,j

! eigenvectors (for lapack)
  real(kind=r8),dimension(:,:),allocatable,save :: vr,vl

! work array (for lapack)
  real(kind=r8),dimension(:),allocatable,save :: work

! parameters for lapack routines
  integer :: lwork,info,ldh,ldvl,ldvr

! flags for lapack -- do we need left/right eigenvectors?
  character(len=1) :: jobvl,jobvr

! record of pivots in decomposition
  integer,dimension(:),allocatable,save :: ipiv

! small enough to look at in a terminal window.
!  nx = 8
!  ny = 8

  call init_geometry(nnx=64,nny=64,chebx = .true.,cheby=.true.)

! y derivative matrix -----------------------------------------
  allocate(b(0:ny,0:ny),bs(ny-1,ny-1))
  allocate(work(ny+1),wrb(ny-1),wib(ny-1))
  allocate(qinv(ny-1,ny-1),q(ny-1,ny-1))
  allocate(vr(ny-1,ny-1),vl(ny-1,ny-1))
  allocate(ipiv(ny-1))


  ! build a (d**2/dx**2 matrix)
  b = 0.0_r8
  do j = 2,ny,2
     b(0,j) = j*j*j/2.0_r8
  enddo
  do i=1,ny-1
     do j = i+2,ny,2
        b(i,j) = j*(j*j-i*i)
     enddo
  enddo

  !eliminate highest order coefficients with bc
  do i=0,ny-2
     b(i,0:ny-2:2) = b(i,0:ny-2:2) - b(i,ny)
     b(i,1:ny-3:2) = b(i,1:ny-3:2) - b(i,ny-1)
  enddo

  bs = b(0:ny-2,0:ny-2)
  bs = transpose(bs)
  b = transpose(b)

!! put the result in the smaller matrix
!  bs = b(0:ny-2,0:ny-2)
!!  bs = b(2:ny,2:ny)

! get eigenvalues,eigenvectors
  ldh = ny-1
  ldvl = ny-1
  ldvr = ny-1

! we want right eigenvectors, don't need left eigenvectors.
  jobvl = "N"
  jobvr = "V"
! get optimal size for work, set lwork=-1
  lwork = -1
  call dgeev(jobvl,jobvr,ny-1,bs,ny-1,wrb,wib,vl,ldvl,vr,ldvr,work,lwork,info)
! check for error condition here?
!  print *, "info: ",info
  lwork = work(1)
  deallocate(work)
  allocate(work(lwork))
! actual call to get eigenvalues, eigenvectors
  call dgeev(jobvl,jobvr,ny-1,bs,ny-1,wrb,wib,vl,ldvl,vr,ldvr,work,lwork,info)
! check for error condition here?
!  print *, "info: ",info
!print *,"max imaj eigenvalue:",maxval(abs(wib))
!print *,"max real eigenvalue:",maxval(abs(wrb))

!   print *,"eigenvalues of b:"
!   do i=1,nx-1
!      print "(2es15.6)",wrb(i),wib(i)
!   enddo

  ! now invert matrix of eigenvectors

  q = vr(1:ny-1,1:ny-1)
  qinv = q

  call dgetrf(ny-1,ny-1,qinv,ny-1,ipiv,info)
  lwork = -1
  call dgetri(ny-1,qinv,ny-1,ipiv,work,lwork,info)
  lwork = work(1)
  deallocate(work)
  allocate(work(lwork))
  call dgetri(ny-1,qinv,ny-1,ipiv,work,lwork,info)

  deallocate(ipiv,vr,vl,work)

! for testing purposes
  bs = b(0:ny-2,0:ny-2)
  bs = matmul(qinv,matmul(bs,q))
! output stuff goes here (if you're not looking at the fw file)
  bs = b(0:ny-2,0:ny-2)

! x derivative matrix --------------------------------------------
  allocate(a(0:nx,0:nx),as(nx-1,nx-1))
  allocate(work(nx+1),wra(nx-1),wia(nx-1))
  allocate(pinv(nx-1,nx-1),p(nx-1,nx-1))
  allocate(vr(nx-1,nx-1),vl(nx-1,nx-1))
  allocate(ipiv(nx-1))

  ! build a (d**2/dx**2 matrix)
  a = 0.0_r8
  do j = 2,nx,2
     a(0,j) = j*j*j/2.0_r8
  enddo
  do i=1,nx-1
     do j = i+2,nx,2
        a(i,j) = j*(j*j-i*i)
     enddo
  enddo


  !eliminate highest order coefficients with bc
  do i=0,nx-2
     a(i,0:nx-2:2) = a(i,0:nx-2:2) - a(i,nx)
     a(i,1:nx-3:2) = a(i,1:nx-3:2) - a(i,nx-1)
  enddo
! insert bc into last two rows
!  a(nx-1:nx,:) = 1.0_r8
!  a(nx-1,1:nx-1:2) = -1.0_r8

! put the result in the smaller matrix
  as = a(0:nx-2,0:nx-2)
!  as = a(2:nx,2:nx)

! get eigenvalues,eigenvectors
  ldh = nx-1
  ldvl = nx-1
  ldvr = nx-1

! we want right eigenvectors, don't need left eigenvectors.
  jobvl = "N"
  jobvr = "V"
! get optimal size for work, set lwork=-1
  lwork = -1
  call dgeev(jobvl,jobvr,nx-1,as,nx-1,wra,wia,vl,ldvl,vr,ldvr,work,lwork,info)
!  print *, "info: ",info
  lwork = work(1)
  deallocate(work)
  allocate(work(lwork))
! actual call to get eigenvalues, eigenvectors
  call dgeev(jobvl,jobvr,nx-1,as,nx-1,wra,wia,vl,ldvl,vr,ldvr,work,lwork,info)

!  print *, "info: ",info
!   print *,"eigenvalues of a:"
!   do i=1,nx-1
!      print "(2es15.6)",wra(i),wia(i)
!   enddo
!   print *,"maximum real eigenvalue: ",maxval(abs(wra)),maxloc(abs(wra))
!   print *,"maximum imaj eigenvalue: ",maxval(abs(wia)),maxloc(abs(wia))
!   print *,"maximum imaj eigenvalue: ",maxval(abs(wib)),maxloc(abs(wib))
  ! now invert matrix of eigenvectors
  p = vr(1:ny-1,1:ny-1)
  pinv = p

  call dgetrf(nx-1,nx-1,pinv,nx-1,ipiv,info)
  lwork = -1
  call dgetri(nx-1,pinv,nx-1,ipiv,work,lwork,info)
  lwork = work(1)
  deallocate(work)
  allocate(work(lwork))
  call dgetri(nx-1,pinv,nx-1,ipiv,work,lwork,info)

! for testing purposes
  as = a(0:nx-2,0:nx-2)

! check P**(-1)AP diagonal:
  as = matmul(pinv,matmul(as,p))
! output stuff goes here (if you're not looking at the fw file)
  as = a(0:nx-2,0:nx-2)


! end x derivative matrix ----------------------------------------------

  ! now, get on with it.
  ! allocate storage for problem
  allocate(u(0:nx,0:ny),f(0:nx,0:ny),forig(0:nx,0:ny))
  allocate(uxx(0:nx,0:ny),uyy(0:nx,0:ny),uorig(0:nx,0:ny))
  allocate(uprime(0:nx-2,0:ny-2),fprime(0:nx-2,0:ny-2))
  

  ! initialize transforms
  call init_trans(nx,ny)

  ! set parameters
  lambda = 0.0e-5_r8

  ! generate U, F in real space
  PI = 4.0_r8*atan(1.0_r8)
  do i=0,nx
     xpt = -cos(i*PI/nx)
     do j=0,ny
        ypt = -cos(j*PI/ny)
        u(i,j) = (xpt**6 - 15*xpt**2 + 14.0_r8) &
                 *(ypt**6 - 15*ypt**2 + 14.0_r8)
        f(i,j) = 30*(xpt**4-1.0_r8)*(ypt**6 - 15*ypt**2 + 14.0_r8) &
               + 30*(ypt**4-1.0_r8)*(xpt**6 - 15*xpt**2 + 14.0_r8) &
               - lambda * u(i,j)

     enddo
  enddo

  dx_base = "eigen"
  call dx_write_positions_2d(grid(1,:),grid(2,:),nx,ny)
  call dx_write_data_2d(u(0:nx,0:ny),0)
  call dx_write_header_2d(nx,ny,0.0,0)
  call dx_write_data_2d(f(0:nx,0:ny),1)
  call dx_write_header_2d(nx,ny,0.0,1)


  ! go to transform space
  call tcxcy(u)
  call tcxcy(f)
  forig = f
  uorig = u

! check matricies against derivatives
  uxx = u
  call d_dx_cheb(uxx)
  call d_dx_cheb(uxx)
  uyy = u
  call d_dy_cheb(uyy)
  call d_dy_cheb(uyy)

  uprime = uorig(0:nx-2,0:ny-2)
  uprime = matmul(as,uprime)
  u = 0.0_r8
  u(0:nx-2,0:ny-2) = uprime
    do i=0,ny-2
       u(nx-1,i) = -sum(u(1:nx-3:2,i))
       u(nx,i) = -sum(u(0:nx-2:2,i))
    enddo
    do i=0,nx-2
       u(i,ny-1) = -sum(u(i,1:ny-3:2))
       u(i,ny) = -sum(u(i,0:ny-2:2))
    enddo
    u(nx-1,ny-1) = -sum(u(1:nx-3:2,ny-1))
    u(nx-1,ny) = -sum(u(1:nx-3:2,ny))
    u(nx,ny-1) = -sum(u(0:nx-2:2,ny-1))
    u(nx,ny) = -sum(u(0:nx-2:2,ny))
  
  where (abs(u) <=  1.0e-10)
     u = 0.0_r8
  endwhere
  

  print *,"l2 norm, u,uxx:",sqrt(sum((u-uxx)**2))
  
  print *,"linf norm, u,uxx:",maxval(abs(u-uxx)),&
          "at",maxloc(abs(u-uxx))
  

  uprime = uorig(0:nx-2,0:ny-2)
  uprime = matmul(uprime,bs)
  u = 0.0_r8
  u(0:nx-2,0:ny-2) = uprime
    do i=0,ny-2
       u(nx-1,i) = -sum(u(1:nx-3:2,i))
       u(nx,i) = -sum(u(0:nx-2:2,i))
    enddo
    do i=0,nx-2
       u(i,ny-1) = -sum(u(i,1:ny-3:2))
       u(i,ny) = -sum(u(i,0:ny-2:2))
    enddo
    u(nx-1,ny-1) = -sum(u(1:nx-3:2,ny-1))
    u(nx-1,ny) = -sum(u(1:nx-3:2,ny))
    u(nx,ny-1) = -sum(u(0:nx-2:2,ny-1))
    u(nx,ny) = -sum(u(0:nx-2:2,ny))
  
  where (abs(u) <=  1.0e-10)
     u = 0.0_r8
  endwhere
  

  print *,"l2 norm, u,uyy:",sqrt(sum((u-uyy)**2))
  
  print *,"linf norm, u,uyy:",maxval(abs(u-uyy)),&
          "at",maxloc(abs(u-uyy))
  

! check diagonalized derivatives
! i.e., PinvAPPinvUQ = PinvUxxQ
!       PinvUQQinvBQ = pinvUyyQ
print *,"------------- diagonalized derivatives -----------"

! first x
  fprime = uxx(0:nx-2,0:ny-2)
  fprime = matmul(pinv,matmul(fprime,q))
  uprime = uorig(0:nx-2,0:ny-2)
! first, get PinvUQ
  uprime = matmul(pinv,matmul(uprime,q))
! next, get PinvAPU'
  uprime = matmul(pinv,matmul(as,matmul(p,uprime)))

  print *,"l2 norm, uprime,fprime:",sqrt(sum((uprime-fprime)**2))
  
  print *,"linf norm, uprime,fprime:",maxval(abs(uprime-fprime)),&
          "at",maxloc(abs(uprime-fprime))
  

! now y
  fprime = uyy(0:nx-2,0:ny-2)
  fprime = matmul(pinv,matmul(fprime,q))
  uprime = uorig(0:nx-2,0:ny-2)
! first, get PinvUQ
  uprime = matmul(pinv,matmul(uprime,q))
! next, get U'QinvBQ
  uprime = matmul(uprime,matmul(qinv,matmul(bs,q)))

  print *,"l2 norm, uprime,fprime:",sqrt(sum((uprime-fprime)**2))
  
  print *,"linf norm, uprime,fprime:",maxval(abs(uprime-fprime)),&
          "at",maxloc(abs(uprime-fprime))
  

print *,"------------- diagonalized eigenvalues -----------"
! check eigenvalues vs. diagonalized as, bs
  as = matmul(pinv,matmul(as,p))
  do i=1,nx-1
     wia(i) = as(i,i)
  enddo
  print *,"l2 norm, wia,wra:",sqrt(sum((wia-wra)**2))
  
  print *,"linf norm, wia,wra:",maxval(abs(wia-wra)),&
          "at",maxloc(abs(wia-wra))
  

  bs = matmul(qinv,matmul(bs,q))
  do i=1,nx-1
     wib(i) = bs(i,i)
  enddo
  print *,"l2 norm, wib,wrb:",sqrt(sum((wib-wrb)**2))
  
  print *,"linf norm, wib,wrb:",maxval(abs(wib-wrb)),&
          "at",maxloc(abs(wib-wrb))
  

print *,"------------- diagonalized equation -----------"
! check the diagonalized equation
!  as = a(0:nx-2,0:nx-2)
!  as = matmul(pinv,matmul(as,p))

!  bs = b(0:nx-2,0:nx-2)
!  bs = matmul(qinv,matmul(bs,q))

  uprime = uorig(0:nx-2,0:ny-2)
  uprime = matmul(pinv,matmul(uprime,q))

  uprime = matmul(as,uprime) + matmul(uprime,bs) - lambda*uprime

  fprime = forig(0:nx-2,0:ny-2)
  fprime = matmul(pinv,matmul(fprime,q))

  print *,"l2 norm, uprime,fprime:",sqrt(sum((uprime-fprime)**2))
  
  print *,"linf norm, uprime,fprime:",maxval(abs(uprime-fprime)),&
          "at",maxloc(abs(uprime-fprime))
  

print *, "----------- invert --------------"
  uprime = matmul(p,matmul(uprime,qinv))
  uprime = matmul(pinv,matmul(uprime,q))

  print *,"l2 norm, uprime,fprime:",sqrt(sum((uprime-fprime)**2))
  
  print *,"linf norm, uprime,fprime:",maxval(abs(uprime-fprime)),&
          "at",maxloc(abs(uprime-fprime))
  

print *, "------------ forig --------------"
  fprime = matmul(p,matmul(uprime,qinv))
  f(0:nx-2,0:ny-2) = fprime
    do i=0,ny-2
       f(nx-1,i) = -sum(f(1:nx-3:2,i))
       f(nx,i) = -sum(f(0:nx-2:2,i))
    enddo
    do i=0,nx-2
       f(i,ny-1) = -sum(f(i,1:ny-3:2))
       f(i,ny) = -sum(f(i,0:ny-2:2))
    enddo
    f(nx-1,ny-1) = -sum(f(1:nx-3:2,ny-1))
    f(nx-1,ny) = -sum(f(1:nx-3:2,ny))
    f(nx,ny-1) = -sum(f(0:nx-2:2,ny-1))
    f(nx,ny) = -sum(f(0:nx-2:2,ny))
  
  where (abs(f) <= 1.0e-10)
     f = 0.0_r8
  endwhere
  
  print *,"l2 norm, forig,f:",sqrt(sum((forig-f)**2))
  
  print *,"linf norm, forig,f:",maxval(abs(forig-f)),&
          "at",maxloc(abs(forig-f))
  

print *, "---------- equation with eigenvalues ------------"
  uprime = uorig(0:nx-2,0:ny-2)
  uprime = matmul(pinv,matmul(uprime,q))

  do i=0,nx-2
     do j=0,ny-2
        fprime(i,j) = (wra(i+1) + wrb(j+1) - lambda)*uprime(i,j)
     enddo
  enddo

  fprime = matmul(p,matmul(fprime,qinv))
  f(0:nx-2,0:ny-2) = fprime
    do i=0,ny-2
       f(nx-1,i) = -sum(f(1:nx-3:2,i))
       f(nx,i) = -sum(f(0:nx-2:2,i))
    enddo
    do i=0,nx-2
       f(i,ny-1) = -sum(f(i,1:ny-3:2))
       f(i,ny) = -sum(f(i,0:ny-2:2))
    enddo
    f(nx-1,ny-1) = -sum(f(1:nx-3:2,ny-1))
    f(nx-1,ny) = -sum(f(1:nx-3:2,ny))
    f(nx,ny-1) = -sum(f(0:nx-2:2,ny-1))
    f(nx,ny) = -sum(f(0:nx-2:2,ny))
  
  where (abs(f) <= 1.0e-10)
     f = 0.0_r8
  endwhere
  
  print *,"l2 norm, forig,f:",sqrt(sum((forig-f)**2))
  
  print *,"linf norm, forig,f:",maxval(abs(forig-f)),&
          "at",maxloc(abs(forig-f))
  

print *, "------------- solving ---------------"
  fprime = forig(0:nx-2,0:ny-2)
  fprime = matmul(pinv,matmul(fprime,q))
  do i=0,nx-2
     do j=0,ny-2
        uprime(i,j) = fprime(i,j)/(wra(i+1) + wrb(j+1) - lambda)
     enddo
  enddo
  uprime = matmul(p,matmul(uprime,qinv))
  u(0:nx-2,0:ny-2) = uprime
  where (abs(u) <= 1.0e-10)
     u = 0.0_r8
  endwhere
  

  print *,"l2 norm, uorig,u:",sqrt(sum((uorig-u)**2))
  
  print *,"linf norm, uorig,u:",maxval(abs(uorig-u)),&
          "at",maxloc(abs(uorig-u))
  

print *, "------------- solving (transmod) ---------------"
  f = forig
  call solve_p(f,u,lambda)

  print *,"l2 norm, uorig,u:",sqrt(sum((uorig-u)**2))
  
  print *,"linf norm, uorig,u:",maxval(abs(uorig-u)),&
          "at",maxloc(abs(uorig-u))
  

  stop
end program eigtest
