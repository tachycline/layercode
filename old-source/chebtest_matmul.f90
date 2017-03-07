program chebtest
  use Types
  use Chebyshev

  real(kind=r8),dimension(:),allocatable :: target,psi,psihat,x,error,orig
  integer :: nx,k,j
  real(kind=r8) :: l2,linf

  nx = 100
  allocate(target(nx),psi(nx),psihat(nx),x(nx),error(nx),orig(nx))

  call init_matmul_inv(nx)
  call init_matmul_for(nx)
  call grid_gauss_cheb(x)
! individual modes

  call plinit()

  do k=0,10
! inverse
     psi = 0.0_r8
     psihat = 0.0_r8
     call cheb_n(target,x,k)
     psihat(k+1) = 1.0_r8
     orig = psihat
     call fromcheb_matmul(psihat,psi)
     error = psi - target
     l2 = sum(error**2)
     print *, "mode",k,"l2:",l2
     linf = maxval(error)
     print *, "mode",k,"linf:",linf
     call plenv(-1.0_r8,1.0_r8,-1.0_r8,1.0_r8,0,0)
     call pllab("","","Inverse Transform")
     call plline(nx,x,error)
     call plcol(2)
     call plline(nx,x,psi)
     call plcol(3)
     call plline(nx,x,target)
     call plcol(1)

! forward
     call tocheb_matmul(psi,psihat)
     error = psihat - orig
     l2 = sum(error**2)
     print *, "mode",k,"l2:",l2
     linf = maxval(error)
     print *, "mode",k,"linf:",linf
     call plenv(-1.0_r8,1.0_r8,-1.0_r8,maxval(error),0,0)
     call pllab("","","Forward Transform")
     call plline(nx,x,error)
     call plcol(2)
     call plline(nx,x,psihat)
     call plcol(3)
!     call plline(nx,x,orig)
     call plcol(1)
  enddo
  call plend()

  stop
end program chebtest
