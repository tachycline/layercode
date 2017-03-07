program chebtest
  use Types
  use Chebyshev

  real(kind=r8),dimension(:),allocatable :: x,t
  integer :: n,nx,j

  nx=128
  allocate(x(nx),t(nx))

  call plinit()
  call plenv(-1.0_r8,1.0_r8,-1.0_r8,1.0_r8,0,0)

  do n=0,5
     call grid_gauss_cheb(x)
     call cheb_n(t,x,n)
     call plcol(2)
     call plline(nx,x,t)

     call grid_gauss_lobatto(x)
     call cheb_n(t,x,n)
     call plcol(3)
     call plline(nx,x,t)
     call plcol(1)
  enddo
  deallocate(x,t)
  call plend()

  stop
end program chebtest
