
program test_trans
  use Types
  use Geometry
  use NewTrans
  use OpenDX

  integer :: i,j

  call init_geometry(nnx=128,nny=128,chebx=.true.,cheby=.true.)
  call init_trans(nx,ny)
  dx_base = "test"
  call dx_write_positions_2d(grid(1,:),grid(2,:),nx,ny)

  print *, nx,ny

!  print *,maxval(grid(1,:)),maxval(grid(2,:))

  psi = 0.0_r8

  do i=0,nx
  do j=0,ny
     psi(i,j) = (1-grid(1,i)**6)*(1-grid(2,j)**6)
     w2(i,j) = -30.0_r8*(grid(1,i)**4)*(1-grid(2,j)**6)
     w3(i,j) = -30.0_r8*(1-grid(1,i)**6)*(grid(2,j)**4)
     w4(i,j) = w3(i,j) + w2(i,j)
  enddo
  enddo

!  psi(6,0) = 1.0_r8
!  call fcxcy(psi)
!  print *,"max difference T6:",maxval(abs(psi-w3))

!  print *,maxval(psi),maxloc(psi)
  call dx_write_data_2d(psi(0:nx,0:ny),0)
  call dx_write_header_2d(nx,ny,0.0,0)

  call dx_write_data_2d(w2(0:nx,0:ny),1)
  call dx_write_header_2d(nx,ny,0.0,1)
  call dx_write_data_2d(w3(0:nx,0:ny),2)
  call dx_write_header_2d(nx,ny,0.0,2)
  call dx_write_data_2d(w4(0:nx,0:ny),3)
  call dx_write_header_2d(nx,ny,0.0,3)

  call tcxcy(psi)
  u = psi
  call d_dx_cheb(u)
  call d_dx_cheb(u)
  v = psi
  call d_dy_cheb(v)
  call d_dy_cheb(v)
  psi = u+v
  call fcxcy(psi)
  call fcxcy(u)
  call fcxcy(v)

!  print *,maxval(psi),maxloc(psi)
!  print *,maxval(w2),maxloc(w2)


  w1 = w4-psi
  print *,"max difference:",maxval(abs(w1)),maxloc(abs(w1))


  call dx_write_data_2d(u(0:nx,0:ny),4)
  call dx_write_header_2d(nx,ny,0.0,4)

  call dx_write_data_2d(v(0:nx,0:ny),5)
  call dx_write_header_2d(nx,ny,0.0,5)

  call dx_write_data_2d(psi(0:nx,0:ny),6)
  call dx_write_header_2d(nx,ny,0.0,6)

  call dx_write_data_2d(w1(0:nx,0:ny),7)
  call dx_write_header_2d(nx,ny,0.0,7)

  stop
end program test_trans
