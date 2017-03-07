
program test_trans
  use Types
  use Geometry
  use NewTrans
  use OpenDX

  integer :: i,j

  call init_geometry(nnx=128,nny=64,chebx=.true.,cheby=.true.)
  call init_trans(nx,ny)

  dx_base = "test"
  call dx_write_positions_2d(grid(1,:),grid(2,:),nx,ny)

  ! set up some initial data
  psi = 0.0_r8

!  psi(10,2) = 1.0_r8
!  call fcxcy(psi)

  do i=0,nx
    do j=0,ny
      psi(i,j) = cos(PI*grid(1,i))*cos(3*PI*grid(2,j))
    enddo
  enddo

  ! take two derivatives
  call tcxcy(psi)
  w1 = psi

  call d_dy_cheb(psi)
  call d_dy_cheb(psi)

  ! do matrix multiply

  w2 = matmul(w1,b)

  ! compare results
  call fcxcy(psi)
  call fcxcy(w1)
  call fcxcy(w2)

  w1 = w2 - psi
  print *, "max difference: ",maxval(w1),"at ",maxloc(w1)

  call dx_write_data_2d(psi,0)
  call dx_write_header_2d(nx,ny,0.0,0)
  call dx_write_data_2d(w2,1)
  call dx_write_header_2d(nx,ny,0.0,1)
  call dx_write_data_2d(w1,2)
  call dx_write_header_2d(nx,ny,0.0,2)

  stop
end program test_trans
