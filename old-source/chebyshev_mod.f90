module Chebyshev
  use Types

  implicit none
  include "fftw3.f"

  real(kind=r8),dimension(:,:),allocatable,private :: tkj
  real(kind=r8),dimension(:,:),allocatable,private :: tkjhat
  integer(kind=i4),private :: forward_plan,inverse_plan
  
  public :: grid_gauss_cheb,grid_gauss_lobatto,cheb_n
  public :: init_matmul_inv,fromcheb_matmul
  public :: init_matmul_for,tocheb_matmul
  

contains

  subroutine grid_gauss_cheb(x)
    real(kind=r8),dimension(:),intent(inout) :: x
    integer :: nx,j
    real(kind=r8) :: pi
  
    pi = 4.0*atan(1.0d0)
  
    nx = size(x)
    do j=0,nx-1
       x(j+1) = cos(pi*(j+0.5)/real(nx))
    enddo
    return
  
  end subroutine grid_gauss_cheb
  
  subroutine grid_gauss_lobatto(x)
    real(kind=r8),dimension(:),intent(inout) :: x
    integer :: nx,j
    real(kind=r8) :: pi
  
    pi = 4.0*atan(1.0d0)
  
    nx = size(x)
    do j=0,nx-1
       x(j+1) = cos(pi*(j)/real(nx))
    enddo
    return
  
  end subroutine grid_gauss_lobatto
  
  subroutine cheb_n(t,x,n)
    real(kind=r8),dimension(:),intent(inout) :: t
    real(kind=r8),dimension(:),intent(in) :: x
    integer,intent(in) :: n
    integer :: nx,i
  
    nx = size(t)
    do i=1,nx
       t(i) = cos(n*acos(x(i)))
    enddo
    return
  end subroutine cheb_n
  
  subroutine init_matmul_inv(nx)
    integer,intent(in) :: nx
    real(kind=r8),dimension(:),allocatable :: t,x
    integer :: j
  
    allocate(tkj(nx,nx),t(nx),x(nx))
    call grid_gauss_cheb(x)
    do j=0,nx-1
      call cheb_n(t,x,j)
      tkj(:,j+1) = t(:)
    enddo
  
    deallocate(t,x)
    return
  end subroutine init_matmul_inv
  
  subroutine fromcheb_matmul(psihat,psi)
    real(kind=r8),dimension(:),intent(inout) :: psihat,psi
  
    psi = matmul(tkj,psihat)
  
    return
  end subroutine fromcheb_matmul
  
  subroutine init_matmul_for(nx)
    integer,intent(in) :: nx
    real(kind=r8),dimension(:),allocatable :: t,x
    integer :: j
    real(kind=r8) :: pi
  
    pi = 4.0_r8*atan(1.0_r8)
  
    allocate(tkjhat(nx,nx),t(nx),x(nx))
    call grid_gauss_cheb(x)
    do j=0,nx-1
      call cheb_n(t,x,j)
      tkjhat(j+1,:) = t(:)*2.0/(nx)
    enddo
    tkjhat(1,:) = tkjhat(1,:)/2.0
  
    deallocate(t,x)
    return
  end subroutine init_matmul_for
  
  subroutine tocheb_matmul(psihat,psi)
    real(kind=r8),dimension(:),intent(inout) :: psihat,psi
  
    psi = matmul(tkjhat,psihat)
  
    return
  end subroutine tocheb_matmul
  
  

end module Chebyshev
