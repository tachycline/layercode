module Geometry
  use Types ! sets r8 to be the right kind for double precision
  private

  public :: init_geometry
  private :: init_grid
  ! number of dimensions
  integer,save,public :: ndim

  ! shape of the fields
  integer,dimension(:),allocatable,save,public :: sizes

  ! shorthands (useful for transforms)
  integer, public :: nx,nx2,ny,ny2,nz,nz2

  ! the axes
  real(kind=r8),dimension(:,:),allocatable,save,public :: grid,deltas
  real(kind=r8),public :: PI
!  real(kind=r8),public :: dx
contains
  subroutine init_geometry(nnx,chebx,nny,cheby,nnz,chebz)
    integer,intent(in) :: nnx
    logical,optional,intent(in) :: chebx
    integer,optional,intent(in) :: nny,nnz
    logical,optional,intent(in) :: cheby,chebz
    !------------------------- end interface
    logical,dimension(:),allocatable :: chebs

    ! set up PI
    PI = 4.0_r8*atan(1.0_r8)

    ! figure out ndim based on the arguments.
    if (present(nnz)) then
       ndim = 3
    else if (present(nny)) then
       ndim = 2
    else
       ndim = 1
    end if

    ! set up sizes, grid
    ! if these have already been allcoated, deallocate them.

    if(allocated(sizes)) then
       deallocate(sizes)
    end if

    if(allocated(chebs)) then
       deallocate(chebs)
    end if

    allocate(sizes(ndim),chebs(ndim))

    if (ndim > 2) then ! do z
       sizes(3) = nnz
       nz = nnz
       nz2 = nz/2
       if (present(chebz)) then
          chebs(3) = chebz
       else
          chebs(3) = .false.
       end if
    end if

    if (ndim > 1) then ! do y
       sizes(2) = nny
       ny = nny
       ny2 = ny/2
       if (present(cheby)) then
          chebs(2) = cheby
       else
          chebs(2) = .false.
       end if
    end if

    ! ndim is always >0, so do x
    sizes(1) = nnx
    nx = nnx
    nx2 = nx/2
    if(present(chebx)) then
       chebs(1) = chebx
    else
       chebs(1) = .false.
    end if

    call init_grid(chebs)

    return
  end subroutine init_geometry

  subroutine init_grid(chebs)
    logical,dimension(:),intent(in) :: chebs
    ! the other necessaries are module data
    integer :: i,j
    ! allocate the grid -- make sure it's big enough
    ! if a grid has already been allocated, deallocate it.
    if(allocated(grid)) then
      deallocate(grid)
      deallocate(deltas)
    end if

    allocate(grid(1:ndim,0:maxval(sizes)))
    allocate(deltas(1:ndim,0:maxval(sizes)))
    ! set up gauss-lobotto points for cheb dimensions, regular grid for fourier
    do i=1,ndim
       if (chebs(i)) then ! chebyshev in this dimension, use gauss-lobotto
          do j=0,sizes(i)
             grid(i,j) = -cos(PI*j/sizes(i))
          end do
       else  ! fourier in this dimension, use regular grid spacing
          do j=0,sizes(i)-1   ! fourier uses one less point; it's periodic
             grid(i,j) = j*2.0_r8/sizes(i) - 1.0_r8
          end do
       end if
    end do

!    print *, grid(1,20)


    do i=1,ndim
       do j=0,sizes(i)/2
          deltas(i,j) = abs(grid(i,j+1) - grid(i,j))
          deltas(i,sizes(i)-j) = deltas(i,j)
       enddo
    enddo

!    print *, grid(1,20)

    ! this is a kludge.  figure out how to do this better.
!    dx = grid(1,1) - grid(1,0)
!    if (dx>grid(2,1)-grid(2,0)) dx = grid(2,1)-grid(2,0)

    return
  end subroutine init_grid
end module Geometry
