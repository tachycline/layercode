
module Diagnostics
  use Types
  use Geometry
!  use Equations
  use newtrans
  private

  public :: compute_energy,compute_ape

  logical,public :: get_energy

  integer,private :: startcount,endcount
  real(kind=r8),public :: minmax1,minmax2,btminmax,bcminmax
  public :: start_timing,report_timing
  real(kind=r8),public :: enstrophy
  
  public :: compute_enstrophy
  

contains
  subroutine compute_energy(psi,q,totale,ave,ske,deltapsi)
     real(kind=r8),dimension(0:,0:),intent(in) :: psi,q
     real(kind=r8),intent(out) :: totale,ave,ske,deltapsi
     integer :: nx,ny,ny2,i,j
     real(kind=r8),dimension(:,:),allocatable :: psihat


!     print *, "computing energy"
     nx = size(psi,1) - 2
     ny = size(psi,2) - 2
     ny2 = ny/2

     allocate(psihat(0:nx-1,0:ny-1))

     do j=0,ny-1
        do i=0,nx-1
           psihat(i,j) = -q(i,j)*psi(i,j)*deltas(1,i)*deltas(2,j)
        enddo
     enddo

     totale = sum(psihat)
     ske = sum(psihat(0:nx,0:ny2))

!     print *, "energy is",totale

! normalize both to full basin size
     ave = totale/(nx*ny)

     deltapsi=(maxval(psi)+minval(psi))/maxval(abs(psi))

     deallocate(psihat)

     return
  end subroutine compute_energy

  subroutine compute_ape(psi1,psi2,lambdasq,delta,ape,avape,sape,avsape)
     real(kind=r8),dimension(0:,0:),intent(in) :: psi1,psi2
     real(kind=r8),intent(in) :: lambdasq,delta
     real(kind=r8),intent(out) :: ape,avape,sape,avsape
     integer :: nx,ny,ny2

     nx = size(psi1,1) - 2
     ny = size(psi1,2) - 2
     ny2 = ny/2

     ape = sum(lambdasq*delta/(2*(1+delta)**2)*(psi2-psi1)**2)
     sape = sum(lambdasq*delta/(2*(1+delta)**2)*&
                 (psi2(0:nx,0:ny2)-psi1(0:nx,0:ny2))**2)


! normalize both to full basin size
     avape = ape/(nx*ny)
     avsape = sape/(nx*ny)

!! get minmax1 and minmax2
     minmax1 = maxval(psi1)+minval(psi1)
     minmax2 = maxval(psi2)+minval(psi2)

     return
  end subroutine compute_ape

  subroutine start_timing()
  
    call system_clock(count=startcount)
  
    return
  end subroutine start_timing
  
  subroutine report_timing(numsteps)
    integer,optional,intent(in) :: numsteps
    integer :: elapsed
    real(kind=r8) :: average
  
    call system_clock(count=endcount)
  
    elapsed = endcount - startcount
    print *,"elapsed time: ",elapsed,"ms"
    if(present(numsteps)) then
      average = elapsed/numsteps
      print *,"average time per step: ",average
    endif
  
    return
  end subroutine report_timing
  
    subroutine compute_enstrophy(q,enstr)
       real(kind=r8),dimension(0:,0:),intent(in) :: q
       real(kind=r8),intent(out) :: enstr
  
       enstr = sum(q*q)/size(q)
  
       return
    end subroutine compute_enstrophy
  

end module Diagnostics
