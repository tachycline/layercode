module NewTrans
  use Types

  private

  
  real(kind=r8),dimension(:),allocatable,save,private :: ax,bx,ay,by,x_2d,y_2d
  integer,private :: nx,ny,nx2,ny2
  integer,public :: nxtrunc,nytrunc
  !private :: ax,bx,ay,by,x_2d,y_2d,nx,ny,nx2,ny2
  
  public :: tcxcy,fcxcy
  
  public :: init_trans
  real(kind=r8),dimension(:),allocatable,save,private :: work_2d_x,table_2d_x
  real(kind=r8),dimension(:),allocatable,save,private :: work_2d_y,table_2d_y
  integer,dimension(0:1),private :: isys
  real(kind=r8),private :: pi,scale_2d,dummy
  !private :: dummy
  !private :: pi,isys,scale_2d,work_2d_x,work_2d_y,table_2d_x,table_2d_y
  
  public :: d_dx_cheb, d_dy_cheb
  
  ! matricies for the similarity transform diagonalizing a,b
    real(kind=r8),dimension(:,:),allocatable,save,public :: &
           trans_q,trans_qinv,p,pinv
  
  ! eigenvalues of a, b
    real(kind=r8),dimension(:),allocatable,save,public:: wra,wia,wrb,wib
  
  ! original derivative matricies (and smaller forms)
    real(kind=r8),dimension(:,:),allocatable,save,private :: b,a,bs,as
  
  ! eigenvectors (for lapack)
    real(kind=r8),dimension(:,:),allocatable,save,private :: vr,vl
  
  ! work array (for lapack)
    real(kind=r8),dimension(:),allocatable,save,private :: work
  
  ! parameters for lapack routines
    integer,private :: lwork,info,ldh,ldvl,ldvr
  
  ! flags for lapack -- do we need left/right eigenvectors?
    character(len=1),private :: jobvl,jobvr
  
  ! record of pivots in decomposition
    integer,dimension(:),allocatable,save,private :: ipiv
  
  !private :: ipiv,jobvl,jobvr,lwork,info,ldh,ldvl,ldvr,work,vr,vl,b,a,as,bs
  !public :: wra,wia,wrb,wib,trans_q,trans_qinv,p,pinv
  
  public :: solve_p, solve_p_nondim,solve_p_nondim_split
  real(kind=r8),dimension(:,:,:),allocatable,save,private :: lhs,lulhs
  real(kind=r8),dimension(:,:),allocatable,save,private :: rhs,terms,wy
  real(kind=r8),dimension(:),allocatable,save,private :: sum1,sum2
  
  !private :: lhs,rhs,lulhs,terms,sum1,sum2,wy
  
contains
  
  
  subroutine tcxcy(psi)
  !tcxcy means To Cheb in X, Cheb in Y
    real(kind=r8),dimension(0:,0:),intent(inout) :: psi
    integer :: i,j
    
    real(kind=r8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b1
    
  
    do j=0,ny
       
       ax(0:nx) = psi(0:nx, j)
       
       
       a0 = ax(nx2-1)+ax(nx2+1)
       a9 = ax(1)
       do i=2,nx2-2,2
          a0 = a0+a9+ax(nx+1-i)
          a1 = ax(  i+1) - a9
          a2 = ax(nx+1-i) - ax(nx-1-i)
          a3 = ax(i)+ax(nx-i)
          a4 = ax(i)-ax(nx-i)
          a5 = a1 - a2
          a6 = a1 + a2
          a7 = x_2d(nx2-i)*a4+x_2d(i)*a6
          a8 = x_2d(nx2-i)*a6-x_2d(i)*a4
          a9 = ax(i+1)
          bx(i  ) = a3+a7
          bx(nx-i) = a3-a7
          bx(i+1  ) = a8+a5
          bx(nx+1-i) = a8-a5
       end do
       bx(1) = ax(0)-ax(nx)
       bx(0) = ax(0)+ax(nx)
       bx(nx2) = 2.0*ax(nx2)
       bx(nx2+1)=2.0*(a9-ax(nx2+1))
       
       
       ! linux version
       call fft(bx(0),bx(1),2,nx2,1)
       
       ! sgi version
       !call zzfft(1,nx2,scale_2d,bx,bx,table_2d_x,work_2d_x,isys)
       
       
       
       a0 = 2.0*a0
       bx(nx) = bx(0)-a0
       bx(0) = bx(0)+a0
       do i=1,nx2-1
          a1 =  0.5           *(bx(i)+bx(nx-i))
          a2 = 0.25/x_2d(nx2+i)*(bx(i)-bx(nx-i))
          bx(i) = a1+a2
          bx(nx-i) = a1-a2
       end do
       
       b1 = 0.5_r8/nx
       bx(0) = bx(0)*b1
       bx(nx) = bx(nx)*b1
       b1 = 2.0_r8*b1
       do i=2,nx-2,2
          bx(i) = bx(i)*b1
       enddo
       b1 = -b1
       do i=1,nx-1,2
          bx(i) = bx(i)*b1
       enddo
       
       
       psi(0:nx,j) = bx(0:nx)
       
    enddo
  
  !  do j=0,nxtrunc-1
    do j=0,nx
       
       ay(0:ny) = psi(j,0:ny)
       
       
       a0 = ay(ny2-1)+ay(ny2+1)
       a9 = ay(1)
       do i=2,ny2-2,2
          a0 = a0+a9+ay(ny+1-i)
          a1 = ay(  i+1) - a9
          a2 = ay(ny+1-i) - ay(ny-1-i)
          a3 = ay(i)+ay(ny-i)
          a4 = ay(i)-ay(ny-i)
          a5 = a1 - a2
          a6 = a1 + a2
          a7 = y_2d(ny2-i)*a4+y_2d(i)*a6
          a8 = y_2d(ny2-i)*a6-y_2d(i)*a4
          a9 = ay(i+1)
          by(i  ) = a3+a7
          by(ny-i) = a3-a7
          by(i+1  ) = a8+a5
          by(ny+1-i) = a8-a5
       end do
       by(1) = ay(0)-ay(ny)
       by(0) = ay(0)+ay(ny)
       by(ny2) = 2.0*ay(ny2)
       by(ny2+1)=2.0*(a9-ay(ny2+1))
       
       
       ! linux version
       call fft(by(0),by(1),2,ny2,1)
       
       ! sgi version
       !call zzfft(1,ny2,scale_2d,by,by,table_2d_y,work_2d_y,isys)
       
       
       
       a0 = 2.0*a0
       by(ny) = by(0)-a0
       by(0) = by(0)+a0
       do i=1,ny2-1
          a1 =  0.5           *(by(i)+by(ny-i))
          a2 = 0.25/y_2d(ny2+i)*(by(i)-by(ny-i))
          by(i) = a1+a2
          by(ny-i) = a1-a2
       end do
       
       b1 = 0.5_r8/ny
       by(0) = by(0)*b1
       by(ny) = by(ny)*b1
       b1 = 2.0_r8*b1
       do i=2,ny-2,2
          by(i) = by(i)*b1
       enddo
       b1 = -b1
       do i=1,ny-1,2
          by(i) = by(i)*b1
       enddo
       
       
       psi(j,0:ny) = by(0:ny)
       
    enddo
  
    psi(nxtrunc-1:,:) = 0.0_r8
    psi(:,nytrunc-1:) = 0.0_r8
  
    return
  end subroutine tcxcy
  
  
  subroutine fcxcy(psi)
  !fcxcy means From Cheb in X, Cheb in Y
    real(kind=r8),dimension(0:,0:),intent(inout) :: psi
    integer :: i,j
    
    real(kind=r8) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,b1
    
  
    psi(nxtrunc-1:,:) = 0.0_r8
    psi(:,nytrunc-1:) = 0.0_r8
  
  !  do j=0,nxtrunc-1
    do j=0,nx
       
       ay(0:ny) = psi(j,0:ny)
       
       
       !ay(0) = ay(0)
       b1 = 0.5_r8
       do i=2,ny-2,2
          ay(i) = b1*ay(i)
       enddo
       b1 = -b1
       do i=1,ny-1,2
          ay(i) = b1*ay(i)
       enddo
       !by(ny) = ay(ny)
       
       a0 = ay(ny2-1)+ay(ny2+1)
       a9 = ay(1)
       do i=2,ny2-2,2
          a0 = a0+a9+ay(ny+1-i)
          a1 = ay(  i+1) - a9
          a2 = ay(ny+1-i) - ay(ny-1-i)
          a3 = ay(i)+ay(ny-i)
          a4 = ay(i)-ay(ny-i)
          a5 = a1 - a2
          a6 = a1 + a2
          a7 = y_2d(ny2-i)*a4+y_2d(i)*a6
          a8 = y_2d(ny2-i)*a6-y_2d(i)*a4
          a9 = ay(i+1)
          by(i  ) = a3+a7
          by(ny-i) = a3-a7
          by(i+1  ) = a8+a5
          by(ny+1-i) = a8-a5
       end do
       by(1) = ay(0)-ay(ny)
       by(0) = ay(0)+ay(ny)
       by(ny2) = 2.0*ay(ny2)
       by(ny2+1)=2.0*(a9-ay(ny2+1))
       
       
       
       ! linux version
       call fft(by(0),by(1),2,ny2,1)
       
       ! sgi version
       !call zzfft(1,ny2,scale_2d,by,by,table_2d_y,work_2d_y,isys)
       
       
       a0 = 2.0*a0
       by(ny) = by(0)-a0
       by(0) = by(0)+a0
       do i=1,ny2-1
          a1 =  0.5           *(by(i)+by(ny-i))
          a2 = 0.25/y_2d(ny2+i)*(by(i)-by(ny-i))
          by(i) = a1+a2
          by(ny-i) = a1-a2
       end do
       
       
       psi(j,0:ny) = by(0:ny)
       
    enddo
  
    do j=0,ny
       
       ax(0:nx) = psi(0:nx, j)
       
       
       !ax(0) = ax(0)
       b1 = 0.5_r8
       do i=2,nx-2,2
          ax(i) = b1*ax(i)
       enddo
       b1 = -b1
       do i=1,nx-1,2
          ax(i) = b1*ax(i)
       enddo
       !bx(nx) = ax(nx)
       
       a0 = ax(nx2-1)+ax(nx2+1)
       a9 = ax(1)
       do i=2,nx2-2,2
          a0 = a0+a9+ax(nx+1-i)
          a1 = ax(  i+1) - a9
          a2 = ax(nx+1-i) - ax(nx-1-i)
          a3 = ax(i)+ax(nx-i)
          a4 = ax(i)-ax(nx-i)
          a5 = a1 - a2
          a6 = a1 + a2
          a7 = x_2d(nx2-i)*a4+x_2d(i)*a6
          a8 = x_2d(nx2-i)*a6-x_2d(i)*a4
          a9 = ax(i+1)
          bx(i  ) = a3+a7
          bx(nx-i) = a3-a7
          bx(i+1  ) = a8+a5
          bx(nx+1-i) = a8-a5
       end do
       bx(1) = ax(0)-ax(nx)
       bx(0) = ax(0)+ax(nx)
       bx(nx2) = 2.0*ax(nx2)
       bx(nx2+1)=2.0*(a9-ax(nx2+1))
       
       
       
       ! linux version
       call fft(bx(0),bx(1),2,nx2,1)
       
       ! sgi version
       !call zzfft(1,nx2,scale_2d,bx,bx,table_2d_x,work_2d_x,isys)
       
       
       a0 = 2.0*a0
       bx(nx) = bx(0)-a0
       bx(0) = bx(0)+a0
       do i=1,nx2-1
          a1 =  0.5           *(bx(i)+bx(nx-i))
          a2 = 0.25/x_2d(nx2+i)*(bx(i)-bx(nx-i))
          bx(i) = a1+a2
          bx(nx-i) = a1-a2
       end do
       
       
       psi(0:nx,j) = bx(0:nx)
       
    enddo
  
    return
  end subroutine fcxcy
  
  
  subroutine init_trans(nnx,nny)
     integer,intent(in) :: nnx,nny
     integer :: i,j
  
     real(kind=r8) :: bk,bkp2,ckm2
     
  
     nx = nnx
     ny = nny
     nx2 = nx/2
     ny2 = ny/2
  
  ! put in logic here, if not dealiasing.
     nxtrunc = 2*nx/3
     nytrunc = 2*ny/3
  
     if(iand(nxtrunc,1) == 1) then ! nxtrunc is odd
       nxtrunc = nxtrunc - 1
     endif
     if(iand(nytrunc,1) == 1) then ! nytrunc is odd
       nytrunc = nytrunc - 1
     endif
  
     ! should I worry about making these allocations safe?
  
     allocate(ax(0:nx),bx(0:nx),x_2d(0:nx))
     allocate(ay(0:ny),by(0:ny),y_2d(0:ny))
  
     pi = 4.0_r8 * atan(1.0_r8)
     do i=0,nx
        x_2d(i) = -cos(i*pi/nx)
     enddo
     do i=0,ny
        y_2d(i) = -cos(i*pi/ny)
     enddo
  
     ! x derivative matrix--------------------------------------------
     
     allocate(a(0:nxtrunc+4,0:nxtrunc+4),as(nxtrunc+3,nxtrunc+3))
     allocate(work(nxtrunc+5),wra(nxtrunc+3),wia(nxtrunc+3))
     allocate(pinv(nxtrunc+3,nxtrunc+3),p(nxtrunc+3,nxtrunc+3))
     allocate(vr(nxtrunc+3,nxtrunc+3),vl(nxtrunc+3,nxtrunc+3))
     allocate(ipiv(nxtrunc+3))
     
     a = 0.0_r8
     do j = 2,nxtrunc+4,2
        a(0,j) = j*j*j/2.0_r8
     enddo
     do i=1,nxtrunc+3
        do j = i+2,nxtrunc+4,2
           a(i,j) = j*(j*j-i*i)
        enddo
     enddo
     
     !eliminate highest order coefficients with bc
     do i=0,nxtrunc+2
        a(i,0:nxtrunc+2:2) = a(i,0:nxtrunc+2:2) - a(i,nxtrunc+4)
        a(i,1:nxtrunc+1:2) = a(i,1:nxtrunc+1:2) - a(i,nxtrunc+3)
     enddo
     
     ! put the result in the smaller matrix
     as = a(0:nxtrunc+2,0:nxtrunc+2)
     
     ! get eigenvalues,eigenvectors
     ldh = nxtrunc+3
     ldvl = nxtrunc+3
     ldvr = nxtrunc+3
     
     ! we want right eigenvectors, don't need left eigenvectors.
     jobvl = "N"
     jobvr = "V"
     ! get optimal size for work, set lwork=-1
     ! this doesn't work on rossby
     !lwork = -1
     !call dgeev(jobvl,jobvr,nx-1,as,nx-1,wra,wia, &
     !                    vl,ldvl,vr,ldvr,work,lwork,info)
     !lwork = work(1)
     
     lwork = 4*(nxtrunc+4)**2
     deallocate(work)
     allocate(work(lwork))
     ! actual call to get eigenvalues, eigenvectors
     call dgeev(jobvl,jobvr,nxtrunc+3,as,nxtrunc+3,wra,wia, &
                         vl,ldvl,vr,ldvr,work,lwork,info)
     p = vr(1:nxtrunc+3,1:nxtrunc+3)
     
     ! now invert matrix of eigenvectors
     pinv = p
     call dgetrf(nxtrunc+3,nxtrunc+3,pinv,nxtrunc+3,ipiv,info)
     
     ! again, the lwork = -1 trick doesn't work on rossby
     ! I think I made it big enough last time around.
     !lwork = -1
     !call dgetri(nx-1,pinv,nx-1,ipiv,work,lwork,info)
     !lwork = work(1)
     !deallocate(work)
     !allocate(work(lwork))
     
     call dgetri(nxtrunc+3,pinv,nxtrunc+3,ipiv,work,lwork,info)
     
     deallocate(ipiv,vr,vl,work)
     
     
     ! y derivative matrix--------------------------------------------
     
     allocate(b(0:nytrunc+4,0:nytrunc+4),bs(nytrunc+3,nytrunc+3))
     allocate(work(nytrunc+5),wrb(nytrunc+3),wib(nytrunc+3))
     allocate(trans_qinv(nytrunc+3,nytrunc+3),trans_q(nytrunc+3,nytrunc+3))
     allocate(vr(nytrunc+3,nytrunc+3),vl(nytrunc+3,nytrunc+3))
     allocate(ipiv(nytrunc+3))
     
     b = 0.0_r8
     do j = 2,nytrunc+4,2
        b(0,j) = j*j*j/2.0_r8
     enddo
     do i=1,nytrunc+3
        do j = i+2,nytrunc+4,2
           b(i,j) = j*(j*j-i*i)
        enddo
     enddo
     
     !eliminate highest order coefficients with bc
     do i=0,nytrunc+2
        b(i,0:nytrunc+2:2) = b(i,0:nytrunc+2:2) - b(i,nytrunc+4)
        b(i,1:nytrunc+1:2) = b(i,1:nytrunc+1:2) - b(i,nytrunc+3)
     enddo
     
     ! put the result in the smaller matrix
     bs = b(0:nytrunc+2,0:nytrunc+2)
     
     b = transpose(b)
     bs = transpose(bs)
     ! get eigenvalues,eigenvectors
     ldh = nytrunc+3
     ldvl = nytrunc+3
     ldvr = nytrunc+3
     
     ! we want right eigenvectors, don't need left eigenvectors.
     jobvl = "N"
     jobvr = "V"
     ! get optimal size for work, set lwork=-1
     ! this doesn't work on rossby
     !lwork = -1
     !call dgeev(jobvl,jobvr,ny-1,bs,ny-1,wrb,wib, &
     !                    vl,ldvl,vr,ldvr,work,lwork,info)
     !lwork = work(1)
     
     lwork = 4*(nytrunc+4)**2
     deallocate(work)
     allocate(work(lwork))
     ! actual call to get eigenvalues, eigenvectors
     call dgeev(jobvl,jobvr,nytrunc+3,bs,nytrunc+3,wrb,wib, &
                         vl,ldvl,vr,ldvr,work,lwork,info)
     trans_q = vr(1:nytrunc+3,1:nytrunc+3)
     
     ! now invert matrix of eigenvectors
     trans_qinv = trans_q
     call dgetrf(nytrunc+3,nytrunc+3,trans_qinv,nytrunc+3,ipiv,info)
     
     ! again, the lwork = -1 trick doesn't work on rossby
     ! I think I made it big enough last time around.
     !lwork = -1
     !call dgetri(ny-1,trans_qinv,ny-1,ipiv,work,lwork,info)
     !lwork = work(1)
     !deallocate(work)
     !allocate(work(lwork))
     
     call dgetri(nytrunc+3,trans_qinv,nytrunc+3,ipiv,work,lwork,info)
     
     deallocate(ipiv,vr,vl,work)
     
     
     if (allocated(lhs)) deallocate(lhs)
     allocate(lhs(1:3,0:nx,0:ny))
     if (allocated(rhs)) deallocate(rhs)
     allocate(rhs(0:nx,0:ny))
     if (allocated(lulhs)) deallocate(lulhs)
     allocate(lulhs(1:3,0:nx,0:ny))
     if (allocated(terms)) deallocate(terms)
     allocate(terms(3,0:nx))
     if (allocated(sum1)) deallocate(sum1)
     allocate(sum1(0:ny))
     if (allocated(sum2)) deallocate(sum2)
     allocate(sum2(0:ny))
     if (allocated(wy)) deallocate(wy)
     allocate(wy(0:nx,0:ny))
     
     !do i=2,nx
     !   bk = 1.0_r8
     !   if(i>nx-2) bk = 0.0_r8
     !   bkp2 = 1.0_r8
     !   if(i>(nx-4)) bkp2 = 0.0_r8
     !   ckm2 = 1.0_r8
     !   if(i==2) ckm2 = 2.0_r8
     
     !   terms(1,i) = ckm2/(4.0_r8*i*(i-1.0_r8))
     !   terms(2,i) =   bk/(2.0_r8*(i*i-1.0_r8))
     !   terms(3,i) = bkp2/(4.0_r8*i*(i+1.0_r8))
     !enddo
     
     !do i=2,nx
     !   lhs(1,i,0:ny-2) =  terms(1,i)*wrb(1:ny-1)
     !   lhs(2,i,0:ny-2) = -terms(2,i)*wrb(1:ny-1) + 1.0_r8
     !   lhs(3,i,0:ny-2) =  terms(3,i)*wrb(1:ny-1)
     !enddo
     !lhs(:,:,ny-1) = 0.0_r8
     !lhs(2,:,ny-1) = 1.0_r8
     
     !lulhs(1,nx,0:ny-1) = lhs(2,nx,0:ny-1)
     !lulhs(1,nx-1,0:ny-1) = lhs(2,nx-1,0:ny-1)
     !lulhs(3,nx,0:ny-1) = 1.0_r8/lulhs(1,nx,0:ny-1)
     !lulhs(3,nx-1,0:ny-1) = 1.0_r8/lulhs(1,nx-1,0:ny-1)
     !do i=nx-2,2,-1
     !   lulhs(2,i,0:ny-1) = lhs(3,i,0:ny-1)/lulhs(1,i+2,0:ny-1)
     !   lulhs(1,i,0:ny-1) = lhs(2,i,0:ny-1) - lulhs(2,i,0:ny-1)&
     !                   *lhs(1,i+2,0:ny-1)
     !   lulhs(3,i,0:ny-1) = (1.0_r8-lulhs(3,i+2,0:ny-1)*lhs(1,i+2,0:ny-1))/&
     !                        lulhs(1,i,0:ny-1)
     !enddo
     !lulhs(1,0,0:ny-1) = 1.0_r8-lulhs(3,2,0:ny-1)*lhs(1,2,0:ny-1)
     !lulhs(1,1,0:ny-1) = 1.0_r8-lulhs(3,3,0:ny-1)*lhs(1,3,0:ny-1)
     
     
  
     allocate(work_2d_x(nx),table_2d_x(nx+256))
     allocate(work_2d_y(ny),table_2d_y(ny+256))
     isys(0) = 1
     scale_2d = 1.0
     !call zzfft(0,nx/2,scale_2d,dummy,dummy,table_2d_x,work_2d_x,isys)
     !call zzfft(0,ny/2,scale_2d,dummy,dummy,table_2d_y,work_2d_y,isys)
     
  
     return
  end subroutine init_trans
  
  
  
  subroutine d_dy_cheb(psi)
     real(kind=r8),dimension(0:,0:),intent(inout) :: psi
     real(kind=r8),dimension(0:size(psi,2)-1) :: a1,a2
     real(kind=r8),dimension(0:size(psi,1)-1,0:size(psi,2)-1) :: psiy
     integer :: i
  
     ! coefficients for psi' are given by c'_{i-1} = c'_{i+1} + 2(i)c_i
     ! since we don't have c_{nx+1},
     a1(:) = psi(:,ny)
     a2(:) = psi(:,ny-1)
     psiy(:,ny) = 0.0_r8
     psiy(:,ny-1) = 2.0_r8*(ny)*a1(:)
     a1 = a2
     a2 = psi(:,ny-2)
     psiy(:,ny-2) = 2.0_r8*(ny-1)*a1(:)
  
     do i=ny-2,2,-1
        a1 = a2
        a2 = psi(:,i-1)
        psiy(:,i-1) = psiy(:,i+1)+2.0_r8*(i)*a1(:)
     end do
  
     psiy(:,0) = 0.5_r8*psiy(:,2)+a2
  
     psi = psiy
  
     return
  end subroutine d_dy_cheb
  
  
  subroutine d_dx_cheb(psi)
     real(kind=r8),dimension(0:,0:),intent(inout) :: psi
     real(kind=r8),dimension(0:size(psi,2)-1) :: a1,a2
     real(kind=r8),dimension(0:size(psi,1)-1,0:size(psi,2)-1) :: psix
     integer :: i
  
     ! coefficients for psi' are given by c'_{i-1} = c'_{i+1} + 2(i)c_i
     ! since we don't have c_{nx+1},
     a1(:) = psi(nx,:)
     a2(:) = psi(nx-1,:)
     psix(nx,:) = 0.0_r8
     psix(nx-1,:) = 2.0_r8*(nx)*a1(:)
     a1 = a2
     a2 = psi(nx-2,:)
     psix(nx-2,:) = 2.0_r8*(nx-1)*a1(:)
  
     do i=nx-2,2,-1
        a1 = a2
        a2 = psi(i-1,:)
        psix(i-1,:) = psix(i+1,:)+2.0_r8*(i)*a1(:)
     end do
  
     psix(0,:) = 0.5_r8*psix(2,:)+a2
  
     psi = psix
  
     return
  end subroutine d_dx_cheb
  
  
  
  subroutine solve_p(f,u,lambda)
    real(kind=r8),dimension(0:,0:),intent(in) :: f
    real(kind=r8),dimension(0:,0:),intent(out) :: u
    real(kind=r8),intent(in) :: lambda
  
    real(kind=r8),dimension(0:nxtrunc+2,&
                       0:nytrunc+2) :: uprime
    real(kind=r8),dimension(0:nxtrunc+2,&
                       0:nytrunc+2) :: fprime
    real(kind=r8) :: a2
    integer :: i,j
  
  !  print *, "size(uprime,1):",size(uprime,1)
  !  print *, "size(trans_q,1):",size(trans_q,1)
  !  print *, "nxtrunc: ",nxtrunc,"size(fprime,1):",size(fprime,1)
  !  print *, "nxtrunc: ",nxtrunc,"size(f,1):",size(f,1)
  
  
    fprime = f(0:nxtrunc+2,0:nytrunc+2)
    fprime = matmul(pinv,matmul(fprime,trans_q))
    do i=0,nxtrunc+2
       do j=0,nytrunc+2
          uprime(i,j) = fprime(i,j)/(wra(i+1) + wrb(j+1) - lambda)
       enddo
    enddo
    uprime = matmul(p,matmul(uprime,trans_qinv))
    u(0:nxtrunc+2,0:nytrunc+2) = uprime
    u(nxtrunc-1:,:) = 0.0_r8
    u(:,nytrunc-1:) = 0.0_r8
  
    return
  end subroutine solve_p
  
  
  subroutine solve_p_nondim(f,u,lambda,alpha)
    real(kind=r8),dimension(0:,0:),intent(in) :: f
    real(kind=r8),dimension(0:,0:),intent(out) :: u
    real(kind=r8),intent(in) :: lambda,alpha
  
    real(kind=r8),dimension(0:nxtrunc+2,&
                       0:nytrunc+2) :: uprime
    real(kind=r8),dimension(0:nxtrunc+2,&
                       0:nytrunc+2) :: fprime
    real(kind=r8) :: a2
    integer :: i,j
  
    a2 = alpha**2
  
  !  print *, "size(uprime,1):",size(uprime,1)
  !  print *, "size(trans_q,1):",size(trans_q,1)
  !  print *, "nxtrunc: ",nxtrunc,"size(fprime,1):",size(fprime,1)
  !  print *, "nxtrunc: ",nxtrunc,"size(f,1):",size(f,1)
  
    fprime = f(0:nxtrunc+2,0:nytrunc+2)
    fprime = matmul(pinv,matmul(fprime,trans_q))
    do i=0,nxtrunc+2
       do j=0,nytrunc+2
          uprime(i,j) = fprime(i,j)/(a2*wra(i+1) + wrb(j+1) - lambda)
       enddo
    enddo
    uprime = matmul(p,matmul(uprime,trans_qinv))
    u(0:nxtrunc+2,0:nytrunc+2) = uprime
    u(nxtrunc-1:,:) = 0.0_r8
    u(:,nytrunc-1:) = 0.0_r8
  
    return
  end subroutine solve_p_nondim
  
  
  subroutine solve_p_nondim_split(f,u,lambda,alpha)
    real(kind=r8),dimension(0:,0:),intent(in) :: f
    real(kind=r8),dimension(0:,0:),intent(out) :: u
    real(kind=r8),dimension(0:size(f,1)-1,0:size(f,2)-1) :: ftmp,utmp
    real(kind=r8),intent(in) :: lambda,alpha
  
    real(kind=r8),dimension(0:2*size(u,1)/3-3,0:2*size(u,2)/3-3) :: uprime
    real(kind=r8),dimension(0:2*size(f,1)/3-3,0:2*size(f,2)/3-3) :: fprime
    real(kind=r8) :: a2
    integer :: i
  
    a2 = alpha**2
    fprime = f(0:nxtrunc-2,0:nytrunc-2)
    fprime = matmul(fprime,trans_q)
    
  
    if(lambda == 0.0_r8) then
      ftmp(0:nxtrunc-2,0:nytrunc-2) = fprime
      ftmp(nxtrunc-1:,:) = 0.0_r8
      ftmp(:,nytrunc-1:) = 0.0_r8
      do i=2,nx-2
         rhs(i,0:ny-1) = ftmp(i-2,0:ny-1)*terms(1,i)-ftmp(i,0:ny-1)*terms(2,i)&
                        +ftmp(i+2,0:ny-1)*terms(3,i)
      enddo
      
      do i=nx-1,nx
         rhs(i,0:ny-1) = ftmp(i-2,0:ny-1)*terms(1,i)-ftmp(i,0:ny-1)*terms(2,i)
      enddo
      
      wy(nx,0:ny-1) = rhs(nx,0:ny-1)
      wy(nx-1,0:ny-1) = rhs(nx-1,0:ny-1)
      do i=nx-2,2,-1
         wy(i,0:ny-1) = rhs(i,0:ny-1)-lulhs(2,i,0:ny-1)*wy(i+2,0:ny-1)
      enddo
      sum1(0:ny-1) = 0.0_r8
      sum2(0:ny-1) = 0.0_r8
      do i=2,nx-2,2
         sum1(0:ny-1) = sum1(0:ny-1)+lulhs(3,i,0:ny-1)*wy(i,0:ny-1)
         sum2(0:ny-1) = sum2(0:ny-1)+lulhs(3,i+1,0:ny-1)*wy(i+1,0:ny-1)
      enddo
      sum1(0:ny-1) = sum1(0:ny-1)+lulhs(3,nx,0:ny-1)*wy(nx,0:ny-1)
      wy(0,0:ny-1) = -sum1(0:ny-1)
      wy(1,0:ny-1) = -sum2(0:ny-1)
      
      utmp = 0.0_r8
      utmp(0,0:ny-1)=wy(0,0:ny-1)/lulhs(1,0,0:ny-1)
      utmp(1,0:ny-1)=wy(1,0:ny-1)/lulhs(1,1,0:ny-1)
      do i=2,nx
         utmp(i,0:ny-1) = (wy(i,0:ny-1)-lhs(1,i,0:ny-1)*utmp(i-2,0:ny-1))/&
                        lulhs(1,i,0:ny-1)
      enddo
      
      
      uprime = utmp(0:nxtrunc-2,0:nytrunc-2)
      
    else
    endif
  
    uprime = matmul(uprime,trans_qinv)
    u(0:nxtrunc-2,0:nytrunc-2) = uprime
    u(nxtrunc-1:,:) = 0.0_r8
    u(:,nytrunc-1:) = 0.0_r8
    
    return
  end subroutine solve_p_nondim_split
  
  
end module NewTrans
