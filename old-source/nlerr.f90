program nlerr
  use Types
  use Geometry
  use NewTrans
  use Control
  use Layers
  use Integrator
  use Initialization
  use Diagnostics

  type(layer) :: mode1,mode2,suml,calc
  real(kind=r8) :: al1,al2,be1,be2,amp1,amp2,thick,dens

  call get_run_info()

  nx = param("nx",128,0)
  ny = param("ny",128,0)
  amp1 = param("amp1",1.0_r8,0)
  amp2 = param("amp2",1.0_r8,0)
  al1 = param("a1",1.0_r8,0)
  al2 = param("a2",1.0_r8,0)
  be1 = param("b1",1.0_r8,0)
  be2 = param("b2",1.0_r8,0)
  thick = param("thick1",1000.0_r8,0)
  dens = param("dens1",1021.0_r8,0)

  call init_geometry(nnx=nx,nny=ny,chebx=.true.,cheby=.true.)
! set up functional pieces
  call init_nondim()
  call init_trans(nx,ny)
  call init_force(ngyres)
  call init_integrator()

  call makelayer(mode1,nx,ny,thick,dens)
  call makelayer(mode2,nx,ny,thick,dens)
  call makelayer(suml,nx,ny,thick,dens)
  call makelayer(calc,nx,ny,thick,dens)

  call init_mode(mode1,al1,al2,amp1)
  call init_mode(mode2,be1,be2,amp2)

  suml%psi = mode1%psi + mode2%psi
  suml%dx = mode1%dx + mode2%dx
  suml%dy = mode1%dy + mode2%dy
  suml%q = mode1%q + mode2%q
  suml%ddx = mode1%ddx + mode2%ddx
  suml%ddy = mode1%ddy + mode2%ddy

! the analytic nonlinear term
  suml%k1 = suml%dx*suml%ddy - suml%dy*suml%ddx

! now, the calculations...
  calc%q = suml%q
  call tcxcy(calc%q)
  call solve_p_nondim(calc%q,calc%psi,0.0_r8,alpha)
  calc%dx = calc%psi
  calc%dy = calc%psi
  call d_dx_cheb(calc%dx)
  call d_dy_cheb(calc%dy)

  calc%ddx = calc%q
  calc%ddy = calc%q
  call d_dx_cheb(calc%ddx)
  call d_dy_cheb(calc%ddy)

  call fcxcy(calc%psi)
  call fcxcy(calc%q)
  call fcxcy(calc%dx)
  call fcxcy(calc%dy)
  call fcxcy(calc%ddx)
  call fcxcy(calc%ddy)

  calc%k1 = calc%dx*calc%ddy - calc%dy*calc%ddx

! use plplot instead for plotting fields

!  call plot_field("calc_psi",calc%psi,0,0.0_r8)
!  call plot_field("suml_psi",suml%psi,0,0.0_r8)
  print *,"l2 norm, calc%psi,suml%psi:",sqrt(sum((calc%psi-suml%psi)**2))
  
  print *,"linf norm, calc%psi,suml%psi:",maxval(abs(calc%psi-suml%psi)),&
          "at",maxloc(abs(calc%psi-suml%psi))
  

!  call plot_field("calc_dx",calc%dx,0,0.0_r8)
!  call plot_field("suml_dx",suml%dx,0,0.0_r8)
  print *,"l2 norm, calc%dx,suml%dx:",sqrt(sum((calc%dx-suml%dx)**2))
  
  print *,"linf norm, calc%dx,suml%dx:",maxval(abs(calc%dx-suml%dx)),&
          "at",maxloc(abs(calc%dx-suml%dx))
  

!  call plot_field("calc_dy",calc%dy,0,0.0_r8)
!  call plot_field("suml_dy",suml%dy,0,0.0_r8)
  print *,"l2 norm, calc%dy,suml%dy:",sqrt(sum((calc%dy-suml%dy)**2))
  
  print *,"linf norm, calc%dy,suml%dy:",maxval(abs(calc%dy-suml%dy)),&
          "at",maxloc(abs(calc%dy-suml%dy))
  

!  call plot_field("calc_q",calc%q,0,0.0_r8)
!  call plot_field("suml_q",suml%q,0,0.0_r8)
  print *,"l2 norm, calc%q,suml%q:",sqrt(sum((calc%q-suml%q)**2))
  
  print *,"linf norm, calc%q,suml%q:",maxval(abs(calc%q-suml%q)),&
          "at",maxloc(abs(calc%q-suml%q))
  

!  call plot_field("calc_ddx",calc%ddx,0,0.0_r8)
!  call plot_field("suml_ddx",suml%ddx,0,0.0_r8)
  print *,"l2 norm, calc%ddx,suml%ddx:",sqrt(sum((calc%ddx-suml%ddx)**2))
  
  print *,"linf norm, calc%ddx,suml%ddx:",maxval(abs(calc%ddx-suml%ddx)),&
          "at",maxloc(abs(calc%ddx-suml%ddx))
  

!  call plot_field("calc_ddy",calc%ddy,0,0.0_r8)
!  call plot_field("suml_ddy",suml%ddy,0,0.0_r8)
  print *,"l2 norm, calc%ddy,suml%ddy:",sqrt(sum((calc%ddy-suml%ddy)**2))
  
  print *,"linf norm, calc%ddy,suml%ddy:",maxval(abs(calc%ddy-suml%ddy)),&
          "at",maxloc(abs(calc%ddy-suml%ddy))
  

!  call plot_field("calc_k1",calc%k1,0,0.0_r8)
!  call plot_field("suml_k1",suml%k1,0,0.0_r8)
  print *,"l2 norm, calc%k1,suml%k1:",sqrt(sum((calc%k1-suml%k1)**2))
  
  print *,"linf norm, calc%k1,suml%k1:",maxval(abs(calc%k1-suml%k1)),&
          "at",maxloc(abs(calc%k1-suml%k1))
  

  !output animation
  !call dx_write_anim()
  close(unit=efnum)
  close(unit=ekenum)
  

  stop
end program nlerr
