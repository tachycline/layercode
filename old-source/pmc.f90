
program pmc
  use Types
  use Geometry
  use NewTrans
  use Control
  use Layers
  use Initialization
  use Diagnostics

  integer :: i,nincrements
  real(kind=r8) :: increment,thick2,dens1,dens2,delta,lambda,deltarho
  real(kind=r8) :: third
  character(len=80) :: exec_path,inc_var
  logical :: decrement

  
  ! get params
  call get_run_info()
  call init_nondim()
  
  exec_path = param("exec_path","/home/cmckay/work/code/build/betabt",0)
  inc_var = param("inc_var","deltai",0)
  increment = param("increment",0.001_r8,0)
  decrement = param("decrement",.false.,0)
  nincrements = param("nincrements",1,0)
  
  if (decrement) then
     increment = - increment
  endif
  
  ! write the initial param file
  nrun = nrun-1
  call writeparams()
  

  
  do i=1,nincrements
     ! get rid of the <basename>.nrun.param file
     print *, "rm big.0.param"
  
  !   call system("rm "//trim(adjustl(paramarchive)))
  
     ! run the executable
    print *, "running exec_path"
    call system(exec_path)
  
  ! read params (this updates nrun, nslab, and several other things)
     call purgeparams()
     call get_run_info()
     call init_nondim()
  
     ! increment the variable (do in terms of special cases)
  
     nrun = nrun - 1 ! to cope with the increment in writeparams()
  
     if (trim(adjustl(inc_var)) == "deltai") then
        lx = param("lx",1.0e6_r8,0)
        ly = param("ly",2.0e6_r8,0)
        beta = param("beta",2.0e-11_r8,0)
        fzero = param("fzero",1.0e-4_r8,0)
        w = param("w",1.0e-6_r8,0)
        thick1 = param("thick1",1.0e3_r8,0)
  
        deltai = sqrt(fzero*w/(beta**2*lx**2*thick1))
        deltai = deltai + increment
        w = deltai**2*beta**2*lx**2*h/fzero
        call modparam("w",w)
     endif
  
     if (trim(adjustl(inc_var)) == "deltam") then
        lx = param("lx",1.0e6_r8,0)
        beta = param("beta",2.0e-11_r8,0)
        kappa = param("kappa",1.0e2_r8,0)
  
        third = 1.0_r8/3.0_r8
        deltam = exp(third*log(kappa/beta))/lx
        deltam = deltam + increment
  
        kappa = deltam**3*beta*lx**3
        call modparam("kappa",kappa)
     endif
  
     if (trim(adjustl(inc_var)) == "ld") then
     nlayers = param("nlayers",1,0)
     select case (nlayers)
     case(1)
       thick1 = param("thick1",1000.0_r8,0)
       deltarho = param("deltarho",10.0_r8,0)
       rhozero = param("rhozero",1022.0_r8,0)
       fzero = param("fzero",1.0e-4_r8,0)
  
       rd = sqrt(9.8_r8*deltarho*thick1/(rhozero*fzero**2))
  
       rd = rd + increment
       deltarho = rd**2*fzero**2*rhozero/(9.8_r8*thick1)
       call modparam("deltarho",deltarho)
     case(2)
       thick1 = param("thick1",1000.0_r8,0)
       thick2 = param("thick2",4000.0_r8,0)
       dens1 = param("dens1",1022.0_r8,0)
       dens2 = param("dens2",1032.0_r8,0)
       fzero = param("fzero",1.0e-4_r8,0)
       rd = sqrt(9.8_r8*(dens2-dens1)*thick1*thick2/&
                     (dens1*fzero**2*(thick1+thick2)))
       rd = rd + increment
       dens2 = dens1+(rd**2*dens1*fzero**2*(thick1+thick2))/&
                                       (9.8_r8*thick1*thick2)
       call modparam("dens2",dens2)
     case default !bail out
       print *, "Wrong number of layers; aborting."
       stop
     end select
     endif
  
     if (trim(adjustl(inc_var)) == "lambda") then
     nlayers = param("nlayers",1,0)
     select case (nlayers)
     case(1)
       thick1 = param("thick1",1000.0_r8,0)
       deltarho = param("deltarho",10.0_r8,0)
       rhozero = param("rhozero",1022.0_r8,0)
       fzero = param("fzero",1.0e-4_r8,0)
       ly = param("ly",1.0e6_r8,0)
  
       rd = sqrt(9.8_r8*deltarho*thick1/(rhozero*fzero**2))
       lambda = ly/rd
  
  
       lambda = lambda + increment
       rd = ly/lambda
  
       deltarho = rd**2*fzero**2*rhozero/(9.8_r8*thick1)
       call modparam("deltarho",deltarho)
     case(2)
       thick1 = param("thick1",1000.0_r8,0)
       thick2 = param("thick2",4000.0_r8,0)
       dens1 = param("dens1",1022.0_r8,0)
       dens2 = param("dens2",1032.0_r8,0)
       fzero = param("fzero",1.0e-4_r8,0)
       rd = sqrt(9.8_r8*(dens2-dens1)*thick1*thick2/&
                     (dens1*fzero**2*(thick1+thick2)))
       rd = rd + increment
       dens2 = dens1+(rd**2*dens1*fzero**2*(thick1+thick2))/&
                                       (9.8_r8*thick1*thick2)
       call modparam("dens2",dens2)
     case default !bail out
       print *, "Wrong number of layers; aborting."
       stop
     end select
     endif
  
     if (trim(adjustl(inc_var)) == "delta") then
  ! figure out current value
       thick1 = param("thick1",1000.0_r8,0)
       thick2 = param("thick2",4000.0_r8,0)
       dens1 = param("dens1",1022.0_r8,0)
       dens2 = param("dens2",1032.0_r8,0)
       fzero = param("fzero",1.0e-4_r8,0)
  
       delta = thick1/thick2
       rd = sqrt(9.8_r8*(dens2-dens1)*thick1*thick2/&
                     (dens1*fzero**2*(thick1+thick2)))
  ! add increment
       delta = delta + increment
  
  ! figure out how to change relevant param(s).
  ! in this case we adjust thicknesses and densities (to hold ld constant)
       thick1 = delta*thick2
       dens2 = dens1+(rd**2*dens1*fzero**2*(thick1+thick2))/&
                                       (9.8_r8*thick1*thick2)
  
       call modparam("thick1",thick1)
       call modparam("dens2",dens2)
     endif
  
  ! write params
     call writeparams()
  enddo
  

  stop
end program pmc
