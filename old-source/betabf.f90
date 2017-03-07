program betabf
  use Types
  use Geometry
  use NewTrans
  use Control
  use Layers
  use Integrator
  use Initialization
  use Diagnostics

  integer :: timestep,i

  call start_timing()
  call init_layer_run()
  print *,"initialization finished"
  call report_timing()
  
  
  !call test_multi()
  !stop
  call start_timing()
  

  call solve_p_nondim(slab(1)%q,slab(1)%psi,0.0_r8,alpha)
  

  do timestep = 1,nsteps
  !print *, "BONK -- beginloop"
    curstep = timestep
    st = 0.0
    nst = 0
      tdatloop: do !while(st < tdat)
        nst = nst +1
        if (st >= tdat) exit tdatloop
  
         if(dt > tdat - st) then
  !         print *, "small:",dt,tdat-st
           dt = tdat - st
           dt3 = dt/3.0_r8
           dt4 = dt/4.0_r8
           dt23= 2.0_r8*dt/3.0_r8
           cnconst = 2*alpha**4/(dt*deltam**3)
           ndt = cnconst
         endif
  
        call layer_step_bf(st)
  
  ! update mean state
         do i=1,nlayers
            meanslab(i)%q = meanslab(i)%q + slab(i)%q*dtold
            meanslab(i)%psi = meanslab(i)%psi + slab(i)%psi*dtold
            meanslab(i)%k1 = slab(i)%psi
            meanslab(i)%k2 = slab(i)%q
            call fcxcy(meanslab(i)%k1)
            call fcxcy(meanslab(i)%k2)
            meanslab(i)%dx = meanslab(i)%dx + meanslab(i)%k1**2*dtold
            meanslab(i)%dy = meanslab(i)%dy + meanslab(i)%k2**2*dtold
         enddo
         meant = meant + dtold
  
        if (st >= tdat) exit tdatloop
      enddo tdatloop
  
      t = t + st
  !    print *,"nst",nst,st,tdat
  !print *, "BONK - after step"
  
    if(modulo(timestep,nframe) == 0) then
      print *, "writing files at time",t/tfac
      call save_slab(slab,curstep)
      call save_mean(meanslab,curstep,meant)
      
      print "(a3,es15.8)","t=",t
      do i=1,nlayers
         print "(a6,es15.8,a8,es15.8)", &
                 "vmax=",slab(i)%vmax,"energy=",slab(i)%tote
      enddo
      
      totalsteps = startsteps + curstep
      call modparam("totalsteps",totalsteps)
      call modparam("t", t)
      call modparam("dt", dt)
      call modparam("nslab", curstep)
      
      call writeparams()
      
    endif
  
      ! get dimensional time
      t_dimen = t/tfac
      
      ! compute the energy
      do i=1,nlayers
         slab(i)%k1 = slab(i)%psi
         ! generate del2psi
         slab(i)%ddx = slab(i)%psi
         slab(i)%ddy = slab(i)%psi
         call d_dx_cheb(slab(i)%ddx)
         call d_dx_cheb(slab(i)%ddx)
         call d_dy_cheb(slab(i)%ddy)
         call d_dy_cheb(slab(i)%ddy)
         slab(i)%del2 = alpha**2*slab(i)%ddx+slab(i)%ddy
         
      
         call fcxcy(slab(i)%k1)
         call fcxcy(slab(i)%del2)
      
         call compute_energy(slab(i)%k1,slab(i)%del2,slab(i)%tote,slab(i)%ave,&
                    slab(i)%ske,slab(i)%deltapsi)
      
         if(saved_mean) then
            smeanslab(i)%k1 = slab(i)%k1-smeanslab(i)%psi
            smeanslab(i)%k2 = slab(i)%del2-smeanslab(i)%del2
            call compute_energy(smeanslab(i)%k1,smeanslab(i)%k2,smeanslab(i)%tote,&
                     smeanslab(i)%ave,smeanslab(i)%ske,smeanslab(i)%deltapsi)
         endif
      enddo
      
      !figure out which energies to record
      select case(which_nrg)
      case(1)  ! single layer
      write(unit=efnum,fmt="(3es20.9)") t_dimen,slab(1)%tote,slab(1)%deltapsi
      if(saved_mean) then
          write(unit=ekenum,fmt="(2es20.9)") t_dimen,smeanslab(1)%tote
      endif
      
      case(2)  ! 2 layer, with everything
      slab(1)%dx = psibt/slab(1)%s
      slab(1)%ddy = psibt/slab(1)%s
      slab(1)%ddx = psibt/slab(1)%s
      call d_dx_cheb(slab(1)%ddx)
      call d_dx_cheb(slab(1)%ddx)
      call d_dy_cheb(slab(1)%ddy)
      call d_dy_cheb(slab(1)%ddy)
      slab(1)%del2 = alpha**2*slab(1)%ddx+slab(1)%ddy
      call fcxcy(slab(1)%dx)
      call fcxcy(slab(1)%del2)
      call compute_energy(slab(1)%dx, slab(1)%del2, bte, btave,btske,btdeltapsi)
      
      slab(2)%dx = psibc
      slab(2)%ddy = psibc
      slab(2)%ddx = psibc
      call d_dx_cheb(slab(2)%ddx)
      call d_dx_cheb(slab(2)%ddx)
      call d_dy_cheb(slab(2)%ddy)
      call d_dy_cheb(slab(2)%ddy)
      slab(2)%del2 = alpha**2*slab(2)%ddx+slab(2)%ddy
      call fcxcy(slab(2)%dx)
      call fcxcy(slab(2)%del2)
      call compute_energy(slab(2)%dx, slab(2)%del2, bce, bcave,bcske,bcdeltapsi)
      
      slab(1)%ddx = slab(1)%psi
      slab(2)%ddx = slab(2)%psi
      call fcxcy(slab(1)%ddx)
      call fcxcy(slab(2)%ddx)
      
      call compute_ape(slab(2)%ddx,slab(1)%ddx,slab(1)%s,gamma,&
                   ape,avape,sape,avsape)
      
      write(unit=efnum,fmt="(13es20.9)") t_dimen,slab(1)%tote,&
             slab(1)%deltapsi,slab(2)%tote,slab(2)%deltapsi,&
             bte,bce,slab(1)%ske,slab(2)%ske,btske,bcske,ape,sape
      if(saved_mean) then
                smeanslab(1)%k1 = slab(1)%k1-smeanslab(1)%dx
                smeanslab(1)%k2 = slab(1)%del2-smeanslab(1)%dy
                call compute_energy(smeanslab(1)%k1,smeanslab(1)%k2,bte,&
                         btave,btske,btdeltapsi)
          
                smeanslab(2)%k1 = slab(2)%k1-smeanslab(2)%dx
                smeanslab(2)%k2 = slab(2)%del2-smeanslab(2)%dy
                call compute_energy(smeanslab(2)%k1,smeanslab(2)%k2,bce,&
                         bcave,bcske,bcdeltapsi)
          
          write(unit=ekenum,fmt="(11es20.9)") t_dimen,smeanslab(1)%tote,&
                smeanslab(1)%deltapsi,smeanslab(2)%tote,smeanslab(2)%deltapsi,&
                bte,bce,btske,bcske,btdeltapsi,bcdeltapsi
      endif
      
      case(3)    ! 2 layer without bt/bc
      slab(1)%ddx = slab(1)%psi
      slab(2)%ddx = slab(2)%psi
      call fcxcy(slab(1)%ddx)
      call fcxcy(slab(2)%ddx)
      
      call compute_ape(slab(2)%ddx,slab(1)%ddx,slab(1)%s,gamma,&
                   ape,avape,sape,avsape)
      
      write(unit=efnum,fmt="(9es20.9)") t_dimen,slab(1)%tote,&
             slab(1)%deltapsi,slab(2)%tote,slab(2)%deltapsi,&
             slab(1)%ske,slab(2)%ske,ape,sape
      if(saved_mean) then
          write(unit=ekenum,fmt="(3es20.9)") t_dimen,smeanslab(1)%tote,smeanslab(2)%tote
      endif
      
      case(4)  ! 2 layer, without ape
      slab(1)%dx = psibt/slab(1)%s
      slab(1)%ddy = psibt/slab(1)%s
      slab(1)%ddx = psibt/slab(1)%s
      call d_dx_cheb(slab(1)%ddx)
      call d_dx_cheb(slab(1)%ddx)
      call d_dy_cheb(slab(1)%ddy)
      call d_dy_cheb(slab(1)%ddy)
      slab(1)%del2 = alpha**2*slab(1)%ddx+slab(1)%ddy
      call fcxcy(slab(1)%dx)
      call fcxcy(slab(1)%del2)
      call compute_energy(slab(1)%dx, slab(1)%del2, bte, btave,btske,btdeltapsi)
      
      slab(2)%dx = psibc
      slab(2)%ddy = psibc
      slab(2)%ddx = psibc
      call d_dx_cheb(slab(2)%ddx)
      call d_dx_cheb(slab(2)%ddx)
      call d_dy_cheb(slab(2)%ddy)
      call d_dy_cheb(slab(2)%ddy)
      slab(2)%del2 = alpha**2*slab(2)%ddx+slab(2)%ddy
      call fcxcy(slab(2)%dx)
      call fcxcy(slab(2)%del2)
      call compute_energy(slab(2)%dx, slab(2)%del2, bce, bcave,bcske,bcdeltapsi)
      
      write(unit=efnum,fmt="(11es20.9)") t_dimen,slab(1)%tote,&
             slab(1)%deltapsi,slab(2)%tote,slab(2)%deltapsi,&
             bte,bce,slab(1)%ske,slab(2)%ske,btske,bcske
      
      if(saved_mean) then
                smeanslab(1)%k1 = slab(1)%k1-smeanslab(1)%dx
                smeanslab(1)%k2 = slab(1)%del2-smeanslab(1)%dy
                call compute_energy(smeanslab(1)%k1,smeanslab(1)%k2,bte,&
                         btave,btske,btdeltapsi)
          
                smeanslab(2)%k1 = slab(2)%k1-smeanslab(2)%dx
                smeanslab(2)%k2 = slab(2)%del2-smeanslab(2)%dy
                call compute_energy(smeanslab(2)%k1,smeanslab(2)%k2,bce,&
                         bcave,bcske,bcdeltapsi)
          
          write(unit=ekenum,fmt="(11es20.9)") t_dimen,smeanslab(1)%tote,&
                smeanslab(1)%deltapsi,smeanslab(2)%tote,smeanslab(2)%deltapsi,&
                bte,bce,btske,bcske,btdeltapsi,bcdeltapsi
      endif
      
      case(5) ! 2 layer, without either ape or bt/bc
      write(unit=efnum,fmt="(7es20.9)") t_dimen,slab(1)%tote,&
             slab(1)%deltapsi,slab(2)%tote,slab(2)%deltapsi,&
             slab(1)%ske,slab(2)%ske
      
      if(saved_mean) then
          write(unit=ekenum,fmt="(3es20.9)") t_dimen,smeanslab(1)%tote,smeanslab(2)%tote
      endif
      
      case(6)  ! 2 layer, with everything plus psimin/max
      slab(1)%dx = psibt/slab(1)%s
      slab(1)%ddy = psibt/slab(1)%s
      slab(1)%ddx = psibt/slab(1)%s
      call d_dx_cheb(slab(1)%ddx)
      call d_dx_cheb(slab(1)%ddx)
      call d_dy_cheb(slab(1)%ddy)
      call d_dy_cheb(slab(1)%ddy)
      slab(1)%del2 = alpha**2*slab(1)%ddx+slab(1)%ddy
      call fcxcy(slab(1)%dx)
      call fcxcy(slab(1)%del2)
      call compute_energy(slab(1)%dx, slab(1)%del2, bte, btave,btske,btdeltapsi)
      
      slab(2)%dx = psibc
      slab(2)%ddy = psibc
      slab(2)%ddx = psibc
      call d_dx_cheb(slab(2)%ddx)
      call d_dx_cheb(slab(2)%ddx)
      call d_dy_cheb(slab(2)%ddy)
      call d_dy_cheb(slab(2)%ddy)
      slab(2)%del2 = alpha**2*slab(2)%ddx+slab(2)%ddy
      call fcxcy(slab(2)%dx)
      call fcxcy(slab(2)%del2)
      call compute_energy(slab(2)%dx, slab(2)%del2, bce, bcave,bcske,bcdeltapsi)
      
      slab(1)%ddx = slab(1)%psi
      slab(2)%ddx = slab(2)%psi
      call fcxcy(slab(1)%ddx)
      call fcxcy(slab(2)%ddx)
      
      call compute_ape(slab(2)%ddx,slab(1)%ddx,slab(1)%s,gamma,&
                   ape,avape,sape,avsape)
      
      slab(1)%dx = psibt
      slab(1)%dy = psibc
      call fcxcy(slab(1)%dx)
      call fcxcy(slab(1)%dy)
      btminmax = minval(slab(1)%dx)+maxval(slab(1)%dx)
      bcminmax = minval(slab(1)%dy)+maxval(slab(1)%dy)
      write(unit=efnum,fmt="(17es20.9)") t_dimen,slab(1)%tote,&
             slab(1)%deltapsi,slab(2)%tote,slab(2)%deltapsi,&
             bte,bce,slab(1)%ske,slab(2)%ske,btske,bcske,ape,sape,&
             minmax1,minmax2,btminmax,bcminmax
      if(saved_mean) then
                smeanslab(1)%k1 = slab(1)%k1-smeanslab(1)%dx
                smeanslab(1)%k2 = slab(1)%del2-smeanslab(1)%dy
                call compute_energy(smeanslab(1)%k1,smeanslab(1)%k2,bte,&
                         btave,btske,btdeltapsi)
          
                smeanslab(2)%k1 = slab(2)%k1-smeanslab(2)%dx
                smeanslab(2)%k2 = slab(2)%del2-smeanslab(2)%dy
                call compute_energy(smeanslab(2)%k1,smeanslab(2)%k2,bce,&
                         bcave,bcske,bcdeltapsi)
          
          write(unit=ekenum,fmt="(11es20.9)") t_dimen,smeanslab(1)%tote,&
                smeanslab(1)%deltapsi,smeanslab(2)%tote,smeanslab(2)%deltapsi,&
                bte,bce,btske,bcske,btdeltapsi,bcdeltapsi
      endif
      
      end select
      
  !print *, "BONK -- endloop",timestep
  enddo
  call report_timing(timestep)
  

  !output animation
  !call dx_write_anim()
  close(unit=efnum)
  close(unit=ekenum)
  

  stop
end program betabf
