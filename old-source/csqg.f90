program stratqg
  use Types
  use Geometry
  use NewTrans
  use Control
  use Strat
  use Initialization
  use Diagnostics

  integer :: timestep,i

  call start_timing()
  call init_strat_run()
  print *,"initialization finished"
  call report_timing()
  
  call start_timing()
  

!print *, "get psi from q"
  do i=1,nmodes
    call solve_p_nondim(slab(i)%q,slab(i)%psi,slab(i)%s,alpha)
  enddo
  

!print *, "strat loop"
  do timestep = 1,nsteps
  !print *, "BONK -- beginloop"
    curstep = timestep
  
    if(modulo(timestep,nrg)==0) then
       slab(:)%get_energy = .true.
       slab(:)%get_enstrophy = .true.
    endif
  !print *, "BONK before step"
      call layer_step_strat_split(curstep,t)
  
  ! make adjustments to mean state
      do i=1,nmodes
         meanslab(i)%psi = meanslab(i)%psi + slab(i)%psi*dtold
         meanslab(i)%q = meanslab(i)%q + slab(i)%q*dtold
      enddo
      meant = meant + dtold
  
  !print *, "BONK - after step"
  
    if((modulo(timestep,nframe) == 0) .or. (timestep==2)) then
       print *, "writing files at step ",timestep
       call save_slab(slab,curstep)
       call save_mean(meanslab,curstep,meant)
  
       print "(a3,es15.8)","t=",t
       do i=1,nmodes
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
  
    if(modulo(timestep,nrg)==0) then
       ! write energy-------------------------------------------------------
       write(unit=efnum,fmt=*) t,slab(:)%tote
       write(unit=skenum,fmt=*) t,slab(:)%ske
       write(unit=ensnum,fmt=*) t,slab(:)%enstrophy
    endif
  !print *, "BONK -- endloop",timestep
  enddo
  call report_timing(timestep)
  

!print *, "cleanup"
  close(unit=skenum)
  close(unit=ensnum)
  !output animation
  !call dx_write_anim()
  close(unit=efnum)
  close(unit=ekenum)
  

  stop
end program stratqg
