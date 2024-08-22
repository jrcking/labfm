module ns_equations
  use kind_parameters
  use common_parameter
  use common_2d
  use analytic_functions
  use nodes
  use derivatives
  use omp_lib
  implicit none

  !! Solve Navier-Stokes equations in weakly-compressible form.
  !! Currently only suitable for periodic domains
  
  real(rkind),dimension(:),allocatable :: RHS_ro,RHS_u,RHS_v,RHS_Y
  
  real(rkind),parameter :: Ma=0.1d0
  real(rkind),parameter :: Re_local=1.0d4
  real(rkind),parameter :: At=0.2d0
  real(rkind),parameter :: Sc=1000.0d0
  real(rkind),parameter :: output_period=0.1d0
  
  real(rkind),parameter :: ooMa2 = one/Ma/Ma
  real(rkind),parameter :: ooRe = one/Re_local
  real(rkind),parameter :: ooReSc = ooRe/Sc

contains
  subroutine solve_ns_equations
     integer(ikind) :: i,j,k,n_out
     real(rkind) :: tmp,x,y
     real(rkind) :: l2_tmp,l2v_tmp
     
     !! Initialise time, u and v
     n_out = 0  
     call set_initial_fields

     !!! Integrate from time to time_end
     do while (time.le.time_end)
        itime = itime + 1

        !! Outputs
        if(time.gt.n_out*output_period) then      !! If we want to output every 0.01 for making an animation
!        if(mod(itime,10).eq.0) then

          n_out = n_out + 1
          call output_uv(n_out)
          write(6,*) itime,100.0*time/time_end,maxval(abs(u(1:npfb))),maxval(abs(v(1:npfb)))
        end if
        
        call set_tstep
         
        call step
!        write(6,*) k,itime,100.0*time/time_end !! indicate progress to screen  

                      
 
     end do
   
     !! Calculate the L2norms

     
     !! Output to screen and file

      
     stop
  end subroutine solve_ns_equations
!! ------------------------------------------------------------------------------------------------
  subroutine set_initial_fields
     integer(ikind) :: i,j
     real(rkind) :: x,y

     time = zero
     time_end = 1.0d2
     itime = 0
    
     
     do i=1,npfb
        x = rp(i,1);y=rp(i,2)
     
        if(four*abs(y).le.one + 0.001*sin(4.0d0*pi*x)) then
           u(i) = half
           v(i) = zero
           Yspec(i) = one
           ro(i) = one
        else
           u(i) = -half
           v(i) = zero
           Yspec(i) = zero
           ro(i) = (one-At)/(one+At)
        end if
     end do
     
     do j=npfb+1,np
        i=irelation(j)
        u(j) = u(i)
        v(j) = v(i)
        ro(j) = ro(i)
        Yspec(j) = Yspec(i)
     end do
  
     return
  end subroutine set_initial_fields  
!! ------------------------------------------------------------------------------------------------
  subroutine step
     !! 3rd order 4step 2 register Runge Kutta
     !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
     !! 'U' and 'S' in comments relate to p31 of Senga2 user guide
     !! Implemented over three registers, because speed is more
     !! important than memory at present.    
     
     !! Register 1 is ro_reg1,rou_reg1,rov_reg1 etc
     !! Register 2 is ro,u,v etc
     !! Register 3 is rhs_ro, rhs_rou, rhs_rov etc
     integer(ikind) :: i,k
     real(rkind) :: time0
     real(rkind),dimension(:),allocatable :: u_reg1,v_reg1,ro_reg1
     real(rkind),dimension(:),allocatable :: Y_reg1
     real(rkind),dimension(3) :: RKa
     real(rkind),dimension(4) :: RKb,RKc
         
     !! Set RKa,RKb with dt (avoids multiplying by dt on per-node basis)
     RKa(:) = dt*rk3_4s_2r_a(:)
     RKb(:) = dt*rk3_4s_2r_b(:)
     RKc(:) = dt*rk3_4s_2r_c(:)     

     allocate(u_reg1(npfb),v_reg1(npfb),ro_reg1(npfb),)
     allocate(rhs_u(npfb),rhs_v(npfb),rhs_ro(npfb))
     allocate(Y_reg1(npfb),rhs_Y(npfb))

     !! Store primary variables in register 1 (w-register)
     !$omp parallel do
     do i=1,npfb
        ro_reg1(i)=ro(i);u_reg1(i)=u(i);v_reg1(i)=v(i)
        Y_reg1(i) = Yspec(i)
     end do
     !$omp end parallel do             

     !! Temporary storage of time
     time0=time

     do k=1,3
        
        !! Set the intermediate time
        time = time0 + RKc(k)         
        
        !! Calculate the RHS
        call calc_all_rhs
     
        !! Set w_i and new u,v
        !$omp parallel do
        do i=1,npfb
           !! Store next U in register 2
           ro(i) = ro_reg1(i) + RKa(k)*rhs_ro(i)
           u(i) = u_reg1(i) + RKa(k)*rhs_u(i)
           v(i) = v_reg1(i) + RKa(k)*rhs_v(i)
           Yspec(i) = Y_reg1(i) + RKa(k)*rhs_Y(i)

           !! Store next S in register 1
           ro_reg1(i) = ro_reg1(i) + RKb(k)*rhs_ro(i)
           u_reg1(i) = u_reg1(i) + RKb(k)*rhs_u(i)
           v_reg1(i) = v_reg1(i) + RKb(k)*rhs_v(i) 
           Y_reg1(i) = Y_reg1(i) + RKb(k)*rhs_Y(i)
        end do
        !$omp end parallel do
              
   
        call reapply_mirror_bcs
                
     end do
     
     !! Final substep: returns solution straight to ro,u,v,E (register 2)
     !! and doesn't update S
     
     !! Set the intermediate time
     time = time0 + RKc(4)      
     
     !! Build the right hand side
     call calc_all_rhs  
     
     !$omp parallel do
     do i=1,npfb
        !! Final values of prim vars
        ro(i) = ro_reg1(i) + RKb(4)*rhs_ro(i)
        u(i) = u_reg1(i) + RKb(4)*rhs_u(i)
        v(i) = v_reg1(i) + RKb(4)*rhs_v(i)
        Yspec(i) = Y_reg1(i) + RKb(4)*rhs_Y(i)
     end do
     !$omp end parallel do  

     !! Deallocation
     deallocate(u_reg1,v_reg1,ro_reg1,Y_reg1)
     deallocate(rhs_u,rhs_v,rhs_ro,rhs_Y)     
     
   
     call reapply_mirror_bcs
              
     !! Filter the solution 
     call calc_filtered_var(ro)
     call calc_filtered_var(u)
     call calc_filtered_var(v)
     call calc_filtered_var(Yspec)

     call reapply_mirror_bcs
      
  
     return
  end subroutine step
!! ------------------------------------------------------------------------------------------------
  subroutine set_tstep
     integer(ikind) :: i
     real(rkind) :: umax,umag,dt_visc,dt_acou,dt_spec,c_ac
     
     c_ac = one/Ma
  
     umax = zero
     do i=1,npfb
        umag = sqrt(u(i)*u(i)+v(i)*v(i))
        umax = max(umax,umag)        
     end do
     
     dt_visc = 0.1*Re*dx*dx
     dt_acou = 0.5*dx/(umax + c_ac)
     dt_spec = 0.1*Re*Sc*dx*dx
     
     dt = min(dt_visc,min(dt_acou,dt_spec))

       
     return
  end subroutine set_tstep  
!! ------------------------------------------------------------------------------------------------
  subroutine reapply_mirror_bcs
     integer(ikind) :: i,j
     !! Update phi in the boundary particles

        !$OMP PARALLEL DO PRIVATE(j)
        do i=npfb+1,np
           j = irelation(i)
           u(i) = u(j)
           v(i) = v(j) 
           ro(i) = ro(j)
           Yspec(i) = Yspec(j)
        end do
        !$OMP END PARALLEL DO

     return
  end subroutine reapply_mirror_bcs
!! ------------------------------------------------------------------------------------------------
  subroutine calc_all_rhs
     integer(ikind) :: i
     real(rkind),dimension(:),allocatable :: lapu,lapv,lapY
     real(rkind),dimension(:,:),allocatable :: gradu,gradv,gradY,gradro
     
     !! Space for derivatives
     allocate(lapu(npfb),lapv(npfb),lapY(npfb))
     allocate(gradu(npfb,dims),gradv(npfb,dims),gradY(npfb,dims),gradro(npfb,dims))
     
     !! Evaluate derivatives
     call calc_laplacian(u,lapu)
     call calc_laplacian(v,lapv)
     call calc_laplacian(Yspec,lapY)          

     call calc_gradient(u,gradu)
     call calc_gradient(v,gradv)
     call calc_gradient(Yspec,gradY)
     call calc_gradient(ro,gradro)               

     !$OMP PARALLEL DO 
     do i=1,npfb
        RHS_ro(i) = -ro(i)*(gradu(i,1)+gradv(i,2)) - u(i)*gradro(i,1) - v(i)*gradro(i,2)
        RHS_u(i) = -u(i)*gradu(i,1) - v(i)*gradu(i,2) &
                   - (ooMa2)*gradro(i,1) + 2.0d0*ooMa2*At*gradY(i,1)/(one+At) &
                   + ooRe*lapu(i) !+ 4.0d0*sin(4.0d0*pi*rp(i,2))
        RHS_v(i) = -u(i)*gradv(i,1) - v(i)*gradv(i,2) &
                   - (ooMa2)*gradro(i,2) + 2.0d0*ooMa2*At*gradY(i,2)/(one+At) &
                   + ooRe*lapv(i)
        RHS_Y(i) = -u(i)*gradY(i,1) - v(i)*gradY(i,2) + ooReSc*lapY(i)
        
        
     end do
     !$OMP END PARALLEL DO
     deallocate(lapu,lapv,lapY)
     deallocate(gradu,gradv,gradro,gradY)

     return
  end subroutine calc_all_rhs
!! ------------------------------------------------------------------------------------------------
end module ns_equations
