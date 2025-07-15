module ns_equations
  use kind_parameters
  use common_parameter
  use common_2d
  use analytic_functions
  use nodes
  use derivatives
  use omp_lib
  implicit none
  
#define ale

  !! Solve Navier-Stokes equations in weakly-compressible form.
  !! Currently only suitable for periodic domains
  
  real(rkind),dimension(:),allocatable :: RHS_ro,RHS_u,RHS_v,RHS_Y!,RHS_Z
  
  real(rkind),parameter :: Ma=0.05d0
  real(rkind),parameter :: Re_local=1.0d2
  real(rkind),parameter :: Ra = 1.0d5
  real(rkind),parameter :: Pr = 1.0d1
  real(rkind),parameter :: output_period=0.1d0
  real(rkind),parameter :: Da = 0.0d0
  
  real(rkind),parameter :: ooMa2 = one/Ma/Ma
  real(rkind),parameter :: momdiff_coeff = sqrt(Pr/Ra)
  real(rkind),parameter :: Ydiff_coeff = one/sqrt(Pr*Ra)
  
  real(rkind),dimension(:,:),allocatable :: utran
!  real(rkind),dimension(:),allocatable :: Zspec
  integer(ikind),parameter :: move_every_nsteps=100

contains
  subroutine solve_ns_equations
     integer(ikind) :: i,j,k,n_out
     real(rkind) :: tmp,x,y
     real(rkind) :: l2_tmp,l2v_tmp
     
     !! Initialise time, u and v
     n_out = 0  
     call set_initial_fields
     call rebuild_particle_maps_etc              
  
     call reapply_mirror_bcs     

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
  subroutine update_utran    
     integer(ikind) :: i,j,k
     real(rkind) ::x,y,dYdx,dYdy,absgradY,htmp
     real(rkind) :: qkd_mag,rad,qq
     real(rkind),dimension(dims) :: dstmp,rij,gradkernel
     real(rkind),dimension(:,:),allocatable :: gradE
     real(rkind),parameter :: swarm_rate = 1.0d0     

     !! Set E=(4Y(1-Y))**(1/4)
     do i=1,np
        E(i) = (4.0d0*Yspec(i)*(one-Yspec(i)))**0.25d0
     end do
     
     !! Calculate the gradient of the test function...
     allocate(gradE(npfb,dims));gradE=zero
     call calc_gradient(E,gradE)
            
     !$omp parallel do private(x,y,dYdy,dYdx,absgradY,j,k &
     !$omp ,dstmp,qkd_mag,rij,rad,qq,gradkernel,htmp)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2)
        
        !! Gradient
        dYdx = gradE(i,1)
        dYdy = gradE(i,2)

        !! Set velocity proportional to gradY
        utran(i,:) = zero

        !! Lagrangian
        utran(i,1) = u(i)
        utran(i,2) = v(i)
        
        utran(i,1) = utran(i,1) + (minval(h(1:npfb))/dt)*swarm_rate*h0*dYdx!*(abs(dYdx))**-0.5
        utran(i,2) = utran(i,2) + (minval(h(1:npfb))/dt)*swarm_rate*h0*dYdy!*(abs(dYdy))**-0.5

        !! Dynamically evolve smoothing length to satisfy partition of unity
        htmp = zero
        do k=1,ij_count(i)
           j=ij_link(i,k)
           rij=rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=rad/h(i)
           htmp = htmp + wab(qq)/(hovdx*hovdx)           
        end do
        qq = one - 0.001d0*(htmp - one)
        h(i) = min(max(0.1d0*h0,qq*h(i)),2.0d0*h0)  !! h relaxes towards PARTITION OF UNITY, with some upper/lower bounds
  
        !! Add a small shifting term
        qkd_mag = 1.5d-1*(minval(h(1:npfb))/dt)*swarm_rate*h(i)/h0!htmp/h0         !1.0d-1


        dstmp = zero
        do k=1,ij_count(i)
           j=ij_link(i,k)
           rij=rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=rad/h(i)
           gradkernel = qkd_mag*((one-half*qq))*rij(:)/max(rad,epsilon(rad))
           if(qq.gt.two) gradkernel(:)=zero
           dstmp = dstmp + gradkernel(:)                      
        end do
    
        !! Add shifting to the transport velocity
        utran(i,1) = utran(i,1) + dstmp(1)
        utran(i,2) = utran(i,2) + dstmp(2)


        !! Limit the transport velocity
        htmp = sqrt(dot_product(utran(i,:),utran(i,:)))
        if(htmp.gt.0.2d0*h(i)/dt) then
           htmp = htmp/(0.2d0*h(i)/dt)
           utran(i,1) = utran(i,1)/htmp
           utran(i,2) = utran(i,2)/htmp
        end if


     end do
     !$omp end parallel do 
   
     deallocate(gradE)
      
     return     
     
  end subroutine update_utran
!! ------------------------------------------------------------------------------------------------
  subroutine set_initial_fields
     integer(ikind) :: i,j
     real(rkind) :: x,y,ftn,y0,dlta

     time = zero
     time_end = 1.0d2
     itime = 0
     
     dlta = 0.04d0
    
!     allocate(Zspec(np))
     
     do i=1,npfb
        x = rp(i,1);y=rp(i,2)
     
        y0 = 0.25d0*(one+0.01*sin(2.0d0*pi*x))
        ftn = -half*(tanh((y-y0)/dlta)-tanh((y+y0)/dlta))
               
     
        u(i) = zero
        v(i) = zero
        Yspec(i) = ftn
!        Zspec(i) = atanh(two*Yspec(i)-one)
        ro(i) = one

     end do
     
     do j=npfb+1,np
        i=irelation(j)
        u(j) = u(i)
        v(j) = v(i)
        ro(j) = ro(i)
        Yspec(j) = Yspec(i)
!        Zspec(j) = Zspec(i)
     end do
     
     allocate(utran(npfb,2))
     utran = zero
  
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
     real(rkind),dimension(:),allocatable :: Y_reg1!,Z_reg1
     real(rkind),dimension(3) :: RKa
     real(rkind),dimension(4) :: RKb,RKc
#ifdef ale
     real(rkind),dimension(:,:),allocatable :: rp_reg1
     allocate(rp_reg1(npfb,2))
#endif     
         
     !! Set RKa,RKb with dt (avoids multiplying by dt on per-node basis)
     RKa(:) = dt*rk3_4s_2r_a(:)
     RKb(:) = dt*rk3_4s_2r_b(:)
     RKc(:) = dt*rk3_4s_2r_c(:)     

     allocate(u_reg1(npfb),v_reg1(npfb),ro_reg1(npfb),)
     allocate(rhs_u(npfb),rhs_v(npfb),rhs_ro(npfb))
     allocate(Y_reg1(npfb),rhs_Y(npfb))
!     allocate(Z_reg1(npfb),rhs_Z(npfb))     

     !! Store primary variables in register 1 (w-register)
     !$omp parallel do
     do i=1,npfb
        ro_reg1(i)=ro(i);u_reg1(i)=u(i);v_reg1(i)=v(i)
        Y_reg1(i) = Yspec(i)
!        Z_reg1(i) = Zspec(i)
#ifdef ale
        rp_reg1(i,:) = rp(i,:)
#endif        
     end do
     !$omp end parallel do             

     !! Temporary storage of time
     time0=time

     do k=1,3
        
        !! Set the intermediate time
        time = time0 + RKc(k)         
        
        !! Calculate the RHS
        call update_utran
        call calc_all_rhs
        
        
#ifdef ale
        if(mod(itime,move_every_nsteps).eq.0) then
           !$omp parallel do
           do i=1,npfb
              rp(i,:) = rp_reg1(i,:) + RKa(k)*utran(i,:)      
              rp_reg1(i,:) = rp_reg1(i,:) + RKb(k)*utran(i,:)
           end do
           !$omp end parallel do        
        end if
#endif     
        !! Set w_i and new u,v
        !$omp parallel do
        do i=1,npfb
           !! Store next U in register 2         
           ro(i) = ro_reg1(i) + RKa(k)*rhs_ro(i)
           u(i) = u_reg1(i) + RKa(k)*rhs_u(i)
           v(i) = v_reg1(i) + RKa(k)*rhs_v(i)
           Yspec(i) = Y_reg1(i) + RKa(k)*rhs_Y(i)
!           Zspec(i) = Z_reg1(i) + RKa(k)*rhs_Z(i)           

           !! Store next S in register 1      
           ro_reg1(i) = ro_reg1(i) + RKb(k)*rhs_ro(i)
           u_reg1(i) = u_reg1(i) + RKb(k)*rhs_u(i)
           v_reg1(i) = v_reg1(i) + RKb(k)*rhs_v(i) 
           Y_reg1(i) = Y_reg1(i) + RKb(k)*rhs_Y(i)
!           Z_reg1(i) = Z_reg1(i) + RKb(k)*rhs_Z(i)           
        end do
        !$omp end parallel do
              
        call rebuild_particle_maps_etc              
   
        call reapply_mirror_bcs
                
     end do
     
     !! Final substep: returns solution straight to ro,u,v,E (register 2)
     !! and doesn't update S
     
     !! Set the intermediate time
     time = time0 + RKc(4)      
     
     !! Build the right hand side
     call update_utran
     call calc_all_rhs  

#ifdef ale
        if(mod(itime,move_every_nsteps).eq.0) then
           !$omp parallel do
           do i=1,npfb
              rp(i,:) = rp_reg1(i,:) + RKb(k)*utran(i,:)      
           end do
           !$omp end parallel do        
        end if
#endif  
     
     !$omp parallel do
     do i=1,npfb
        !! Final values of prim vars
        ro(i) = ro_reg1(i) + RKb(4)*rhs_ro(i)
        u(i) = u_reg1(i) + RKb(4)*rhs_u(i)
        v(i) = v_reg1(i) + RKb(4)*rhs_v(i)
        Yspec(i) = Y_reg1(i) + RKb(4)*rhs_Y(i)
!        Zspec(i) = Z_reg1(i) + RKb(4)*rhs_Z(i)        
     end do
     !$omp end parallel do  

     !! Deallocation
     deallocate(u_reg1,v_reg1,ro_reg1,Y_reg1)!,Z_reg1)
     deallocate(rhs_u,rhs_v,rhs_ro,rhs_Y)!,rhs_Z) 
#ifdef ale
     deallocate(rp_reg1)
#endif         
     
     call rebuild_particle_maps_etc     
     !! Add in here some switch to ensure Y\in[0,1]
     do i=1,npfb
        Yspec(i) = max(zero,min(Yspec(i),one))
     end do
     !! Evolution of Z=atanh(2Y-1) is a bit unstable...

   
     call reapply_mirror_bcs
              
     !! Filter the solution 
!     call calc_filtered_var(ro)
!     call calc_filtered_var(u)
!     call calc_filtered_var(v)
!     call calc_filtered_var(Yspec)
!     call calc_filtered_var(Zspec)     
!     call reapply_mirror_bcs
     
      
  
     return
  end subroutine step
!! ------------------------------------------------------------------------------------------------
  subroutine set_tstep
     integer(ikind) :: i
     real(rkind) :: umax,umag,dt_visc,dt_acou,dt_spec,c_ac,dt_move
     
     c_ac = one/Ma
  
     umax = zero
     do i=1,npfb
        umag = sqrt(u(i)*u(i)+v(i)*v(i))
        umax = max(umax,umag)        
     end do
     
     dt_visc = 0.1*dx*dx/momdiff_coeff
     dt_acou = 0.5*dx/(umax + c_ac)
     dt_spec = 0.1*dx*dx/Ydiff_coeff

     
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
!        Zspec(i) = Zspec(j)        
     end do
     !$OMP END PARALLEL DO

     return
  end subroutine reapply_mirror_bcs
!! ------------------------------------------------------------------------------------------------
  subroutine rebuild_particle_maps_etc
     use nodes
     use neighbours
     use moments
#ifdef ale
     integer(ikind) :: i,j
     
     if(mod(itime,move_every_nsteps).eq.0) then     

        !! Replace particles which have left domain
        !$omp parallel do
        do i=1,npfb
           if(rp(i,1).le.xmin) rp(i,1) = rp(i,1) + (xmax-xmin)
           if(rp(i,1).gt.xmax) rp(i,1) = rp(i,1) - (xmax-xmin)
           if(rp(i,2).le.ymin) rp(i,2) = rp(i,2) + (ymax-ymin)
           if(rp(i,2).gt.ymax) rp(i,2) = rp(i,2) - (ymax-ymin)
        end do
        !$omp end parallel do
   
        !! Rebuild mirrors
        deallocate(irelation) 
        call create_mirror_particles
     
        deallocate(ij_count,ij_link)
        call find_neighbours
          
        if(allocated(ij_w_grad)) deallocate(ij_w_grad,ij_w_lap)
        if(allocated(ij_w_hyp)) deallocate(ij_w_hyp,ij_w_hyp2)
        if(allocated(filter_coeff)) deallocate(filter_coeff)              
        call calc_interparticle_weights_nofilt_o2
!        call calc_interparticle_weights_sph             
!        call filter_coefficients           
        
     end if

#endif
     return
  end subroutine rebuild_particle_maps_etc
!! ------------------------------------------------------------------------------------------------
  subroutine calc_all_rhs
     integer(ikind) :: i
     real(rkind),dimension(:),allocatable :: lapu,lapv,lapY
     real(rkind),dimension(:,:),allocatable :: gradu,gradv,gradY,gradro!,gradZ
     
     
     !! Space for derivatives
     allocate(lapu(npfb),lapv(npfb),lapY(npfb))
     allocate(gradu(npfb,dims),gradv(npfb,dims),gradY(npfb,dims),gradro(npfb,dims))
!     allocate(gradZ(npfb,dims))
     
     !! Evaluate derivatives
     call calc_laplacian(u,lapu)
     call calc_laplacian(v,lapv)
     call calc_laplacian(Yspec,lapY)          

     call calc_gradient(u,gradu)
     call calc_gradient(v,gradv)
     call calc_gradient(Yspec,gradY)
!     call calc_gradient(Zspec,gradZ)     
     call calc_gradient(ro,gradro)               

     !$OMP PARALLEL DO 
     do i=1,npfb
        RHS_ro(i) = -ro(i)*(gradu(i,1)+gradv(i,2)) - u(i)*gradro(i,1) - v(i)*gradro(i,2)
        RHS_u(i) = -u(i)*gradu(i,1) - v(i)*gradu(i,2) &
                   - (ooMa2)*gradro(i,1) + momdiff_coeff*lapu(i) !+ 4.0d0*sin(4.0d0*pi*rp(i,2))
        RHS_v(i) = -u(i)*gradv(i,1) - v(i)*gradv(i,2) &
                   - (ooMa2)*gradro(i,2) + momdiff_coeff*lapv(i) + two*(Yspec(i)-half)
        RHS_Y(i) = -u(i)*gradY(i,1) - v(i)*gradY(i,2) + Ydiff_coeff*lapY(i) + Da*4.0*Yspec(i)*(one-Yspec(i))
        
!        RHS_Z(i) = -u(i)*gradZ(i,1) - v(i)*gradZ(i,2) + half*(one/Yspec(i) + one/(one-Yspec(i)))* &
!                                                        (Ydiff_coeff*lapY(i) + Da*4.0*Yspec(i)*(one-Yspec(i)))
        
     end do
     !$OMP END PARALLEL DO
     
#ifdef ale     
     if(mod(itime,move_every_nsteps).eq.0) then
        !$omp parallel do
        do i=1,npfb
           RHS_ro(i) = RHS_ro(i) + utran(i,1)*gradro(i,1) + utran(i,2)*gradro(i,2)
           RHS_u(i) = RHS_u(i) + utran(i,1)*gradu(i,1) + utran(i,2)*gradu(i,2)
           RHS_v(i) = RHS_v(i) + utran(i,1)*gradv(i,1) + utran(i,2)*gradv(i,2)
           RHS_Y(i) = RHS_Y(i) + utran(i,1)*gradY(i,1) + utran(i,2)*gradY(i,2)
!           RHS_Z(i) = RHS_Z(i) + utran(i,1)*gradZ(i,1) + utran(i,2)*gradZ(i,2)           
        end do
        !$omp end parallel do
     end if
#endif
     
     deallocate(lapu,lapv,lapY)
     deallocate(gradu,gradv,gradro,gradY)!,gradZ)

     return
  end subroutine calc_all_rhs
!! ------------------------------------------------------------------------------------------------
end module ns_equations
