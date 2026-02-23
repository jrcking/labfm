module ns_equations_cylinder
  use kind_parameters
  use common_parameter
  use common_2d
  use analytic_functions
  use nodes
  use derivatives
  use omp_lib
  implicit none
  
  private
  public :: solve_ns_equations_cylinder
  

  !! Solve Navier-Stokes equations in weakly-compressible form.
  !! Currently only suitable for periodic domains
  
  real(rkind),dimension(:),allocatable :: RHS_ro,RHS_u,RHS_v,RHS_Y
  real(rkind),dimension(:),allocatable :: vort
  
  real(rkind),parameter :: Ma=0.1d0
  real(rkind),parameter :: Ra = 1.0d2
  real(rkind),parameter :: Pr = 1.0d0
  real(rkind),parameter :: output_period=0.01d0
  real(rkind),parameter :: Da = 1.0d0
  
  real(rkind),parameter :: ooMa2 = one/Ma/Ma
  real(rkind),parameter :: momdiff_coeff = sqrt(Pr/Ra)
  real(rkind),parameter :: Ydiff_coeff = one/sqrt(Pr*Ra)

  !! Random coeffs for initial conditions
  integer(ikind),parameter :: ncoef=10
  real(rkind),dimension(ncoef) :: Ymag_coef,Yphase_coef
  
  real(rkind),dimension(:,:),allocatable :: utran
  integer(ikind),parameter :: move_every_nsteps=1
  real(rkind),parameter :: alpha=0.5d0  !! Narrowness parameter for mapping Y->Z    
  real(rkind),parameter :: dltah0 = three

contains
  subroutine solve_ns_equations_cylinder
     integer(ikind) :: i,j,k,n_out
     real(rkind) :: tmp,x,y
     real(rkind) :: l2_tmp,tot_vol,vol_local,meanu,meanv

     !! Simulation time     
     time_end = 200.0d0
     
     !! Initialise time, u and v
     n_out = 0  
     call set_initial_fields
     call rebuild_particle_maps_etc   
     allocate(vort(npfb))           
  
     call reapply_mirror_bcs     

     !!! Integrate from time to time_end
     do while (time.le.time_end)
        itime = itime + 1

        !! Outputs
        if(itime.eq.1.or.time.gt.n_out*output_period) then      !! If we want to output every 0.01 for making an animation
!        if(mod(itime,1).eq.0) then

          n_out = n_out + 1
          call output_uv(n_out)
          write(6,*) itime,time,sqrt(maxval(u(1:npfb)**two + v(1:npfb)**two)), &
                     maxval(h(1:npfb))/h0,minval(h(1:npfb))/h0

          !! Save something to file                       
          l2_tmp = zero
          tot_vol = zero
          meanu = zero
          meanv = zero
          do i=1,npfb
             vol_local = h(i)*h(i)/(hovdx*hovdx)  !! Local particle volume
             l2_tmp = l2_tmp + ro(i)*(u(i)*u(i) + v(i)*v(i))*vol_local
             tot_vol = tot_vol + vol_local
             meanu = meanu + u(i)*vol_local
             meanv = meanv + v(i)*v(i)*vol_local
          end do
          l2_tmp = half*l2_tmp/tot_vol
          !! Volume averaged K.E.
          meanu = meanu/tot_vol
          meanv = sqrt(meanv/tot_vol)
          !! Iteration time, dt, mean kinetic energy,mean velocity components
          write(212,*) itime,time,dt,l2_tmp,meanu,meanv,tot_vol
          flush(212)

          if(time.gt.4.0d0) then
             write(313,*) rp(1:npfb,1)
             write(314,*) rp(1:npfb,2)
             write(315,*) h(1:npfb)
             write(316,*) u(1:npfb)
             write(317,*) v(1:npfb)
             write(318,*) utran(1:npfb,1)
             write(319,*) utran(1:npfb,2)
             flush(313)
             flush(314)
             flush(315)
             flush(316)
             flush(317)
             flush(318)
             flush(319)                                                                 
           
          end if


        end if
        
        call set_tstep
         
        call step_rk3
        
                
      
     end do
   
     !! Output some velocity profiles
     call output_slices

     
     !! Output to screen and file

      
     stop
  end subroutine solve_ns_equations_cylinder
!! ------------------------------------------------------------------------------------------------
  subroutine update_utran    
     integer(ikind) :: i,j,k
     real(rkind) ::x,y,dYdx,dYdy,absgradY,htmp
     real(rkind) :: qkd_mag,rad,qq
     real(rkind),dimension(dims) :: dstmp,rij,gradkernel
     real(rkind),dimension(:,:),allocatable :: gradE,gradY
     real(rkind),parameter :: swarm_rate = dble(move_every_nsteps)*one
     real(rkind) :: max_vel,norm_gradY

     !! Set E=(4Y(1-Y))**(1/4)
     max_vel = zero
     do i=1,np

        E(i) = (4.0d0*Yspec(i)*(one-Yspec(i)))**alpha
#if cyl==1        
        !! Hard-code a small highly resolved clump...
        x=rp(i,1);y=rp(i,2);rad=sqrt(x*x+y*y)
        qq = half*(tanh((rad-cyl_radius)/(dltah0*h0)) +one)
        htmp = (4.0d0*qq*(one-qq))**alpha           
        E(i) = max(E(i),htmp)
#endif      
#if line==1
        x=rp(i,1);
        qq = half*(tanh(x/(dltah0*h0))+one)
        htmp = (4.0d0*qq*(one-qq))**alpha
        E(i) = max(E(i),htmp)  
#endif        
  
        max_vel = max(max_vel,u(i)*u(i) + v(i)*v(i))
     end do
     max_vel = sqrt(max_vel)
     
     !! Calculate the gradient of the test function...
     allocate(gradE(npfb,dims));gradE=zero
     call calc_gradient(E,gradE)
     
     !! Calculate the gradient of Y
     allocate(gradY(npfb,dims));gradY=zero
     call calc_gradient(Yspec,gradY)     
     
            
     !$omp parallel do private(x,y,dYdy,dYdx,absgradY,j,k &
     !$omp ,dstmp,qkd_mag,rij,rad,qq,gradkernel,htmp,norm_gradY)
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


        utran(i,1) = utran(i,1) + (minval(h(1:npfb))/dt)*swarm_rate*h0*dYdx
        utran(i,2) = utran(i,2) + (minval(h(1:npfb))/dt)*swarm_rate*h0*dYdy

        !! Normalised gradient of resolution
        !norm_gradY = h0*sqrt(gradY(i,1)**two + gradY(i,2)**two)        

        !! Dynamically evolve smoothing length to satisfy partition of unity
#if ale==1
        htmp = zero
        do k=1,ij_count(i)
           j=ij_link(i,k)
           rij=rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=rad/h(i)
           htmp = htmp + wab(qq)/(hovdx*hovdx)           
        end do
        qq = one - 0.01d0*(htmp - one)
        h(i) = min(max(0.1d0*h0,qq*h(i)),2.0d0*h0)  !! h relaxes towards PARTITION OF UNITY, with some upper/lower bounds
#endif   
        !! Add a small shifting term
        qkd_mag = (minval(h(1:npfb))/dt)*swarm_rate*0.2d0*(one+(4.0d0*Yspec(i)*(one-Yspec(i)))**alpha)
#if cyl==1
!LAG
!    qkd_mag = (minval(h(1:npfb))/dt)*swarm_rate*0.4d0
        
        x=rp(i,1);y=rp(i,2);rad=sqrt(x*x+y*y)
        qq = half*(tanh((rad-cyl_radius)/(dltah0*h0)) +one)
        htmp = max((4.0d0*qq*(one-qq))**alpha,(4.0d0*Yspec(i)*(one-Yspec(i)))**alpha)        
        qkd_mag = (minval(h(1:npfb))/dt)*swarm_rate*0.2d0*(one+htmp)                        
#endif
#if line==1
        x=rp(i,1);
        qq = half*(tanh(x/(dltah0*h0))+one)
        htmp = (4.0d0*qq*(one-qq))**alpha
        qkd_mag = (minval(h(1:npfb))/dt)*swarm_rate*0.2d0*(one+htmp)                                
#endif        
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
     
     deallocate(gradE,gradY)    
           
     return     
     
  end subroutine update_utran
!! ------------------------------------------------------------------------------------------------
  subroutine set_initial_fields
     integer(ikind) :: i,j
     real(rkind) :: x,y,ftn,y0,dlta,eta,rr
     real(rkind) :: ur,uth,th

     time = zero
     itime = 0
     
     dlta = 0.04d0         

     !! Initialise coefficients     
     do i=1,ncoef
        Ymag_coef(i) = rand()
        Yphase_coef(i) = two*pi*rand()
     end do     

     
    
     
     do i=1,npfb
        x = rp(i,1);y=rp(i,2)
        rr = sqrt(x*x+y*y)
     
        !! Multi-mode perturbation
        y0 = ymin + 0.4*(ymax-ymin)
        eta = zero
        do j=1,ncoef
           eta = eta + 0.005*sin( two*dble(j)*pi*x + Yphase_coef(j))!*Ymag_coef(j)
        end do
        ftn = -half*(tanh((y-y0+eta)/dlta)-one)        
    
        u(i) = zero
        v(i) = zero
        Yspec(i) = zero!ftn 
        ro(i) = one! - Ma*Ma*dlta*log(cosh((y-y0+eta)/dlta))
        
        th = acos(x/rr) !! theta coordinate
        ur = (one - (cyl_radius/rr)**two)*cos(th)
        uth = -(one + (cyl_radius/rr)**two)*sin(th)
        
        !! Potential flow past cylinder - just to speed up initial transients
!        u(i) = ur*cos(th) - uth*sin(th)
!        v(i) = ur*sin(th) + uth*cos(th)        
!        ro(i) = one + Ma*Ma*(cos(two*th)*(cyl_radius/rr)**two - half*(cyl_radius/rr)**4.0d0)
        

        !! Taylor Green vortex initial conditions
!        u(i) = -cos(2.0*pi*x)*sin(2.0*pi*y)!*cos(z)!*oosqrt2
!        v(i) = sin(2.0*pi*x)*cos(2.0*pi*y)!*cos(z)    !!c c                                         
!        ro(i) = one -Ma*Ma*(1.00d0/4.0d0)*(cos(two*two*pi*x)+cos(two*two*pi*y))!*(two+cos(two*z))  

     end do
     
     call reapply_mirror_bcs
     call adjust_for_symmetry_bcs

     
     allocate(utran(npfb,2))
     utran = zero
  
     return
  end subroutine set_initial_fields
!! ------------------------------------------------------------------------------------------------    
  subroutine adjust_for_symmetry_bcs
     integer(ikind) :: i
     real(rkind) :: dro


     return
  end subroutine adjust_for_symmetry_bcs
!! ------------------------------------------------------------------------------------------------
!! ------------------------------------------------------------------------------------------------
  subroutine set_tstep
     integer(ikind) :: i
     real(rkind) :: umax,umag,dt_visc,dt_acou,dt_spec,c_ac,dt_move,smin
          
     c_ac = one/Ma
  
     umax = zero
     do i=1,npfb
        umag = sqrt(u(i)*u(i)+v(i)*v(i))
        umax = max(umax,umag)        
     end do
     smin = minval(h(1:npfb))/hovdx
     
     dt_visc = 0.2*smin*smin/momdiff_coeff
     dt_acou = 0.8*smin/(umax + c_ac)
     dt_spec = 0.2*smin*smin/Ydiff_coeff
     
     dt = min(dt_visc,min(dt_acou,dt_spec))            
       
     return
  end subroutine set_tstep  
!! ------------------------------------------------------------------------------------------------
  subroutine rebuild_particle_maps_etc
     use nodes
     use neighbours
     use moments
#if ale==1
     integer(ikind) :: i,j
     
     if(mod(itime,move_every_nsteps).eq.0) then     

        !! Replace particles which have left domain
        !$omp parallel do
        do i=1,npfb
           if(xbcond.eq.1) then
              if(rp(i,1).le.xmin) rp(i,1) = rp(i,1) + (xmax-xmin)
              if(rp(i,1).gt.xmax) rp(i,1) = rp(i,1) - (xmax-xmin)
           else if(xbcond.eq.2) then
              if(rp(i,1).le.xmin) rp(i,1) = rp(i,1) + two*(xmin-rp(i,1));utran(i,1)=-utran(i,1);u(i)=-u(i)
              if(rp(i,1).gt.xmax) rp(i,1) = rp(i,1) - two*(rp(i,1)-xmax);utran(i,1)=-utran(i,1);u(i)=-u(i)
           end if
           if(ybcond.eq.1) then
              if(rp(i,2).le.ymin) rp(i,2) = rp(i,2) + (ymax-ymin)
              if(rp(i,2).gt.ymax) rp(i,2) = rp(i,2) - (ymax-ymin)
           else if(ybcond.eq.2) then
              if(rp(i,2).le.ymin) rp(i,2) = rp(i,2) + two*(ymin-rp(i,2));utran(i,1)=-utran(i,2);v(i)=-v(i)
              if(rp(i,2).gt.ymax) rp(i,2) = rp(i,2) - two*(rp(i,2)-ymax);utran(i,1)=-utran(i,2);v(i)=-v(i)           
           end if              
        end do
        !$omp end parallel do
   
        !! Rebuild mirrors
        call create_mirror_particles
     
        deallocate(ij_count,ij_link)
        call find_neighbours
          
        if(allocated(ij_w_grad)) deallocate(ij_w_grad,ij_w_lap)
        if(allocated(ij_w_hyp)) deallocate(ij_w_hyp,ij_w_hyp2)
        if(allocated(filter_coeff)) deallocate(filter_coeff)              
!        call calc_interparticle_weights      
        call calc_interparticle_weights_nofilt_o2
!        call calc_interparticle_weights_nofilt_sph        
        
     end if
#endif
     return
  end subroutine rebuild_particle_maps_etc
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
     
            

     !! Forcing::
     !! Re=1, f=40.0d0*momdiff_coeff
     !! Re=10, f= 41.5
     !! Re=100, f= 49 ish
     !! Re=1000, f= 84ish


     !$OMP PARALLEL DO 
     do i=1,npfb
        RHS_ro(i) = -ro(i)*(gradu(i,1)+gradv(i,2)) - u(i)*gradro(i,1) - v(i)*gradro(i,2) 
        RHS_u(i) = -u(i)*gradu(i,1) - v(i)*gradu(i,2) &
                   - (ooMa2)*gradro(i,1) + momdiff_coeff*lapu(i) + 16.0d0*momdiff_coeff
        RHS_v(i) = -u(i)*gradv(i,1) - v(i)*gradv(i,2) &
                   - (ooMa2)*gradro(i,2) + momdiff_coeff*lapv(i) !+ two*(Yspec(i)-half)
        RHS_Y(i) = zero!-u(i)*gradY(i,1) - v(i)*gradY(i,2) + Ydiff_coeff*lapY(i) + Da*4.0*Yspec(i)*(one-Yspec(i))
        
        !! Calculate vorticity...
        vort(i) = gradv(i,1)-gradu(i,2)
        w(i) = vort(i)
     end do
     !$OMP END PARALLEL DO
     
#if ale==1     
     if(mod(itime,move_every_nsteps).eq.0) then
        !$omp parallel do
        do i=1,npfb
           RHS_ro(i) = RHS_ro(i) + utran(i,1)*gradro(i,1) + utran(i,2)*gradro(i,2)
           RHS_u(i) = RHS_u(i) + utran(i,1)*gradu(i,1) + utran(i,2)*gradu(i,2)
           RHS_v(i) = RHS_v(i) + utran(i,1)*gradv(i,1) + utran(i,2)*gradv(i,2)
           RHS_Y(i) = zero!RHS_Y(i) + utran(i,1)*gradY(i,1) + utran(i,2)*gradY(i,2)
        end do
        !$omp end parallel do
     end if
#endif

     
     deallocate(lapu,lapv,lapY)
     deallocate(gradu,gradv,gradro,gradY)

     return
  end subroutine calc_all_rhs
!! ------------------------------------------------------------------------------------------------
  subroutine step_RK3
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
#if ale==1
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


     !! Store primary variables in register 1 (w-register)
     !$omp parallel do
     do i=1,npfb
        ro_reg1(i)=ro(i);u_reg1(i)=u(i);v_reg1(i)=v(i)
        Y_reg1(i) = Yspec(i)
#if ale==1
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
        
        
#if ale==1
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
           Yspec(i) = max(zero,min(Yspec(i),one))

           !! Store next S in register 1      
           ro_reg1(i) = ro_reg1(i) + RKb(k)*rhs_ro(i)
           u_reg1(i) = u_reg1(i) + RKb(k)*rhs_u(i)
           v_reg1(i) = v_reg1(i) + RKb(k)*rhs_v(i) 
           Y_reg1(i) = Y_reg1(i) + RKb(k)*rhs_Y(i)           
        end do
        !$omp end parallel do
              
        call rebuild_particle_maps_etc              
  
        call reapply_mirror_bcs
        call adjust_for_symmetry_bcs   
                
     end do
     
     !! Final substep: returns solution straight to ro,u,v,E (register 2)
     !! and doesn't update S
     
     !! Set the intermediate time
     time = time0 + RKc(4)      
     
     !! Build the right hand side
     call update_utran
     call calc_all_rhs  


#if ale==1
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
        Yspec(i) = max(zero,min(Yspec(i),one))
     end do
     !$omp end parallel do  

     !! Deallocation
     deallocate(u_reg1,v_reg1,ro_reg1,Y_reg1)
     deallocate(rhs_u,rhs_v,rhs_ro,rhs_Y) 
#if ale==1
     deallocate(rp_reg1)
#endif         
     
     call rebuild_particle_maps_etc   
     
     call correct_mass
         
     call reapply_mirror_bcs
     call adjust_for_symmetry_bcs     
              
     !! Filter the solution 
!     call calc_filtered_var(ro)
!     call calc_filtered_var(u)
!     call calc_filtered_var(v)
!!     call calc_filtered_var(Yspec)
!     call reapply_mirror_bcs
!     call adjust_for_symmetry_bcs           
  
     return
  end subroutine step_RK3
!! ------------------------------------------------------------------------------------------------ 
  subroutine correct_mass
     integer(ikind) :: i
     real(rkind) :: tm,tv,vol_local,dro
     !! Calculate average density deviation (from rho_char)
     tm = zero;tv = zero
     
     
     
     !$omp parallel do private(vol_local) reduction(+:tm,tv)
     do i=1,npfb
        vol_local = h(i)*h(i)/(hovdx*hovdx)  !! Local particle volume
        tm = tm + (ro(i)-one)*vol_local  !! N.B. there is no need to scale to make dimensional
        tv = tv + vol_local
     end do
     !$omp end parallel do         
     dro = tm/tv   
     
     !! Adjust density uniformly
     do i=1,npfb
        ro(i) = ro(i) - dro
     end do

     return
  end subroutine correct_mass
!! ------------------------------------------------------------------------------------------------
  subroutine output_slices
     !! Interpolate properties at 3 vertical slices
     integer(ikind),parameter :: nsy=100
     integer(ikind) i,j,k,jlow,jhigh
     real(rkind),dimension(nsy) :: Sy,uy,vy,vorty
     real(rkind) :: dy,Sx,dst,qq,waq

     !! Position of slice in X
     Sx = 0.375

     dy = one/dble(nsy)
     do i=1,nsy
        Sy(i) = -half + half*dy + dble(i-1)*dy
     end do
     
     uy=zero
     vy=zero
     vorty=zero
     
     !! N.B. larger h/s for interpolation
     
     do i=1,np  !! Loop all particles
        if(abs(rp(i,1)-Sx).gt.two*two*h(i)) cycle !! Skip particles not near transect
        jlow = floor((rp(i,2)-two*two*h(i)-half*dy+half)/dy + one)
        jhigh = floor((rp(i,2)+two*two*h(i)-half*dy+half)/dy + one) + 1 
        jlow = max(jlow,1)
        jhigh = min(jhigh,nsy) 

        !! Loop over nearby elements in slice and add contribution
        do j=jlow,jhigh
                
           dst = sqrt((rp(i,2)-Sy(j))**two + (rp(i,1)-Sx)**two)

           if(dst/(two*h(i)).le.two)then 
              qq=dst/(two*h(i))
              waq = wab(qq)/(two*two*hovdx*hovdx)
              uy(j) = uy(j) + waq*u(i)
              vy(j) = vy(j) + waq*v(i)              
           endif                   
        end do
     end do
     

     
     

     do i=1,nsy
        write(219,*) Sy(i),uy(i),vy(i)
     end do
     flush(219)
     
     return
  end subroutine output_slices
!! ------------------------------------------------------------------------------------------------   
end module ns_equations_cylinder
