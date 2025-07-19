module swarm
  use kind_parameters
  use common_parameter
  use common_2d
  use analytic_functions
  use nodes
  use derivatives
  use omp_lib
  use sphtools
  implicit none
  
  !! Define scalar field Y as some analytic function
  !! Let this move across the domain
  !! Currently only suitable for periodic domains
  
  real(rkind),parameter :: Ma=0.02d0
  real(rkind),parameter :: Re_local=1.0d4
  real(rkind),parameter :: At=0.4d0
  real(rkind),parameter :: Sc=1000.0d0 
  real(rkind),parameter :: output_period=0.02d0
  real(rkind),parameter :: swarm_rate = 1.0d0
   
  real(rkind) :: x0,y0,u0,v0
  real(rkind) :: l2_error
   
  integer(ikind),parameter :: ncoef=100
  real(rkind),dimension(ncoef) :: umag_coef,uphase_coef
  real(rkind),dimension(ncoef) :: vmag_coef,vphase_coef 
  real(rkind),dimension(ncoef) :: shape_coef_x,shape_coef_y 
  real(rkind),dimension(:),allocatable :: dYdx_exact

contains
  subroutine swarm_test
     integer(ikind) :: i,j,k,n_out
     real(rkind) :: tmp,x,y
     real(rkind) :: l2_tmp,l2v_tmp
     
!     !! Initial swarming
!     call swarm_on_initial_conditions
     
     !! Initialise time, u and v
     n_out = 0  
     time = zero
     time_end = 1.0d2
     itime = 0     
     x0=zero;y0=zero
     u0=zero;v0=zero

     do i=1,ncoef
        umag_coef(i) = rand()
        uphase_coef(i) = two*pi*rand()
        vmag_coef(i) = rand()
        vphase_coef(i) = two*pi*rand()
        shape_coef_x(i) = rand()
        shape_coef_y(i) = rand()
     end do

     !! Re-do neighbours and so on.
     call rebuild_particle_maps_etc

     allocate(dYdx_exact(npfb))
     call set_analytic_Yspec
     do i=1,npfb
        if(E(i).gt.half) then
           ro(i) = one
        else
           ro(i) = zero
        end if
     end do

     !!! Integrate from time to time_end
     do while (time.le.time_end)
        itime = itime + 1
       
        call set_tstep
      
        !! Move forward in time
        time = time + dt         
        call set_analytic_Yspec

        !! Update the transport velocity. This is the bit which will be the automata
        call update_utran                      

        !! Outputs
        if(time.gt.n_out*output_period) then      !! If we want to output every 0.01 for making an animation
!        if(mod(itime,10).eq.0) then

          n_out = n_out + 1
          call output_uv(n_out)
          write(6,*) itime,100.0*time/time_end,l2_error
        end if
        
        !! Move the nodes (i.e. swarm)
        call move_nodes
        
        !! Re-do neighbours and so on.
        call rebuild_particle_maps_etc
 
     end do
   
     !! Calculate the L2norms

     
     !! Output to screen and file

      
     stop
  end subroutine swarm_test
!! ------------------------------------------------------------------------------------------------  
  subroutine swarm_on_initial_conditions
     integer(ikind) :: i,j,k,n_out
     real(rkind) :: tmp,x,y
     real(rkind) :: l2_tmp,l2v_tmp
     logical :: keepgoing
     
     !! Initialise time, u and v
     n_out = 0  
     time = zero
     itime = 0     
     x0=zero;y0=zero
     u0=zero;v0=zero

     !! Re-do neighbours and so on.
     call rebuild_particle_maps_etc

     allocate(dYdx_exact(npfb))
     call set_analytic_Yspec

     !!! Swarm until some threshold is met
     keepgoing = .true.
     do while (keepgoing)
        itime = itime + 1
       
        call set_tstep  
        call set_analytic_Yspec

        !! Update the transport velocity. This is the bit which will be the automata
        call update_utran      
        
        !! What is the l2 of velocity?
        l2_tmp = zero
        do i=1,npfb
           l2_tmp = l2_tmp + u(i)*u(i) + v(i)*v(i)
        end do
        l2_tmp = sqrt(l2_tmp/npfb)   
        if(l2_tmp.le.0.2d0) keepgoing = .false.             

        !! Outputs
        if(itime.eq.1.or.mod(itime,100).eq.0) then

          n_out = n_out + 1
          write(6,*) itime,l2_error,l2_tmp
        end if
        
        !! Move the nodes (i.e. swarm)
        call move_nodes
        
        !! Re-do neighbours and so on.
        call rebuild_particle_maps_etc
 
     end do
     deallocate(dYdx_exact)
     write(6,*) "Initial shifting completed"
     write(6,*) "Moving to main shifting test"

      
     return
  end subroutine swarm_on_initial_conditions  
!! ------------------------------------------------------------------------------------------------
  subroutine update_utran
     integer(ikind) :: i,j,k
     real(rkind) ::x,y,dYdx,dYdy,absgradY,htmp
     real(rkind) :: qkd_mag,rad,qq
     real(rkind),dimension(dims) :: dstmp,rij,gradkernel
     real(rkind),dimension(:,:),allocatable :: gradY,gradE
    
     allocate(gradY(npfb,dims));gradY=zero
     call calc_gradient(Yspec,gradY)
     
     !! Calculate the gradient of the test function...
     allocate(gradE(npfb,dims));gradE=zero
     call calc_gradient(E,gradE)
     
     !! Calculate L2 norm of error
     qq = zero
     do i=1,npfb
        htmp = gradE(i,1) - dYdx_exact(i) !! Local error
        qq = qq + htmp*htmp  
!        ro(i) = htmp
     end do
     qq=sqrt(qq/npfb) !! L2
     l2_error = qq
     write(212,*) itime,time,l2_error,maxval(abs(ro(1:npfb)))
     flush(212)
     
     
    
     !$omp parallel do private(x,y,dYdy,dYdx,absgradY,j,k &
     !$omp ,dstmp,qkd_mag,rij,rad,qq,gradkernel,htmp)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2)
        
        !! Gradient
        dYdx = gradY(i,1)!*0.25d0*(4.0d0*E(i)*(1.0d0-E(i)))**(-3.0d0/4.0d0)
        dYdy = gradY(i,2)!*0.25d0*(4.0d0*E(i)*(1.0d0-E(i)))**(-3.0d0/4.0d0)


        !! Set velocity proportional to gradY
        u(i) = zero;v(i)=zero
!! Lagrangian
!u(i) = u0;v(i) = v0        
        u(i) = u(i) + (minval(h(1:npfb))/dt)*swarm_rate*h0*dYdx!*(abs(dYdx))**-0.5
        v(i) = v(i) + (minval(h(1:npfb))/dt)*swarm_rate*h0*dYdy!*(abs(dYdy))**-0.5

        !! Set smoothing length to approx 24!
        htmp = zero
        do k=1,ij_count(i)
           j=ij_link(i,k)
           rij=rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=rad/h(i)
           htmp = htmp + wab(qq)/(hovdx*hovdx)           
        end do
        qq = one - 0.001d0*(htmp - one)
        h(i) = min(max(0.1d0*h0,qq*h(i)),2.0d0*h0)  !! h relaxes towards PARTITION OF UNITY, with some upper/lower bounds
  
        !! Add a small shifting term
        qkd_mag = 1.0d-1*(minval(h(1:npfb))/dt)*swarm_rate*h(i)/h0!htmp/h0        


        dstmp = zero
        do k=1,ij_count(i)
           j=ij_link(i,k)
           rij=rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=rad/h(i)
           gradkernel = qkd_mag*((one-half*qq))*rij(:)/max(rad,epsilon(rad))

!           if(ro(i).eq.1.and.ro(j).eq.0) gradkernel=gradkernel*1.2d0

           if(qq.gt.two) gradkernel(:)=zero
           dstmp = dstmp + gradkernel(:)                      
        end do
    
        !! Add shifting to the transport velocity
        u(i) = u(i) + dstmp(1)
        v(i) = v(i) + dstmp(2)


        !! Limit the transport velocity
        htmp = sqrt(u(i)*u(i) + v(i)*v(i))
        if(htmp.gt.0.2d0*h(i)/dt) then
           htmp = htmp/(0.2d0*h(i)/dt)
           u(i) = u(i)/htmp
           v(i) = v(i)/htmp
        end if


     end do
     !$omp end parallel do 
   
     deallocate(gradY)
  
     return
  end subroutine update_utran
!! ------------------------------------------------------------------------------------------------
  subroutine move_nodes
     integer(ikind) :: i
    
     !$omp parallel do
     do i=1,npfb
        rp(i,1) = rp(i,1) + dt*u(i)
        rp(i,2) = rp(i,2) + dt*v(i)
     end do
     !$omp end parallel do 
  
     return
  end subroutine move_nodes
!! ------------------------------------------------------------------------------------------------
  subroutine set_analytic_Yspec
     integer(ikind) :: i,j
     real(rkind) :: x,y,rr,dlta
     real(rkind) :: rrxp,rrxm,rryp,rrym,rrxpyp,rrxmyp,rrxpym,rrxmym
     real(rkind) :: sh_x,sh_y
     real(rkind) :: AFscl !! Length-scale of analytic function     

     dlta = 0.04 !- 0.01*sin(13.24d0*time) !! 0.075

     AFscl = 0.25 !+ 0.05*sin(12.453d0*time)

     u0 = zero
     v0 = zero
     do i=1,ncoef
!        u0 = u0 + 0.1d0*umag_coef(i)*cos(two*pi*dble(i)*time/dble(ncoef) +uphase_coef(i))
!        v0 = v0 + 0.1d0*vmag_coef(i)*cos(two*pi*dble(i)*time/dble(ncoef) +vphase_coef(i))        
        
     end do   
     
     sh_x = one !+ 0.15*cos(4.98*time)
     sh_y = one !- 0.15*cos(4.98*time)
!     u0=zero;v0=zero
     
     x0 = x0 + u0*dt
     y0 = y0 + v0*dt

     if(x0.le.xmin) x0 = x0 + (xmax-xmin)
     if(x0.gt.xmax) x0 = x0 - (xmax-xmin)
     if(y0.le.ymin) y0 = y0 + (ymax-ymin)
     if(y0.gt.ymax) y0 = y0 - (ymax-ymin)

     
     do i=1,npfb
        x = rp(i,1);y=rp(i,2)
             
        rr    = sqrt(sh_x*(x-x0)**two + sh_y*(y-y0)**two)       
        rrxp  = sqrt(sh_x*(x-x0+one)**two + sh_y*(y-y0)**two)  
        rrxm  = sqrt(sh_x*(x-x0-one)**two + sh_y*(y-y0)**two)  
        rryp  = sqrt(sh_x*(x-x0)**two + sh_y*(y-y0+one)**two)  
        rrym  = sqrt(sh_x*(x-x0)**two + sh_y*(y-y0-one)**two)  
        rrxpyp= sqrt(sh_x*(x-x0+one)**two + sh_y*(y-y0+one)**two)  
        rrxpym= sqrt(sh_x*(x-x0+one)**two + sh_y*(y-y0-one)**two)  
        rrxmyp= sqrt(sh_x*(x-x0-one)**two + sh_y*(y-y0+one)**two)  
        rrxmym= sqrt(sh_x*(x-x0-one)**two + sh_y*(y-y0-one)**two)                                          

        if(.true.) then    
        Yspec(i) = -half*(tanh((rr-AFscl)/dlta)-tanh((rr+AFscl)/dlta)) &          !! tanh
                 - half*(tanh((rrxp-AFscl)/dlta)-tanh((rrxp+AFscl)/dlta)) &
                 - half*(tanh((rrxm-AFscl)/dlta)-tanh((rrxm+AFscl)/dlta)) &
                 - half*(tanh((rryp-AFscl)/dlta)-tanh((rryp+AFscl)/dlta)) &
                 - half*(tanh((rrym-AFscl)/dlta)-tanh((rrym+AFscl)/dlta)) &
                 - half*(tanh((rrxpyp-AFscl)/dlta)-tanh((rrxpyp+AFscl)/dlta)) &
                 - half*(tanh((rrxpym-AFscl)/dlta)-tanh((rrxpym+AFscl)/dlta)) &
                 - half*(tanh((rrxmyp-AFscl)/dlta)-tanh((rrxmyp+AFscl)/dlta)) &
                 - half*(tanh((rrxmym-AFscl)/dlta)-tanh((rrxmym+AFscl)/dlta))
                 
        !! Analytic dY/dx
        dYdx_exact(i) = sh_x*(x-x0)*( (cosh((rr+AFscl)/dlta))**-two - (cosh((rr-AFscl)/dlta))**-two )/(two*dlta*rr) &
              + sh_x*(x-x0+one)*( (cosh((rrxp+AFscl)/dlta))**-two - (cosh((rrxp-AFscl)/dlta))**-two )/(two*dlta*rrxp) &
              + sh_x*(x-x0-one)*( (cosh((rrxm+AFscl)/dlta))**-two - (cosh((rrxm-AFscl)/dlta))**-two )/(two*dlta*rrxm) &    
              + sh_x*(x-x0)*( (cosh((rryp+AFscl)/dlta))**-two - (cosh((rryp-AFscl)/dlta))**-two )/(two*dlta*rryp) &
              + sh_x*(x-x0)*( (cosh((rrym+AFscl)/dlta))**-two - (cosh((rrym-AFscl)/dlta))**-two )/(two*dlta*rrym) &
              + sh_x*(x-x0+one)*( (cosh((rrxpyp+AFscl)/dlta))**-two - (cosh((rrxpyp-AFscl)/dlta))**-two )/(two*dlta*rrxpyp) &
              + sh_x*(x-x0+one)*( (cosh((rrxpym+AFscl)/dlta))**-two - (cosh((rrxpym-AFscl)/dlta))**-two )/(two*dlta*rrxpym) &
              + sh_x*(x-x0-one)*( (cosh((rrxmyp+AFscl)/dlta))**-two - (cosh((rrxmyp-AFscl)/dlta))**-two )/(two*dlta*rrxmyp) &
              + sh_x*(x-x0-one)*( (cosh((rrxmym+AFscl)/dlta))**-two - (cosh((rrxmym-AFscl)/dlta))**-two )/(two*dlta*rrxmym) 
        !! Analytic dY/dy                                                                   
!        dYdx_exact(i) = (y-y0)*( (cosh((rr+AFscl)/dlta))**-two - (cosh((rr-AFscl)/dlta))**-two )/(two*dlta*rr) &
!              + (y-y0)*( (cosh((rrxp+AFscl)/dlta))**-two - (cosh((rrxp-AFscl)/dlta))**-two )/(two*dlta*rrxp) &
!              + (y-y0)*( (cosh((rrxm+AFscl)/dlta))**-two - (cosh((rrxm-AFscl)/dlta))**-two )/(two*dlta*rrxm) &    
!              + (y-y0+one)*( (cosh((rryp+AFscl)/dlta))**-two - (cosh((rryp-AFscl)/dlta))**-two )/(two*dlta*rryp) &
!              + (y-y0-one)*( (cosh((rrym+AFscl)/dlta))**-two - (cosh((rrym-AFscl)/dlta))**-two )/(two*dlta*rrym) &
!              + (y-y0+one)*( (cosh((rrxpyp+AFscl)/dlta))**-two - (cosh((rrxpyp-AFscl)/dlta))**-two )/(two*dlta*rrxpyp) &
!              + (y-y0-one)*( (cosh((rrxpym+AFscl)/dlta))**-two - (cosh((rrxpym-AFscl)/dlta))**-two )/(two*dlta*rrxpym) &
!              + (y-y0+one)*( (cosh((rrxmyp+AFscl)/dlta))**-two - (cosh((rrxmyp-AFscl)/dlta))**-two )/(two*dlta*rrxmyp) &
!              + (y-y0-one)*( (cosh((rrxmym+AFscl)/dlta))**-two - (cosh((rrxmym-AFscl)/dlta))**-two )/(two*dlta*rrxmym)   
        
        !! Line in domain
!        Yspec(i) = half*((tanh((x)/dlta))+one)
        !! Analytic dY/dx
!        dYdx_exact(i) = (half/dlta)*(cosh(x/dlta))**-two
        
        E(i) = Yspec(i)
        Yspec(i) = (4.0d0*Yspec(i)*(one-Yspec(i)))**0.25d0
                 
        else
        Yspec(i) = exp(-(rr/AFscl)**two) &                !! Gaussian
                 + exp(-(rrxp/AFscl)**two) &
                 + exp(-(rrxm/AFscl)**two) &
                 + exp(-(rryp/AFscl)**two) &
                 + exp(-(rrym/AFscl)**two) &
                 + exp(-(rrxpyp/AFscl)**two) &
                 + exp(-(rrxpym/AFscl)**two) &
                 + exp(-(rrxmyp/AFscl)**two) &
                 + exp(-(rrxmym/AFscl)**two)                                                    
        end if                 
     end do
     
     call reapply_mirror_bcs
       
     return
  end subroutine set_analytic_Yspec 
!! ------------------------------------------------------------------------------------------------
  subroutine set_tstep
     integer(ikind) :: i
     real(rkind) :: umax,umag,dt_visc,dt_acou,dt_spec,c_ac,dt_adv
     real(rkind) :: smin
     
     c_ac = one/Ma
  
     umax = zero
     do i=1,npfb
        umag = sqrt(u(i)*u(i)+v(i)*v(i))
        umax = max(umax,umag)        
     end do
     umax = sqrt(u0*u0+v0*v0)
     smin = minval(h(1:npfb))/hovdx
     
     dt_visc = 0.5*Re_local*smin*smin
     dt_acou = 0.5*smin/(umax + c_ac)
!     dt_adv = 0.05*dx/sqrt(u0*u0+v0*v0)     
     dt_spec = dt_visc*Sc
     
     dt = min(dt_visc,min(dt_acou,dt_spec))
!     write(6,*) dt_visc,dt_acou,dt_spec

       
     return
  end subroutine set_tstep  
!! ------------------------------------------------------------------------------------------------
  subroutine rebuild_particle_maps_etc
     use nodes
     use neighbours
     use moments
     integer(ikind) :: i,j
     
     !! Replace particles which have left domain
     !$omp parallel do
     do i=1,npfb
        if(xbcond.eq.1) then
           if(rp(i,1).le.xmin) rp(i,1) = rp(i,1) + (xmax-xmin)
           if(rp(i,1).gt.xmax) rp(i,1) = rp(i,1) - (xmax-xmin)
        else if(xbcond.eq.2) then
           if(rp(i,1).le.xmin) rp(i,1) = rp(i,1) + two*(xmin-rp(i,1))
           if(rp(i,1).gt.xmax) rp(i,1) = rp(i,1) - two*(rp(i,1)-xmax)
        end if
        if(ybcond.eq.1) then
           if(rp(i,2).le.ymin) rp(i,2) = rp(i,2) + (ymax-ymin)
           if(rp(i,2).gt.ymax) rp(i,2) = rp(i,2) - (ymax-ymin)
        else if(ybcond.eq.2) then
           if(rp(i,2).le.ymin) rp(i,2) = rp(i,2) + two*(ymin-rp(i,2))
           if(rp(i,2).gt.ymax) rp(i,2) = rp(i,2) - two*(rp(i,2)-ymax)           
        end if              
     end do
     !$omp end parallel do
   
     !! Rebuild mirrors
     call create_mirror_particles
  
     deallocate(ij_count,ij_link)
     call find_neighbours
          
     if(allocated(ij_w_grad)) deallocate(ij_w_grad)
     if(allocated(ij_w_lap)) deallocate(ij_w_lap,ij_w_hyp,ij_w_hyp2)
     if(allocated(filter_coeff)) deallocate(filter_coeff)              
!     call calc_interparticle_weights
     call calc_interparticle_weights_nofilt_o2
!     call filter_coefficients           
        

     return
  end subroutine rebuild_particle_maps_etc
!! ------------------------------------------------------------------------------------------------
end module swarm
