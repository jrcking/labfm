module burgers_equation
  use kind_parameters
  use common_parameter
  use common_2d
  use analytic_functions
  use nodes
  use derivatives
  use omp_lib
  implicit none

  !! Solve a Burgers equation. Note, because it's a bit messy, we won't use the generic time integration 
  !! utilities here, but specific burgers ones...
  !! if irelation is allocated, then we're doing some periodic burgers (simple), 
  !! if its unallocated, we're doing one with analytic BCs (given by ftype 3 in analytic functions...

contains
  subroutine solve_burgers_equation
     integer(ikind) :: i,j,k,n_out
     real(rkind) :: tmp,x,y
     real(rkind) :: l2_tmp,l2v_tmp

     !! Set the equation type indicator
     Re = 100.0        !! Set the Reynolds number
     
     !! Initialise time, u and v
     time = 0.0d0
     dt = min(0.1*hmin/1.0,0.05*hmin*hmin*Re)  !! Constraint is most stringent of advection or diffusion (taking max(U)=1 )
     time_end = 1.0  !! Burgers eqn is already non dimensional..
     itime = 0
     n_out = 0

     !! Set initial conditions
     if(allocated(irelation))then         !! Periodic Burgers
        !$OMP PARALLEL DO PRIVATE(x,y)
        do i=1,np
           x = rp(i,1);y=rp(i,2)
           u(i) = sin(2.0*pi*x/lambda)!ftn(x,y)
           v(i) = 0.0d0!-cos(2.0*pi*x/lambda)*cos(2.0*pi*y/lambda)   !! set v(i)=0 and modify u(i) if want 1D sawtooth    

           !! Pencil code test
!           u(i)=0.0d0-tanh((x-0.5d0)*Re/2.0d0)
!           v(i)=0.0d0               
        end do
        !$OMP END PARALLEL DO
     else                                  !! Aperiodic Burgers - shock development is balanced by viscous dissipation
        !$OMP PARALLEL DO PRIVATE(x,y)      
        do i=1,np                         !! so constant thickness smoothed shock travels diagonally...
           x=rp(i,1);y=rp(i,2)
           u(i) = 0.75 - ftn(x,y)
           v(i) = 0.75 + ftn(x,y) 
        end do
        !$OMP END PARALLEL DO
     end if

     !!! Integrate from time to time_end
     do while (time.le.time_end)
        itime = itime + 1
        if(time.gt.n_out*0.01) then      !! If we want to output every 0.01 for making an animation

          call tanhwiggle

          n_out = n_out + 1
          call output_uv(n_out)
          write(6,*) itime,100.0*time/time_end,maxval(u(1:npfb))
        end if
         
        call integrate_dt_RK4_burgers
!        write(6,*) k,itime,100.0*time/time_end !! indicate progress to screen  

        !! Filter the solution?
!        call calc_filtered_var(u)
!        call calc_filtered_var(v)                        
 
     end do

     if(allocated(irelation)) then
        call sawtooth   !! If sawtooth test case, call this to check convergence against anlytical solution
!call tanhwiggle        
     end if
   
     !! Calculate the L2norms
     l2_tmp = 0.0d0;l2v_tmp=0.0d0     
     !$OMP PARALLEL DO REDUCTION(+:l2_tmp,l2v_tmp) PRIVATE(tmp,x,y)
     do i=1,npfb
        x = rp(i,1);y = rp(i,2)
        tmp = w(i) 
!        tmp = u(i) - 0.75d0 + ftn(x,y) !! if irelation not allocated
        l2_tmp = l2_tmp + tmp*tmp
        tmp = w(i) 
!        tmp = v(i) - 0.75 - ftn(x,y)   !! if irelation not allocated
        l2v_tmp = l2v_tmp + tmp*tmp
     end do
     !$OMP END PARALLEL DO
     !! sqrt norms
     l2_tmp = sqrt(l2_tmp/npfb);l2v_tmp = sqrt(l2v_tmp/npfb)
     
     !! Output to screen and file
     write(6,*) npfb,itime,time,l2_tmp,l2v_tmp
     write(1,*) h0/lambda,l2_tmp,l2v_tmp
     flush(1)

     return
  end subroutine solve_burgers_equation
!! ------------------------------------------------------------------------------------------------
  subroutine integrate_dt_RK4_burgers
     integer(ikind) :: i,k
     real(rkind) :: time0
     real(rkind),dimension(:),allocatable :: u0,v0
     real(rkind),dimension(:,:),allocatable :: RKu,RKv
     real(rkind),dimension(4) :: RK_coeffs     

     allocate(u0(npfb),v0(npfb));u0=0.0d0;v0=0.0d0
     allocate(RKu(npfb,4),RKv(npfb,4));RKu=0.0d0;RKv=0.0d0
     RK_coeffs = (/ 1.0d0,0.5d0,0.5d0,1.0d0 /)
        
     !! Temporary storage of u and time
     time0=time
     !$OMP PARALLEL DO
     do i=1,npfb
        u0(i) = u(i)
        v0(i) = v(i)
     end do
     !$OMP END PARALLEL DO
       
     !! RK step 1
     call calc_rhs_burgers(1,RKu(:,1));call calc_rhs_burgers(2,RKv(:,1))

     !! RK steps 2,3,4
     do k=2,4
        !! Set the time for this substep
        time = time0 + RK_coeffs(k)*dt
           
        !! Find the temporary u for next substep
        !$OMP PARALLEL DO
        do i=1,npfb
           u(i) = u0(i) + RK_coeffs(k)*RKu(i,k-1)*dt
           v(i) = v0(i) + RK_coeffs(k)*RKv(i,k-1)*dt
        end do
        !$OMP END PARALLEL DO
        
        !! Apply BCs to intermediate vars
        call reapply_bcs_burgers
   
        !! calc the RHS and store
        call calc_rhs_burgers(1,RKu(:,k));call calc_rhs_burgers(2,RKv(:,k))
     end do

     !! Final weighted sum of substeps
     time = time0 + dt
     !$OMP PARALLEL DO
     do i=1,npfb
        u(i) = u0(i) + dt*(1.0/6.0)*(RKu(i,1) + 2.0*RKu(i,2) + 2.0*RKu(i,3) + RKu(i,4))
        v(i) = v0(i) + dt*(1.0/6.0)*(RKv(i,1) + 2.0*RKv(i,2) + 2.0*RKv(i,3) + RKv(i,4))
     end do
     !$OMP END PARALLEL DO


     call reapply_bcs_burgers      

     deallocate(RKu,RKv,u0,v0)
     return
  end subroutine integrate_dt_RK4_burgers
!! ------------------------------------------------------------------------------------------------
  subroutine reapply_bcs_burgers
     integer(ikind) :: i,j
     !! Update phi in the boundary particles

     if(allocated(irelation))then             !! Periodic Burgers
        !$OMP PARALLEL DO PRIVATE(j)
        do i=npfb+1,np
           j = irelation(i)
           u(i) = u(j)
           v(i) = v(j) 

!           u(i) = -tanh(((rp(i,1)-0.5d0))*Re/2.0d0)
!           v(i) = 0.0d0           
        end do
        !$OMP END PARALLEL DO
     else                                     !! Aperiodic burgers 
        !$OMP PARALLEL DO                     
        do i=npfb+1,np
           u(i) = 0.75 - ftn(rp(i,1),rp(i,2))       
           v(i) = 0.75 + ftn(rp(i,1),rp(i,2))          
        end do
        !$OMP END PARALLEL DO
     end if
     return
  end subroutine reapply_bcs_burgers
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_burgers(i_component,rhs)
     integer(ikind),intent(in) :: i_component   !! sets the component we're calculating now (i.e. u or v for vector eqns)
     real(rkind),dimension(:),intent(out) :: rhs
     integer(ikind) :: i
     real(rkind),dimension(dims) :: tmp_vec
     real(rkind) :: tmp_scal

     allocate(lap(npfb));lap=0.0d0
     allocate(grad(npfb,dims));grad=0.0d0
     if(i_component.eq.1)then  !! u component
        i_uv_burgers = 1    !! this flag tells d2fdx2 and d2fdy2 to change signs between f=u and f=v
        call calc_laplacian(u,lap)
        call calc_gradient(u,grad)
     else if(i_component.eq.2)then  !! v component
        i_uv_burgers = 2
        call calc_laplacian(v,lap)
        call calc_gradient(v,grad)
     end if
     !$OMP PARALLEL DO PRIVATE(tmp_vec,tmp_scal)
     do i=1,npfb
        tmp_vec(1) = u(i);tmp_vec(2) = v(i)
        tmp_vec(1) = -1.0d0;tmp_vec(2) = 0.0d0
        tmp_scal = dot_product(tmp_vec,grad(i,:))
        rhs(i) = lap(i)/Re - tmp_scal
     end do
     !$OMP END PARALLEL DO
     deallocate(lap,grad)

     return
  end subroutine calc_rhs_burgers
!! ------------------------------------------------------------------------------------------------
  subroutine sawtooth
     !! Subroutine to calculate the analytic velocity profile for sawtooth, and give comparison with numerical!
     !! only valid for 1D intitial conditions (v(x,y,t)=0 all x,y,t)
     real(rkind),dimension(:),allocatable :: A
     real(rkind) :: A0,sum_num,sum_denom,tmp_exp,x
     integer(ikind) :: i,n,nn

     !! How many terms in the sums?
     n = 29

     !! Load the bessel functions 
     !!(far easier than calculating them, and they're =0 by I29. NB, they need changing if Re.ne.100...
     allocate(A(n));
     open(unit=30,file='./data_out/bessels')
     read(30,*) A0
     do i=1,n
        read(30,*) A(i)
     end do
     close(30)
     
     !$OMP PARALLEL DO PRIVATE(sum_num,sum_denom,nn,x,tmp_exp)
     do i=1,npfb
        !! Calculate the analytical solution for each particle
        sum_num = 0.0d0;sum_denom = A0;
        x = rp(i,1)
        do nn=1,n
           tmp_exp = exp(-4.0*nn*nn*pi*pi*time/(Re*lambda**2.0))
           sum_num = sum_num + nn*A(nn)*sin(2.0*nn*pi*x/lambda)*tmp_exp
           sum_denom = sum_denom + A(nn)*cos(2.0*nn*pi*x/lambda)*tmp_exp
        end do
        !! Set w(i) with u(i)_numerical - u_analytical
        w(i) = u(i) - (4.0*pi/(lambda*Re))*sum_num/sum_denom
     end do
     !$OMP END PARALLEL DO
    
     deallocate(A)


  end subroutine sawtooth
!! ------------------------------------------------------------------------------------------------
  subroutine tanhwiggle
     !! Subroutine to calculate the analytic velocity profile for sawtooth, and give comparison with numerical!
     !! only valid for 1D intitial conditions (v(x,y,t)=0 all x,y,t)
     real(rkind) :: x,tmp
     integer(ikind) :: i,n,nn

    
     !$OMP PARALLEL DO PRIVATE(tmp,x)
     do i=1,npfb
        !! Calculate the analytical solution for each particle
        x = rp(i,1)
        tmp = -tanh(((x-0.5d0))*Re/2.0d0)
        w(i) = u(i)-tmp
     end do
     !$OMP END PARALLEL DO
    

  end subroutine tanhwiggle  

end module burgers_equation
