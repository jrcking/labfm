module basic_convergence_studies
  use kind_parameters
  use common_parameter
  use common_2d
  use analytic_functions
  use derivatives
  use omp_lib
  implicit none

contains
  subroutine gradient_convergence_test
     integer(ikind) :: i,k,j
     real(rkind) :: tmp,x,y,tmp2
     real(rkind) :: l2_x,l2_y,e2_x,e2_y

     time = 0.0d0 !! No time integration here

     !! Find the gradient...
     allocate(grad(npfb,dims))
     call calc_gradient(u,grad)
!     call calc_gradientlog(u,grad)     

     !! Calculate the L2norms
     l2_x=0.0d0;l2_y=0.0d0;e2_x=0.0d0;e2_y=0.0d0
     !$OMP PARALLEL DO REDUCTION(+:l2_x,l2_y,e2_x,e2_y) PRIVATE(tmp,x,y,tmp2)
     do i=1,npfb
        x = rp(i,1);y = rp(i,2)
        !! Calculate the L2 norms
        tmp2 = dfdx(x,y);tmp = grad(i,1) - tmp2 
        l2_x = l2_x + tmp*tmp
        e2_x = e2_x + tmp2*tmp2
v(i)=tmp          
        tmp2 = dfdy(x,y);tmp = grad(i,2) - tmp2
        l2_y = l2_y + tmp*tmp
        e2_y = e2_y + tmp2*tmp2
!write(6,*) dfdx(x,y)/grad(i,1)        
     end do
     !$OMP END PARALLEL DO
     !! sqrt norms
     l2_x = sqrt(l2_x)/max(sqrt(e2_x),epsilon(e2_x))
     l2_y = sqrt(l2_y)/max(sqrt(e2_y),epsilon(e2_y))
     deallocate(grad)
     !! Output and flush
     write(1,*) h0/lambda,hovdx_av,l2_x,l2_y! !h0/lambda,hovdx_av  !kl/(nx*pi)
     write(6,*) npfb,hovdx_av,l2_x,l2_y
     flush(1)


     i=npfb/2
     do k=1,ij_count(i)
        j=ij_link(i,k)
        w(j) = 2.0d0
     end do
     w(i) = 1.0d0


     return
  end subroutine gradient_convergence_test
!! ------------------------------------------------------------------------------------------------
  subroutine laplacian_convergence_test
     integer(ikind) :: i,k
     real(rkind) :: tmp,x,y,tmp2
     real(rkind) :: l2_l,e2_l

     time = 0.0d0
 
     !! Find the laplacian... for three types of correction: none, BL, HO
     allocate(lap(npfb))
     call calc_laplacian(u,lap)
!     call calc_laplacian_log(u,lap)     

     !! Calculate the L2norms
     l2_l = 0.0d0;e2_l=0.0d0
     !$OMP PARALLEL DO REDUCTION(+:l2_l,e2_l) PRIVATE(tmp,x,y,tmp2)
     do i=nb+1,npfb
        x = rp(i,1);y=rp(i,2) 
        !! Calculate the L2 norm
        tmp2 = d2fdx2(x,y) + d2fdy2(x,y);tmp = lap(i) - tmp2
        l2_l = l2_l + tmp*tmp
v(i)=lap(i)
        e2_l = e2_l + tmp2*tmp2
     end do
     !$OMP END PARALLEL DO
     !! sqrt norms
     l2_l = sqrt(l2_l/npfb)/sqrt(e2_l/npfb)
     deallocate(lap)


     !! Output to file and flush
     write(1,*) h0/lambda,hovdx_av,l2_l !h0/lambda,hovdx_av !kl/(nx*pi)
     write(6,*) npfb,hovdx_av,l2_l
     flush(1)

     return
  end subroutine laplacian_convergence_test
!! ------------------------------------------------------------------------------------------------
  subroutine vortex_resolve_test
     !! Testing the ability of LABFM to resolve vortices, based on arguments in
     !! Mayrhofer et al. 2015
     integer(ikind) :: i,k
     real(rkind) :: tmp,x,y,tmp2
     real(rkind) :: l2_l,e2_l

     time = 0.0d0
 
     !! Set the initial conditions...
     allocate(grad(npfb,dims),v(np))
     allocate(w(np));w=0.0d0
     !$OMP PARALLEL DO PRIVATE(x,y)
     do i=1,np
        x=rp(i,1);y=rp(i,2)
        u(i) = sin(2.0*pi*x/lambda)*cos(2.0*pi*y/lambda)
        w(i) = 250.0d0*(cos(4.0*pi*x/lambda)+cos(4.0*pi*y/lambda))
     end do
     !$OMP END PARALLEL DO

     call calc_gradient(w,grad)
     !$OMP PARALLEL DO PRIVATE(x,y)
     do i=1,npfb
        v(i) = -2.0d0*u(i)*grad(i,1)
        x=rp(i,1);y=rp(i,2)       
        w(i) = 2.0d0*pi*1.0d3*sin(2.0*pi*x/lambda)*sin(4.0*pi*x/lambda)*cos(2.0*pi*y/lambda) !! analytic
     end do
     !$OMP END PARALLEL DO

     !! Calculate the L2norms
     l2_l = 0.0d0;e2_l=0.0d0
     !$OMP PARALLEL DO REDUCTION(+:l2_l,e2_l)
     do i=1,npfb
        !! Calculate the L2 norm
        u(i) = v(i) - w(i)
        l2_l = l2_l + u(i)*u(i)
        e2_l = e2_l + w(i)*w(i)
     end do
     !$OMP END PARALLEL DO
     !! sqrt norms
     l2_l = sqrt(l2_l/npfb)/sqrt(e2_l/npfb)
     deallocate(grad)


     !! Output to file and flush
     write(1,*) h0/lambda,hovdx_av,l2_l !h/lambda
     write(6,*) lambda*0.5d0/dx,l2_l !h/lambda
     flush(1)

     return
  end subroutine vortex_resolve_test
!! ------------------------------------------------------------------------------------------------
  subroutine stability_test
    integer(ikind) :: i,j,k,ij,nsize
    real(rkind) :: temp
    real(rkind),dimension(:,:),allocatable :: dm,vl
    real(rkind),dimension(:),allocatable :: lambda_re,lambda_im,work

    !! Build the global operator for d/dx
    if(allocated(irelation).eqv..false.)then !! Analytic boundaries...
    nsize=np;allocate(dm(nsize,nsize));dm=0.0d0
    do i=1,npfb
       do k=1,ij_count(i)
          j=ij_link(i,k)
          temp = -ij_w_grad(i,k,1)
!          temp = -ij_w_lap(i,k)
          if(i.ne.j)then
             dm(i,i) = dm(i,i) + temp
             dm(i,j) = dm(i,j) - temp
          end if
       end do
    end do
  
    else                                           !! Periodic boundaries...
    !! Build the LHS
    nsize=npfb;allocate(dm(nsize,nsize));dm=0.0d0
    do i=1,npfb
       do k=1,ij_count(i)
          j=ij_link(i,k)
          temp = -ij_w_grad(i,k,1)
          if(i.ne.j)then
             dm(i,i) = dm(i,i) + temp
             if(j.gt.npfb+nb_n)then    ! If j is mirror, then
                ij=irelation(j)    ! look for parent of j, ij
                if(i.ne.ij)then
                   dm(i,ij)=dm(i,ij)-temp
                else
                   dm(i,i)=dm(i,i)-temp
                end if
             elseif(j.le.npfb+nb_n)then   ! if j is not a mirror or a boundary node
                dm(i,j)=dm(i,j)-temp
             end if
          end if
       end do
    end do
    end if
write(6,*) "built DM"

    call openblas_set_num_threads(8)

    allocate(lambda_re(nsize),lambda_im(nsize),work(4*nsize));k=1
    allocate(vl(nsize,nsize));vl=0.0d0
    call dgeev("V","N",nsize,dm,nsize,lambda_re,lambda_im,vl,nsize,temp,1,work,4*nsize,k)
write(6,*) "found eigens"

    do i=1,npfb
       write(14,*) lambda_re(i),lambda_im(i)
       u(i) = abs(vl(i,1))
    end do

    deallocate(lambda_re,lambda_im,dm,vl,work)
    call openblas_set_num_threads(1)
stop
    return
  end subroutine stability_test
!! ------------------------------------------------------------------------------------------------  
  subroutine freq_response_test
     !! Calculate the frequency response of grad and lap operators
     !! ends with "stop" for single run. 
     !! if want multiple runs, comment out "stop", and ensure in main labfm file 
     !! that nx=nx*2 is commented out.
     !! plot from plots directory in Octave/matlab with "modal_response_G" and "modal_response_L"
     integer(ikind) :: i,k,j,kk
     real(rkind) :: tmp,x,y,fji,gsum,lsum,fsum
     real(rkind) :: l2_x,l2_l
     real(rkind),dimension(dims) :: rij

     do kk=1,floor(nx*pi),2  !! For each wavenumber (in steps of 2)...

        !! Calculate the L2norms
        l2_x=0.0d0;l2_l=0.0d0
        !$OMP PARALLEL DO PRIVATE(gsum,lsum,k,j,rij,fji,tmp,x,y) REDUCTION(+:l2_x,l2_l) 
        do i=1,npfb
           gsum=0.0d0;lsum=0.0d0
           do k=1,ij_count(i)
              j = ij_link(i,k) 
              rij = rp(j,:)-rp(i,:);x=rij(1);y=rij(2)
              fji = sin(kk*x)
              gsum = gsum + fji*ij_w_grad(i,k,1)
!              fji = 0.0625-0.0625*cos(kk*x)*cos(kk*y)
              fji = 0.5-0.5*cos(kk*x)*cos(kk*y)                    
              lsum = lsum + fji*ij_w_lap(i,k)
           end do
           tmp = gsum/(nx*pi)  !! For modal response
           l2_x = l2_x + tmp*tmp
           tmp = lsum/(nx*pi)**2
           l2_l = l2_l + tmp*tmp
           write(1,*) kk/(nx*pi),gsum/(nx*pi),lsum/(nx*pi)**2
        end do
        !$OMP END PARALLEL DO
        !! sqrt norms
        l2_x = sqrt(l2_x/npfb);l2_l=sqrt(l2_l/npfb)

        !! Output and flush
        write(1,*) kk/(nx*pi),l2_x,l2_l
        write(6,*) kk/(nx*pi),l2_x,l2_l
        flush(1)
     
     end do
     stop         
     return
  end subroutine freq_response_test
!! ------------------------------------------------------------------------------------------------  
  subroutine filter_test
     !! Calculate the frequency response a filter (via the hyperviscosity operator)
     !! ends with "stop" for single run. 
     !! if want multiple runs, comment out "stop", and ensure in main labfm file 
     !! that nx=nx*2 is commented out.
     
     !!!! Note that ij_w_hyp(i,k) is the hyperviscosity operator weight for the k-th neighbour of particle i     
     
     integer(ikind) :: i,k,j,kk
     real(rkind) :: tmp,x,y,fji,gsum,lsum,fsum,lscal
     real(rkind) :: l2_f
     real(rkind),dimension(dims) :: rij
     real(rkind) :: filter_coeff
     real(rkind),dimension(:),allocatable :: alpha_filter
     
     !! Determine the filter coefficients a priori
     allocate(alpha_filter(npfb))
     !$omp parallel do private(lsum,k,j,rij,fji,tmp,x,y,lscal)
     do i=1,npfb !! Loop over all particles
        lscal = 2.0*h(i)*sqrt(pi/dble(ij_count(i)))  !! Define lengthscale (effectively the particle spacing)
        tmp = (2.0/3.0)*pi/(lscal) !! (2/3)*Nyquist wavenumber
        lsum = 0.0d0
        do k=1,ij_count(i)              !! Loop over all neighbours
           j = ij_link(i,k)
           rij = rp(j,:)-rp(i,:);x=rij(1);y=rij(2)  

           fji = 1.0d0 - cos(tmp*x)*cos(tmp*y)     !! Test function is sinusoidal in x and y
           lsum = lsum + fji*ij_w_hyp(i,k)   
        end do
        alpha_filter(i) = (1.0/3.0)/lsum      !! Use response to test function to scale filter coefficient
     end do
     !$omp end parallel do
    
     
     !! Evaluate wavenumber response of filter
     do kk=1,floor(nx*pi),2  !! Loop over all wavenumbers from 1 to Nyquist, in steps of 2

        l2_f=0.0d0
        !$OMP PARALLEL DO PRIVATE(lsum,k,j,rij,fji,tmp,x,y,filter_coeff) REDUCTION(+:l2_f) 
        do i=1,npfb   !! Loop over all particles
           filter_coeff = alpha_filter(i)
           lsum=0.0d0
           do k=1,ij_count(i)    !! Loop over all neighbours
              j = ij_link(i,k) 
              rij = rp(j,:)-rp(i,:);x=rij(1);y=rij(2)

              fji = 1.0d0 - cos(kk*x)*cos(kk*y)   !! Test function is sinusoidal with wavenumber kk
              lsum = lsum + fji*ij_w_hyp(i,k)
           end do
           tmp = 1.0-lsum*filter_coeff    
           l2_f = l2_f + tmp*tmp
           write(1,*) kk/(nx*pi),tmp    !! Write response of THIS particle at THIS wavenumber to file
        end do
        !$OMP END PARALLEL DO
        !! sqrt norms
        l2_f=sqrt(l2_f/npfb)

        !! Output and flush
!        write(1,*) kk/(nx*pi),l2_f
        write(6,*) kk/(nx*pi),l2_f   ! Write L2 norm of response of all particles at THIS wavenumber to file.
        flush(1)
     
     end do
     
     deallocate(alpha_filter)
     stop         
     return
  end subroutine filter_test 
end module basic_convergence_studies
