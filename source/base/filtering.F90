module filtering
  use kind_parameters
  use common_parameter
  use common_2d
  use analytic_functions
  use derivatives
  use omp_lib
  use nodes
  implicit none

  real(rkind),dimension(:),allocatable :: alpha_filter

contains
!! ------------------------------------------------------------------------------------------------  
  subroutine filter_some_noise
     integer(ikind) :: ii,i
     !! Generate a noisy field and filter it, plotting some useful things along the way.
     
     call filter_config
     
     call create_noisy_field
  
     do ii=1,50
        call calc_some_stats(ii)
        call apply_filter_once
        
  
        call output_uv(ii)
     end do
  
  
     stop
  end subroutine filter_some_noise
!! ------------------------------------------------------------------------------------------------       
  subroutine create_noisy_field
     integer(ikind) :: i,j
     real(rkind) :: rand1,rand2,rnum
     
     do i=1,npfb
     
        !! Normally distributed
        rand1 = rand()
        rand2 = rand()
        rnum = sqrt(-2.0d0*log(rand1))*cos(2.0d0*pi*rand2)

        u(i) = 12.0d0 + 16.0d0*rnum!rand()-0.5d0
        w(i) = u(i)
        v(i) = 0.0d0
     end do
  
     do j=npfb+1,np
        u(j) = u(irelation(j))
     end do
  
     return
  end subroutine create_noisy_field
!! ------------------------------------------------------------------------------------------------       
  subroutine calc_some_stats(iter)
     integer(ikind),intent(in) :: iter
     integer(ikind) :: i,j
     real(rkind) :: meanu,varu,l2change

     meanu = 0.0d0
     do i=1,npfb
        meanu = meanu + u(i)
     end do
     meanu = meanu/dble(npfb)
     
     varu = 0.0d0
     do i=1,npfb
        varu = varu + (u(i)-meanu)**2.0d0
     end do
     varu = varu/dble(npfb)
     
     l2change = 0.0d0
     do i=1,npfb
        l2change = l2change + (w(i)-u(i))**2.0d0       
     end do
     l2change = sqrt(l2change/dble(npfb))
     
     
     write(6,*) iter,meanu,varu,l2change,maxval(u(1:npfb)),minval(u(1:npfb))
     write(15,*) iter,meanu,varu,l2change,maxval(u(1:npfb)),minval(u(1:npfb))     
   
  
     return
  end subroutine calc_some_stats
!! ------------------------------------------------------------------------------------------------    
  subroutine filter_config
     !! Calculate the frequency response a filter (via the hyperviscosity operator)
     !! ends with "stop" for single run. 
     !! if want multiple runs, comment out "stop", and ensure in main labfm file 
     !! that nx=nx*2 is commented out.
     
     !!!! Note that ij_w_hyp(i,k) is the hyperviscosity operator weight for the k-th neighbour of particle i     
     
     integer(ikind) :: i,k,j,kk
     real(rkind) :: tmp,x,y,fji,gsum,lsum,fsum,lscal,lmax
     real(rkind) :: l2_f
     real(rkind),dimension(dims) :: rij
     real(rkind) :: filter_coeff
     
     !! Determine the filter coefficients a priori
     allocate(alpha_filter(npfb))
if(.true.)then
     !$omp parallel do private(lsum,kk,k,j,rij,fji,tmp,x,y,lscal,lmax)
     do i=1,npfb !! Loop over all particles
        lscal = 2.0*h(i)*sqrt(pi/dble(ij_count(i)))  !! Define lengthscale (effectively the particle spacing)
 
        !! Loop over a range of wavenumbers and evaluate the effective wavenumber at each
        lmax = 0.0d0       
        do kk = 1,ceiling(pi/lscal),2
           tmp = kk
           lsum = 0.0d0
           do k=1,ij_count(i)              !! Loop over all neighbours
              j = ij_link(i,k)
              rij = rp(j,:)-rp(i,:);x=rij(1);y=rij(2)  

              fji = 1.0d0 - cos(tmp*x)*cos(tmp*y)     !! Test function is sinusoidal in x and y
              lsum = lsum + fji*ij_w_hyp(i,k)   
           end do
           lmax = max(lmax,lsum)
        end do

        !! Test specific wavenumber
        tmp = (5.0d0/6.0d0)*pi/lscal
        lsum = 0.0d0
        do k=1,ij_count(i)              !! Loop over all neighbours
           j = ij_link(i,k)
           rij = rp(j,:)-rp(i,:);x=rij(1);y=rij(2)  

           fji = 1.0d0 - cos(tmp*x)*cos(tmp*y)     !! Test function is sinusoidal in x and y
           lsum = lsum + fji*ij_w_hyp(i,k)   
        end do


!   if(lmax.gt.lsum) then  !! response at k=1 is not strongest
!      alpha_filter(i) = 0.5d0*(lmax+lsum)
!   else
!      alpha_filter(i) = lsum
!   end if 
   
   alpha_filter(i) = lsum
   alpha_filter(i) = alpha_filter(i)!/(two/three)     
        
!        alpha_filter(i) = lmax/(two/three)      !! Use response to test function to scale filter coefficient
!        write(6,*) h(i),lmax
     end do
     !$omp end parallel do
else
     alpha_filter(:)=1.0d0
endif     
                  
     return
  end subroutine filter_config 
!! ------------------------------------------------------------------------------------------------  
  subroutine apply_filter_once
     real(rkind) :: lsum,fji,x,y
     real(rkind),dimension(dims) :: rij
     integer(ikind) :: i,j,k
    
     !$OMP PARALLEL DO PRIVATE(lsum,k,j,rij,fji,x,y) 
     do i=1,npfb   !! Loop over all particles                   
         
        !! Store unfiltered u in w
        w(i) = u(i)   
           
        !! Evaluate filter           
        lsum = 0.0d0
        do k=1,ij_count(i)    !! Loop over all neighbours
           j = ij_link(i,k) 
           rij = rp(j,:)-rp(i,:);x=rij(1);y=rij(2)

           fji = u(i)-u(j)   !! Test function is sinusoidal with wavenumber kk
           lsum = lsum + fji*ij_w_hyp(i,k)/alpha_filter(i)
        end do
        
        !! Testing filters (basic is lsum = lsum)
        lsum = 1.0d0*lsum!*(abs(lsum/u(i)))**1.0d0
                
        !! Store filtered u in v           
        v(i) = u(i) - lsum
             
     end do
     !$OMP END PARALLEL DO
     
     !! Pass back across to u
     !$omp parallel do
     do i=1,npfb
        u(i) = v(i)
     end do
     !$omp end parallel do
     
     !! Mirrors
     do j=npfb+1,np
        u(j) = u(irelation(j))
     end do
     
     return
  end subroutine apply_filter_once  
end module filtering
