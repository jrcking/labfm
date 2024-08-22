module derivatives
!! This module contains routines to calculate the derivative of a function.
  use kind_parameters
  use common_parameter
  use common_2d
  use sphtools
  use analytic_functions
  use omp_lib
  implicit none
contains

  subroutine calc_gradientlog(phi,gradphi)
    !! Calculate the gradient of scalar phi.
    real(rkind),dimension(:),intent(in) :: phi
    real(rkind),dimension(:,:),intent(inout) :: gradphi
    real(rkind),dimension(:),allocatable :: logphi
    integer :: i,j,k
    real(rkind) :: fji
    real(rkind),dimension(dims) :: gradtmp
   
    allocate(logphi(np))
    !$omp parallel do
    do i=1,np
       logphi(i) = log(phi(i))
    end do
    !$omp end parallel do
 
    !$OMP PARALLEL DO PRIVATE(j,k,fji,gradtmp)
    do i=1,npfb
       gradtmp=0.0d0;
       do k=1,ij_count(i)
          j = ij_link(i,k) 
          fji = (logphi(j)-logphi(i))
          gradtmp(:) = gradtmp(:) + fji*ij_w_grad(i,k,:)
       end do
       gradphi(i,:) = phi(i)*gradtmp(:)
    end do
    !$OMP END PARALLEL DO

    deallocate(logphi)
    
    return
  end subroutine calc_gradientlog
!! ------------------------------------------------------------------------------------------------
  subroutine calc_gradient(phi,gradphi)
    !! Calculate the gradient of scalar phi.
    real(rkind),dimension(:),intent(in) :: phi
    real(rkind),dimension(:,:),intent(inout) :: gradphi
    integer(ikind) :: i,j,k
    real(rkind) :: fji
    real(rkind),dimension(dims) :: gradtmp
   
 
    !$OMP PARALLEL DO PRIVATE(j,k,fji,gradtmp)
    do i=1,npfb
       gradtmp=0.0d0;
       do k=1,ij_count(i)
          j = ij_link(i,k) 
          fji = (phi(j)-phi(i))
          gradtmp(:) = gradtmp(:) + fji*ij_w_grad(i,k,:)
       end do
       gradphi(i,:) = gradtmp(:)
    end do
    !$OMP END PARALLEL DO

    return
  end subroutine calc_gradient
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_integral(phi,intphi)
     !! Following the scheme of Hopkins, calculate the integral of a function
     real(rkind),dimension(:),intent(in) :: phi
     real(rkind),intent(inout) :: intphi
     integer(ikind) ::i,j,k 
     real(rkind) :: fji
     real(rkind) :: inttmp

     inttmp = 0.0d0  
     !$omp parallel do reduction(+:inttmp)
     do i=1,np
        inttmp = inttmp + phi(i)*hqw(i)
     end do
     !$omp end parallel do
     
     intphi = inttmp
  
  
     return
  end subroutine calc_integral    
!! ------------------------------------------------------------------------------------------------
  subroutine calc_laplacian(phi,lapphi)
    !! Calculate the Laplacian of a scalar phi
    real(rkind),dimension(:),intent(in) :: phi
    real(rkind),dimension(:),intent(inout) :: lapphi
    integer i,j,k
    real(rkind) :: fji,lap_tmp

    !$OMP PARALLEL DO PRIVATE(j,k,fji,lap_tmp)
    do i=1,npfb
       lap_tmp = 0.0d0
       do k=1,ij_count(i)
          j = ij_link(i,k) 
          fji = phi(j)-phi(i)
          lap_tmp = lap_tmp + fji*ij_w_lap(i,k)
       end do
       lapphi(i) = lap_tmp
    end do
    !$OMP END PARALLEL DO

  end subroutine calc_laplacian
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_laplacian_log(phi,lapphi)
    !! Calculate the Laplacian of a scalar phi
    real(rkind),dimension(:),intent(in) :: phi
    real(rkind),dimension(:),intent(inout) :: lapphi
    integer i,j,k
    real(rkind) :: fji,lap_tmp
    real(rkind),dimension(dims) :: grad_tmp

    !$OMP PARALLEL DO PRIVATE(j,k,fji,lap_tmp,grad_tmp)
    do i=1,npfb
       lap_tmp = 0.0d0;grad_tmp=0.0d0
       do k=1,ij_count(i)
          j = ij_link(i,k) 
          fji = log(phi(j)/phi(i))          
          lap_tmp = lap_tmp + fji*ij_w_lap(i,k)
          grad_tmp = grad_tmp + fji*ij_w_grad(i,k,:)
       end do
       lapphi(i) = lap_tmp*phi(i) + phi(i)*dot_product(grad_tmp,grad_tmp)
    end do
    !$OMP END PARALLEL DO

  end subroutine calc_laplacian_log
!! ------------------------------------------------------------------------------------------------ 
  subroutine calc_filtered_var(phi)
    !! Calculate the hyperviscosity filtered phi
    real(rkind),dimension(:),intent(inout) :: phi
    real(rkind),dimension(:),allocatable :: filtphi
    integer i,j,k
    real(rkind) :: fji,hyp_tmp

    !! Allocate temporary store
    allocate(filtphi(npfb))
    
    !! Calculate filtered phi
    !$OMP PARALLEL DO PRIVATE(j,k,fji,hyp_tmp)
    do i=1,npfb
       hyp_tmp = 0.0d0
       do k=1,ij_count(i)
          j = ij_link(i,k) 
          fji = phi(j)-phi(i)
          hyp_tmp = hyp_tmp + fji*ij_w_hyp(i,k)
       end do
       filtphi(i) = phi(i) + hyp_tmp                                                                  
    end do
    !$OMP END PARALLEL DO

    !! Pass filtered phi back to phi
    phi(1:npfb)=filtphi(1:npfb)

  end subroutine calc_filtered_var 
!! ------------------------------------------------------------------------------------------------  
end module derivatives
