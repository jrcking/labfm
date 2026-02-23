module moments
  use kind_parameters
  use common_parameter
  use common_2d
  use sphtools
  use omp_lib
  use svd_lib
  implicit none
!! Choice of ABF type:: 1=original, 2=Hermite polynomials, 3=Legendre, 4=Laguerre, 5 = Taylor monomials
!! ABFs 2, 3 and 4 are multiplied by an RBF (Wab(qq) set in sphtools).
#define ABF 6

contains
!! ------------------------------------------------------------------------------------------------
!! ================================================================================================
  subroutine calc_interparticle_weights
     integer(ikind) :: i,j,k
     real(rkind) :: rad,qq,x,y,xx,yy,det,tmp,kval,rad2,rad3
     real(rkind),dimension(dims) :: gradw,rij

     !! Linear system to find ABF coefficients
     real(rkind),dimension(:,:),allocatable :: amatGx,amatGy,amatL,amathyp,amathyp2
     real(rkind),dimension(:),allocatable :: gvec,rvec
     real(rkind),dimension(:),allocatable :: xvec,bvecGx,bvecGy,bvecL,bvechyp,bvechyp2
     integer(ikind),dimension(:),allocatable :: ipiv
     integer(ikind) :: i1,i2,nsize,nsizeG
     real(rkind) :: ff1,hh,xs,ys,hyp_scaling
     real(rkind) :: ic,res1,res2
     
     real(rkind),dimension(27*4) :: work
     integer(ikind),dimension(27) :: iwork
     character*1 :: iflag
     real(rkind) :: cnum
     iflag = '1'
     
     ic = (0.0d0,1.0d0)

     !! Set desired order::
#if order==2
     k=2
#elif order==3
     k=3
#elif order==4
     k=4
#elif order==5
     k=5
#elif order==6
     k=6
#elif order==7
     k=7
#elif order==8
     k=8
#elif order==9
     k=9
#elif order==10
     k=10
#elif order==11
     k=11
#elif order==12
     k=12     
#endif
     nsizeG=(k*k+3*k)/2   !!  5,9,14,20,27,35,44... for k=2,3,4,5,6,7,8...

     !! Left hand sides and arrays for interparticle weights
     allocate(ij_w_grad(npfb,nplink,2),ij_w_lap(npfb,nplink),ij_w_hyp(npfb,nplink),ij_w_hyp2(npfb,nplink))
     ij_w_grad=0.0d0;ij_w_lap=0.0d0;ij_w_hyp=0.0d0
     allocate(amatGx(nsizeG,nsizeG),amatGy(nsizeG,nsizeG),amatL(nsizeG,nsizeG),amathyp(nsizeG,nsizeG),amathyp2(nsizeG,nsizeG))
     amatGx=0.0d0;amatGy=0.0d0;amatL=0.0d0;amathyp=0.0d0;amathyp2=0.0d0
 
     !! Right hand sides, vectors of monomials and ABFs
     allocate(bvecGx(nsizeG),bvecGy(nsizeG),bvecL(nsizeG),gvec(nsizeG),xvec(nsizeG),bvechyp(nsizeG),bvechyp2(nsizeG))
     bvecGx=0.0d0;bvecGy=0.0d0;bvecL=0.0d0;gvec=0.0d0;xvec=0.0d0;bvechyp=0.0d0;bvechyp2=0.0d0
     allocate(ipiv(nsizeG));ipiv=0
     allocate(rvec(nsizeG));rvec=0.0d0

#ifdef ob
     !! No parallelism for individual linear system solving...
     call openblas_set_num_threads(1)
#endif

     cnum = 0.0d0
     res1=zero
     res2=zero
     !$OMP PARALLEL DO PRIVATE(nsize,amatGx,amathyp,k,j,rij,rad,qq,x,y,xx,yy, &
     !$OMP ff1,gvec,xvec,i1,i2,amatL,amatGy,bvecGx,bvecGy,bvecL,bvechyp,hh,xs,ys,hyp_scaling, &
     !$omp work,iwork,ipiv,rvec,amathyp2,bvechyp2) &
     !$OMP reduction(+:cnum,res1,res2)
     do i=1,npfb
        nsize = nsizeG
        amatGx=0.0d0
        hh=h(i)
                
        do k=1,ij_count(i)
           j = ij_link(i,k) 
           rij(:) = rp(i,:) - rp(j,:)         
           x = -rij(1);y = -rij(2)
           

           !! Different types of ABF need different arguments (xx,yy)
           !! to account for domain of orthogonality
           rad = sqrt(x*x + y*y)/hh;qq=rad
#if ABF==1   
           ff1 = fac(qq)/h0/hh
           xx=x;yy=y
#elif ABF==2     
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Hermite   
#elif ABF==3
           ff1 = Wab(qq)
           xx=x/hh/ss;yy=y/hh/ss  !! Legendre
#elif ABF==4
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Laguerre   
#elif ABF==5
           ff1 = Wab(qq)      !! Taylor monomials
           xx = x/hh;yy=y/hh          
#elif ABF==6
           ff1 = Wab(qq)        !! Fourier
           xx=pi*x/hh/ss;yy=pi*y/hh/ss
#endif      
           !! Populate the ABF and monomial arrays
#if order>=2
           gvec(1:5) = abfs2(rad,xx,yy);
           xvec(1:5) = monomials2(x/hh,y/hh)
#endif     
#if order>=3
           gvec(6:9) = abfs3(rad,xx,yy);
           xvec(6:9) = monomials3(x/hh,y/hh)
#endif
#if order>=4
           gvec(10:14) = abfs4(rad,xx,yy);
           xvec(10:14) = monomials4(x/hh,y/hh)
#endif                       
#if order>=5
           gvec(15:20) = abfs5(rad,xx,yy);
           xvec(15:20) = monomials5(x/hh,y/hh)
#endif                       
#if order>=6
           gvec(21:27) = abfs6(rad,xx,yy);
           xvec(21:27) = monomials6(x/hh,y/hh)
#endif                       
#if order>=7
           gvec(28:35) = abfs7(rad,xx,yy);
           xvec(28:35) = monomials7(x/hh,y/hh)
#endif                       
#if order>=8
           gvec(36:44) = abfs8(rad,xx,yy);
           xvec(36:44) = monomials8(x/hh,y/hh)
#endif                       
#if order>=9
           gvec(45:54) = abfs9(rad,xx,yy);
           xvec(45:54) = monomials9(x/hh,y/hh) 
#endif   
#if order>=10
           gvec(55:65) = abfs10(rad,xx,yy);
           xvec(55:65) = monomials10(x/hh,y/hh) 
#endif            
#if order>=11
           gvec(66:77) = abfs11(rad,xx,yy);
           xvec(66:77) = monomials11(x/hh,y/hh) 
#endif            
#if order>=12
           gvec(78:90) = abfs12(rad,xx,yy);
           xvec(78:90) = monomials12(x/hh,y/hh) 
#endif            
                          
           !! Multiply ABFs by base RBF                         
           gvec(:) = gvec(:)*ff1  
                                             
           !! Build the LHS - it is the same for all three operators
           do i1=1,nsize
              amatGy(i1,:) = xvec(i1)*gvec(:)   !! Contribution to LHS for this interaction
           end do   
           amatGx(:,:) = amatGx(:,:) + amatGy(:,:)                                        
        end do   

!call dgecon(iflag,nsize,amatGx,nsize,qq,rad,work,iwork,i1)
!write(98,*) rp(i,1),1.0d0/rad
!        cnum = cnum + rad

        amatGy = amatGx;amatL=amatGx;amathyp = amatGx !! copying LHS
        amathyp2 = amathyp
        
        !! Build RHS for ddx and ddy
        bvecGx=0.0d0;bvecGx(1)=1.0d0!/hh          
        bvecGy=0.0d0;bvecGy(2)=1.0d0!/hh             
        
        !! Solve system for grad coefficients
        i1=0;i2=0                       
#ifdef ob
        call dgesv(nsize,1,amatGx,nsize,i1,bvecGx,nsize,i2)   
#else        
        call svd_solve(amatGx,nsize,bvecGx)       
#endif        

        i1=0;i2=0;nsize=nsizeG    
#ifdef ob
        call dgesv(nsize,1,amatGy,nsize,i1,bvecGy,nsize,i2)    
#else        
        call svd_solve(amatGy,nsize,bvecGy)       
#endif


        !! Solve system for Lap coefficients
        bvecL(:)=0.0d0;bvecL(3)=1.0d0;bvecL(5)=1.0d0;i1=0;i2=0;nsize=nsizeG
#ifdef ob
        call dgesv(nsize,1,amatL,nsize,i1,bvecL,nsize,i2)
#else        
        call svd_solve(amatL,nsize,bvecL)
#endif        
        
        !! Solve system for Hyperviscosity (regular viscosity if ORDER<4)
#if order<=3
        bvechyp(:)=0.0d0;bvechyp(3)=1.0d0;bvechyp(5)=1.0d0
        hyp_scaling = one/hh/hh
#elif order<=5
        bvechyp(:)=0.0d0;bvechyp(10)=-1.0d0;bvechyp(12)=-2.0d0;bvechyp(14)=-1.0d0
        hyp_scaling = one/hh/hh/hh/hh
#elif order<=7
        bvechyp(:)=0.0d0;bvechyp(21)=1.0d0;bvechyp(23)=3.0d0;bvechyp(25)=3.0d0;bvechyp(27)=1.0d0
        hyp_scaling = one/hh/hh/hh/hh/hh/hh        
#elif order<=9
!        bvechyp(:)=0.0d0;bvechyp(36)=-1.0d0;bvechyp(38)=-4.0d0;bvechyp(40)=-6.0d0;bvechyp(42)=-4.0d0;bvechyp(44)=-1.0d0
        bvechyp(:)=0.0d0;bvechyp(36)=-1.0d0;bvechyp(38)=-0.0d0;bvechyp(40)=-1.0d0;bvechyp(42)=-0.0d0;bvechyp(44)=-1.0d0
        hyp_scaling = one/hh/hh/hh/hh/hh/hh/hh/hh        


        bvechyp2(:)=0.0d0;bvechyp2(21)=1.0d0;bvechyp2(23)=3.0d0;bvechyp2(25)=3.0d0;bvechyp2(27)=1.0d0
!        bvechyp2(:)=0.0d0;bvechyp2(21)=1.0d0;bvechyp2(23)=3.0d0;bvechyp2(25)=3.0d0;bvechyp2(27)=1.0d0
#else
        bvechyp(:)=0.0d0;bvechyp(55)=1.0d0;bvechyp(57)=5.0d0;bvechyp(59)=10.0d0;bvechyp(61)=10.0d0
        bvechyp(63)=5.0d0;bvechyp(65)=1.0d0
        hyp_scaling = one/hh/hh/hh/hh/hh/hh/hh/hh/hh/hh                
#endif
        i1=0;i2=0;nsize=nsizeG 
#ifdef ob        
        call dgesv(nsize,1,amathyp,nsize,i1,bvechyp,nsize,i2)   
#else        
        call svd_solve(amathyp,nsize,bvechyp)             
#endif        
        i1=0;i2=0;nsize=nsizeG 
#ifdef ob        
        call dgesv(nsize,1,amathyp2,nsize,i1,bvechyp2,nsize,i2)   
#else        
        call svd_solve(amathyp2,nsize,bvechyp2)             
#endif        
        !! Another loop of neighbours to calculate interparticle weights
        do k=1,ij_count(i)
           j = ij_link(i,k) 
           rij(:) = rp(i,:) - rp(j,:)          
           x=-rij(1);y=-rij(2)
          
           !! Calculate arguments (diff ABFs need args over diff ranges)
           rad = sqrt(x*x + y*y)/hh;qq=rad           
#if ABF==1   
           ff1 = fac(qq)/h0/hh
           xx=x;yy=y
#elif ABF==2     
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Hermite   
#elif ABF==3
           ff1 = Wab(qq)
           xx=x/hh/ss;yy=y/hh/ss  !! Legendre
#elif ABF==4
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Laguerre   
#elif ABF==5
           ff1 = Wab(qq)      !! Taylor monomials
           xx = x/hh;yy=y/hh    
#elif ABF==6
           ff1 = Wab(qq)       !! Fourier
           xx=pi*x/hh/ss;yy=pi*y/hh/ss
#endif          
           !! Populate the ABF array        
#if order>=2
           gvec(1:5) = abfs2(rad,xx,yy)
#endif     
#if order>=3
           gvec(6:9) = abfs3(rad,xx,yy)
#endif
#if order>=4
           gvec(10:14) = abfs4(rad,xx,yy)
#endif                       
#if order>=5
           gvec(15:20) = abfs5(rad,xx,yy)
#endif                       
#if order>=6
           gvec(21:27) = abfs6(rad,xx,yy)
#endif                       
#if order>=7
           gvec(28:35) = abfs7(rad,xx,yy)
#endif                       
#if order>=8
           gvec(36:44) = abfs8(rad,xx,yy)
#endif                       
#if order>=9
           gvec(45:54) = abfs9(rad,xx,yy)
#endif   
#if order>=10
           gvec(55:65) = abfs10(rad,xx,yy)
#endif 
#if order>=11
           gvec(66:77) = abfs11(rad,xx,yy)
#endif 
#if order>=12
           gvec(78:90) = abfs12(rad,xx,yy)
#endif 

           !! Multiply ABFs by base RBF                         
           gvec(:) = gvec(:)*ff1 

           !! Weights for operators        
           ij_w_grad(i,k,1) = dot_product(bvecGx,gvec)/hh
           ij_w_grad(i,k,2) = dot_product(bvecGy,gvec)/hh
           ij_w_lap(i,k) = dot_product(bvecL,gvec)/hh/hh 
           ij_w_hyp(i,k) = dot_product(bvechyp,gvec)!*hyp_scaling   
           ij_w_hyp2(i,k) = dot_product(bvechyp2,gvec)!/hyp_scaling
                 
        end do
#ifdef ob        
        write(2,*) ""  !! Temporary fix for a weird openblas/lapack/openmp bug 
        !! The bug may be due to a race condition or something, but not really sure
        !! writing nothing to fort.2 ( I guess) is akin to a pause. Discovered by trial
        !! and error. Results in a big fort.2 file, and should probably find a proper solution 
        !! sometime.
#endif        

     end do
     !$OMP END PARALLEL DO
         
         
!     write(6,*) 1.0d0/(cnum/dble(npfb))

     deallocate(amatGx,amatGy,amatL,bvecGx,bvecGy,bvecL,gvec,xvec,bvechyp,amathyp)     
     return
  end subroutine calc_interparticle_weights
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_interparticle_weights_nofilt_o2
     integer(ikind) :: i,j,k
     real(rkind) :: rad,qq,x,y,xx,yy,det,tmp,kval,rad2,rad3
     real(rkind),dimension(dims) :: gradw,rij

     !! Linear system to find ABF coefficients
     real(rkind),dimension(:,:),allocatable :: amatGx,amatGy,amatL
     real(rkind),dimension(:),allocatable :: gvec,rvec
     real(rkind),dimension(:),allocatable :: xvec,bvecGx,bvecGy,bvecL
     integer(ikind),dimension(:),allocatable :: ipiv
     integer(ikind) :: i1,i2,nsize,nsizeG
     real(rkind) :: ff1,hh,xs,ys
     real(rkind) :: ic,res1,res2

     !! Set desired order::
     k=2

     nsizeG=(k*k+3*k)/2   !!  5,9,14,20,27,35,44... for k=2,3,4,5,6,7,8...

     !! Left hand sides and arrays for interparticle weights
     allocate(ij_w_grad(npfb,nplink,2),ij_w_lap(npfb,nplink))
     ij_w_grad=0.0d0;ij_w_lap=0.0d0
     allocate(amatGx(nsizeG,nsizeG),amatGy(nsizeG,nsizeG),amatL(nsizeG,nsizeG))
     amatGx=0.0d0;amatGy=0.0d0;amatL=0.0d0
 
     !! Right hand sides, vectors of monomials and ABFs
     allocate(bvecGx(nsizeG),bvecGy(nsizeG),bvecL(nsizeG),gvec(nsizeG),xvec(nsizeG))
     bvecGx=0.0d0;bvecGy=0.0d0;bvecL=0.0d0;gvec=0.0d0;xvec=0.0d0

     !$OMP PARALLEL DO PRIVATE(nsize,amatGx,k,j,rij,rad,qq,x,y,xx,yy, &
     !$OMP ff1,gvec,xvec,i1,i2,amatL,amatGy,bvecGx,bvecGy,bvecL,hh,xs,ys) 

     do i=1,npfb
        nsize = nsizeG
        amatGx=0.0d0
        hh=h(i)
                
        do k=1,ij_count(i)
           j = ij_link(i,k) 
           rij(:) = rp(i,:) - rp(j,:)         
           x = -rij(1);y = -rij(2)
           

           !! Different types of ABF need different arguments (xx,yy)
           !! to account for domain of orthogonality           
           rad = sqrt(x*x + y*y)/hh;qq=rad
#if ABF==1   
           ff1 = fac(qq)/hh/hh
           xx=x;yy=y
#elif ABF==2     
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Hermite   
#elif ABF==3
           ff1 = Wab(qq)
           xx=x/hh/ss;yy=y/hh/ss  !! Legendre
#elif ABF==4
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Laguerre   
#elif ABF==5
           ff1 = Wab(qq)      !! Taylor monomials
           xx = x/hh;yy=y/hh    
#elif ABF==6
           ff1 = Wab(qq)       !! Fourier
           xx=pi*x/hh/ss;yy=pi*y/hh/ss
#endif            

           !! Populate the ABF and monomial arrays
           gvec(1:5) = abfs2(rad,xx,yy);
           xvec(1:5) = monomials2(x/hh,y/hh)

                                  
           !! Multiply ABFs by base RBF                         
           gvec(:) = gvec(:)*ff1  
                                             
           !! Build the LHS - it is the same for all three operators
           do i1=1,nsize
              amatGy(i1,:) = xvec(i1)*gvec(:)   !! Contribution to LHS for this interaction
           end do   
           amatGx(:,:) = amatGx(:,:) + amatGy(:,:)                                        
        end do   


        amatGy = amatGx;amatL=amatGx!! copying LHS

        
        !! Build RHS for ddx and ddy
        bvecGx=0.0d0;bvecGx(1)=1.0d0!/hh          
        bvecGy=0.0d0;bvecGy(2)=1.0d0!/hh             
        
        !! Solve system for grad coefficients
        i1=0;i2=0                          
        call svd_solve(amatGx,nsize,bvecGx)       

        i1=0;i2=0;nsize=nsizeG         
        call svd_solve(amatGy,nsize,bvecGy)       


        !! Solve system for Lap coefficients
        bvecL(:)=0.0d0;bvecL(3)=1.0d0;bvecL(5)=1.0d0;i1=0;i2=0;nsize=nsizeG     
        call svd_solve(amatL,nsize,bvecL)
             
        !! Another loop of neighbours to calculate interparticle weights
        do k=1,ij_count(i)
           j = ij_link(i,k) 
           rij(:) = rp(i,:) - rp(j,:)          
           x=-rij(1);y=-rij(2)
          
           !! Calculate arguments (diff ABFs need args over diff ranges)
           rad = sqrt(x*x + y*y)/hh;qq=rad     
#if ABF==1   
           ff1 = fac(qq)/hh/hh
           xx=x;yy=y
#elif ABF==2     
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Hermite   
#elif ABF==3
           ff1 = Wab(qq)
           xx=x/hh/ss;yy=y/hh/ss  !! Legendre
#elif ABF==4
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Laguerre   
#elif ABF==5
           ff1 = Wab(qq)      !! Taylor monomials
           xx = x/hh;yy=y/hh    
#elif ABF==6
           ff1 = Wab(qq)       !! Fourier
           xx=pi*x/hh/ss;yy=pi*y/hh/ss
#endif                  
           !! Populate the ABF array        
           gvec(1:5) = abfs2(rad,xx,yy)


           !! Multiply ABFs by base RBF                         
           gvec(:) = gvec(:)*ff1 

           !! Weights for operators        
           ij_w_grad(i,k,1) = dot_product(bvecGx,gvec)/hh
           ij_w_grad(i,k,2) = dot_product(bvecGy,gvec)/hh
           ij_w_lap(i,k) = dot_product(bvecL,gvec)/hh/hh 

                 
        end do     

     end do
     !$OMP END PARALLEL DO
     deallocate(amatGx,amatGy,amatL,bvecGx,bvecGy,bvecL,gvec,xvec)          
     
   
     return
  end subroutine calc_interparticle_weights_nofilt_o2 
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_interparticle_weights_nofilt_sph
     integer(ikind) :: i,j,k
     real(rkind) :: rad,qq,x,y,xx,yy,det,tmp,kval,rad2,rad3
     real(rkind),dimension(dims) :: gradw,rij

     !! Linear system to find ABF coefficients
     real(rkind),dimension(:,:),allocatable :: amatGx,amatGy,amatL
     real(rkind),dimension(:),allocatable :: gvec,rvec
     real(rkind),dimension(:),allocatable :: xvec,bvecGx,bvecGy,bvecL
     integer(ikind),dimension(:),allocatable :: ipiv
     integer(ikind) :: i1,i2,nsize,nsizeG
     real(rkind) :: ff1,hh,xs,ys
     real(rkind) :: ic,res1,res2

     !! Bonet & Lok
     nsizeG=2   !!  5,9,14,20,27,35,44... for k=2,3,4,5,6,7,8...

     !! Left hand sides and arrays for interparticle weights
     allocate(ij_w_grad(npfb,nplink,2),ij_w_lap(npfb,nplink))
     ij_w_grad=0.0d0;ij_w_lap=0.0d0
     allocate(amatGx(nsizeG,nsizeG),amatGy(nsizeG,nsizeG),amatL(nsizeG,nsizeG))
     amatGx=0.0d0;amatGy=0.0d0;amatL=0.0d0
 
     !! Right hand sides, vectors of monomials and ABFs
     allocate(bvecGx(nsizeG),bvecGy(nsizeG),bvecL(nsizeG),gvec(nsizeG),xvec(nsizeG))
     bvecGx=0.0d0;bvecGy=0.0d0;bvecL=0.0d0;gvec=0.0d0;xvec=0.0d0

     !$OMP PARALLEL DO PRIVATE(nsize,amatGx,k,j,rij,rad,qq,x,y,xx,yy, &
     !$OMP ff1,gvec,xvec,i1,i2,amatL,amatGy,bvecGx,bvecGy,bvecL,hh,xs,ys,gradw) 

     do i=1,npfb
        nsize = nsizeG
        amatGx=zero
        hh=h(i)
                
        do k=1,ij_count(i)
           j = ij_link(i,k) 
           rij(:) = rp(i,:) - rp(j,:)         
           x = rij(1);y = rij(2)
           

           !! Different types of ABF need different arguments (xx,yy)
           !! to account for domain of orthogonality
           rad = sqrt(x*x + y*y);qq=rad/hh
           rad = max(rad,hh*1e-5)

           ff1 = fac(qq)/(hh*hovdx*hovdx)
           gradw(:) = -rij(:)*ff1/rad   !uncorrected kernel gradient           

           !! Bonet & Lok correction tensor
           amatGx(1,:) = amatGx(1,:) + gradw(:)*x
           amatGx(2,:) = amatGx(2,:) + gradw(:)*y           

           !! Uncorrected gradient
           ij_w_grad(i,k,:) = -gradw(:)
           ij_w_lap(i,k) = -two*dot_product(rij,ij_w_grad(i,k,:))/rad/rad

        end do

        !! Inverse of tensor
        qq = amatGx(1,1)*amatGx(2,2) - amatGx(1,2)*amatGx(2,1)  
                                      
        amatGy(1,1) = amatGx(2,2)/qq
        amatGy(2,2) = amatGx(1,1)/qq
        amatGy(1,2) =-amatGx(1,2)/qq
        amatGy(2,1) =-amatGx(2,1)/qq
           
        do k=1,ij_count(i)
           j = ij_link(i,k) 
           
           gradw(:) = ij_w_grad(i,k,:)
                      
           !! Correct gradient:
           ij_w_grad(i,k,:) = matmul(amatGy(:,:),gradw(:))
                  
        end do   

     end do
     !$OMP END PARALLEL DO
!     write(6,*) 1.0d0/(cnum/dble(npfb))

     deallocate(amatGx,amatGy,amatL,bvecGx,bvecGy,bvecL,gvec,xvec)     
     return
  end subroutine calc_interparticle_weights_nofilt_sph   
!! ------------------------------------------------------------------------------------------------
  subroutine filter_coefficients
     !! Determine filter coefficients for hyperviscosity a priori, requiring
     !! 2/3 damping at wavenumbers 2/3 of Nyquist...
     integer(ikind) :: i,j,k
     real(rkind) :: fji,lsum,x,y,tmp,lscal
     real(rkind),dimension(dims) :: rij

     !! Allocate the coefficients
     allocate(filter_coeff(npfb))
     
     !$omp parallel do private(lsum,k,j,rij,fji,tmp,x,y,lscal)
     do i=1,npfb
        !! Set the particle length-scale (=~dx, but we *shouldn't* have access to dx)
        lscal = 2.0*h(i)*sqrt(pi/dble(ij_count(i)))
        !! Set the target wavenumber (2/3 of Nyquist)
        tmp = (5.0d0/6.0d0)*pi/lscal !1.5
        
        !! Calculate hyperviscosity operator of a signal at this wavenumber
        lsum = 0.0d0
        do k=1,ij_count(i)
           j = ij_link(i,k)
           rij = rp(j,:)-rp(i,:);x=rij(1);y=rij(2)
           fji = 1.0d0 - cos(tmp*x)*cos(tmp*y) 
           lsum = lsum + fji*ij_w_hyp(i,k)
        end do
        
        !! Set the filter coefficient
        filter_coeff(i) = (5.0d0/6.0d0)/lsum   !2/3
        
        ij_w_hyp(i,:) = ij_w_hyp(i,:)*filter_coeff(i)
     end do
     !$omp end parallel do
     
     
    
     return
  end subroutine filter_coefficients
!! ------------------------------------------------------------------------------------------------
  subroutine adapt_stencils
#ifndef ob  
     !! This subroutine refines the stencil sizes. For each node, "smoothing length" is reduced by 
     !! 2% per iteration, and this continues until several checks are failed:
     !! 1) residual of LABFM laplacian system too big.
     !! 2) Amplification factor of d/dx,d/dy,Laplacian > 1+tol, for N wavenumbers up to Nyquist
     !! 3) h is smaller than X% of initial h.
     !! 
     !! Currently does this calculation for first layer of nodes, then if 3D, copies new smoothing 
     !! lengths to other layers, and builds new stencils.
     integer(ikind) :: i,j,k,ii,jj,nk,jjj
     real(rkind) :: rad,qq,x,y,xx,yy,hchecksum,hchecksumL,hchecksumX,hchecksumY
     real(rkind),dimension(dims) :: rij
     real(rkind),dimension(:,:),allocatable :: amat,mmat
     real(rkind),dimension(:),allocatable :: bvecL,bvecX,bvecY,gvec,xvec
     integer(ikind) :: i1,i2,nsize,nsizeG,n_to_reduce,full_j_count_i
     real(rkind) :: ff1,hh,reduction_factor,res_tol,amp_tol,sumcheck
     logical :: reduce_h
     integer(ikind),dimension(:),allocatable :: full_j_link_i
     real(rkind),dimension(dims) :: grad_s
     real(rkind) :: grad_s_mag,hfactor,radmax
     real(rkind),dimension(:),allocatable :: neighbourcountreal
 
#if order==4
     k=4
#elif order==6
     k=6
#elif order==8
     k=8
#elif order==10
     k=10
#endif
     nsizeG=(k*k+3*k)/2   !!  5,9,14,20,27,35,44... for k=2,3,4,5,6,7,8...

     !! Left hand sides 
     allocate(amat(nsizeG,nsizeG),mmat(nsizeG,nsizeG))
     amat=zero;mmat=zero
 
     !! Right hand sides, vectors of monomials and ABFs
     allocate(bvecL(nsizeG),bvecX(nsizeG),bvecY(nsizeG),gvec(nsizeG),xvec(nsizeG))
     bvecL=zero;bvecX=zero;bvecY=zero;gvec=zero;xvec=zero

    
     !! Temporary neighbour lists...
     allocate(full_j_link_i(nplink));full_j_link_i=0
     allocate(neighbourcountreal(npfb));neighbourcountreal=zero

     !! Set parameters of h-reduction
     reduction_factor = 0.98 
#if order==4
     res_tol = 1.0d-3*dble(nsizeG**4)*epsilon(hchecksum)/dble(k)   !! For 6th order
#elif order==6
     res_tol = 5.0d-3*dble(nsizeG**4)*epsilon(hchecksum)/dble(k)   !! For 6th order
#elif order==8
     res_tol = 4.0d-2*dble(nsizeG**4)*epsilon(hchecksum)/dble(k)   !! For 8th order    
#elif order==10     
     res_tol = 1.0d+0*dble(nsizeG**4)*epsilon(hchecksum)/dble(k)   !! For 10th order
#endif     
#if ABF==6
     res_tol = res_tol*1.0d-3
#endif

     nk = 32   !! How many wavenumbers between 1 and Nyquist to check... ! 16
     amp_tol = 1.0001d0   !! Maximum allowable amplification

   
     sumcheck = zero
!     !$OMP PARALLEL DO PRIVATE(nsize,amat,k,j,rij,rad,qq,x,y,xx,yy, &
!     !$OMP ff1,gvec,xvec,i1,i2,bvecL,bvecX,bvecY,hh,full_j_count_i,full_j_link_i, &
!     !$OMP hchecksum,reduce_h,ii,mmat,hchecksumL,hchecksumX,hchecksumY) &
!     !$omp reduction(+:sumcheck)
     do i=1,npfb
        ii=0
        reduce_h=.true.
        do while (reduce_h)
           !! Reduce h (temporarily stored in hh
           hh=h(i)*reduction_factor 
           
           !! Build temporary neighbour lists
           full_j_count_i=0
           full_j_link_i(:)=0
           do k=1,ij_count(i)
              j=ij_link(i,k)
              rij(:) = rp(i,:) - rp(j,:)              
              rad = sqrt(dot_product(rij,rij))
              if(rad.le.hh*ss)then               !! Only for nodes within new support radius
                 full_j_count_i = full_j_count_i + 1
                 full_j_link_i(full_j_count_i) = j
              end if              
           end do
 
           !! Build linear system 
           nsize = nsizeG
           amat=zero
           do k=1,full_j_count_i
              j = full_j_link_i(k) 
              rij(:) = rp(i,:) - rp(j,:)
              x = -rij(1);y = -rij(2)

              !! Different types of ABF need different arguments (xx,yy)
              !! to account for domain of orthogonality
#if ABF==1
              !! Find dW/dr of fundamental RBF
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1 = fac(qq)/h0/hh
              xx=x;yy=y
#elif ABF==2
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
              xx=x/hh;yy=y/hh    !! Hermite          
#elif ABF==3
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
              xx=x/hh/ss;yy=y/hh/ss  !! Legendre  
#elif ABF==6
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function              
              xx=pi*x/hh/ss;yy=pi*y/hh/ss               
#endif     
           !! Populate the ABF and monomial arrays
#if order>=2
              gvec(1:5) = abfs2(rad,xx,yy);
              xvec(1:5) = monomials2(x/hh,y/hh)
#endif     
#if order>=3
              gvec(6:9) = abfs3(rad,xx,yy);
              xvec(6:9) = monomials3(x/hh,y/hh)
#endif
#if order>=4
              gvec(10:14) = abfs4(rad,xx,yy);
              xvec(10:14) = monomials4(x/hh,y/hh)
#endif                       
#if order>=5
              gvec(15:20) = abfs5(rad,xx,yy);
              xvec(15:20) = monomials5(x/hh,y/hh)
#endif                       
#if order>=6
              gvec(21:27) = abfs6(rad,xx,yy);
              xvec(21:27) = monomials6(x/hh,y/hh)
#endif                       
#if order>=7
              gvec(28:35) = abfs7(rad,xx,yy);
              xvec(28:35) = monomials7(x/hh,y/hh)
#endif                       
#if order>=8
              gvec(36:44) = abfs8(rad,xx,yy);
              xvec(36:44) = monomials8(x/hh,y/hh)
#endif                       
#if order>=9
              gvec(45:54) = abfs9(rad,xx,yy);
              xvec(45:54) = monomials9(x/hh,y/hh) 
#endif   
#if order>=10
              gvec(55:65) = abfs10(rad,xx,yy);
              xvec(55:65) = monomials10(x/hh,y/hh) 
#endif            
#if order>=11
              gvec(66:77) = abfs11(rad,xx,yy);
              xvec(66:77) = monomials11(x/hh,y/hh) 
#endif            
#if order>=12
              gvec(78:90) = abfs12(rad,xx,yy);
              xvec(78:90) = monomials12(x/hh,y/hh) 
#endif    
              gvec(1:nsizeG) = gvec(1:nsizeG)*ff1                            


              !! Build the LHS - it is the same for all three diff operators (ddx,ddy,Lap)
              do i1=1,nsize
                 amat(i1,:) = amat(i1,:) + xvec(i1)*gvec(:)
              end do
           end do
           mmat = amat

           !! Solve system for Laplacian
           bvecL(:)=zero;bvecL(3)=one/hh/hh;bvecL(5)=one/hh/hh;i1=0;i2=0
           call svd_solve(amat,nsize,bvecL)       

           !! Solve system for d/dx           
           amat=mmat
           bvecX(:)=zero;bvecX(1)=one/hh;i1=0;i2=0;nsize=nsizeG           
           call svd_solve(amat,nsize,bvecX)       

           !! Solve system for d/dy
           amat=mmat
           bvecY(:)=zero;bvecY(2)=one/hh;i1=0;i2=0;nsize=nsizeG           
           call svd_solve(amat,nsize,bvecY)       
           
           !! First (main) check for h-reduction: the residual of the linear system solution (Laplacian)
           !! Multiply the solution vector by the LHS, and calculate the residual
           xvec=zero;xvec(3)=-one/hh/hh;xvec(5)=-one/hh/hh !! Initialise with -C^{L} (LABFM paper notation)
           hchecksum = zero
           do i1=1,nsizeG
              do i2=1,nsizeG
                 xvec(i1) = xvec(i1) + mmat(i1,i2)*bvecL(i2)
              end do
              hchecksum = hchecksum + xvec(i1)**two
           end do
           hchecksum = hchecksum*hh*hh
           hchecksum = sqrt(hchecksum/dble(nsizeG))
           
           !! Second check for h-reduction: amplification of wavenumbers below Nyquist
           i1=0
           if(hchecksum.lt.res_tol/hh) then
              do while(i1.le.nk.and.hchecksum.lt.res_tol/hh)
                 i1 = i1 + 1
                 hchecksumL = zero
                 hchecksumX = zero
                 hchecksumY = zero
                 do k=1,full_j_count_i
                    j = full_j_link_i(k) 
                    rij(:) = rp(i,:) - rp(j,:)
                    x = -rij(1);y = -rij(2)

                    !! Different types of ABF need different arguments (xx,yy)
                    !! to account for domain of orthogonality
#if ABF==1
                 !! Find dW/dr of fundamental RBF
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1 = fac(qq)/h0/hh
                    xx=x;yy=y
#elif ABF==2
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
                    xx=x/hh;yy=y/hh    !! Hermite          
#elif ABF==3
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
                    xx=x/hh/ss;yy=y/hh/ss  !! Legendre   
#elif ABF==6                    
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function                    
                    xx=pi*x/hh/ss;yy=pi*y/hh/ss                                   
#endif     
                    !! Populate the ABF and monomial arrays
#if order>=2
                    gvec(1:5) = abfs2(rad,xx,yy);
#endif     
#if order>=3
                    gvec(6:9) = abfs3(rad,xx,yy);
#endif
#if order>=4
                    gvec(10:14) = abfs4(rad,xx,yy);
#endif                       
#if order>=5
                    gvec(15:20) = abfs5(rad,xx,yy);
#endif                       
#if order>=6
                    gvec(21:27) = abfs6(rad,xx,yy);
#endif                       
#if order>=7
                    gvec(28:35) = abfs7(rad,xx,yy);
#endif                       
#if order>=8
                    gvec(36:44) = abfs8(rad,xx,yy);
#endif                       
#if order>=9
                    gvec(45:54) = abfs9(rad,xx,yy);
#endif   
#if order>=10
                    gvec(55:65) = abfs10(rad,xx,yy);
#endif            
#if order>=11
                    gvec(66:77) = abfs11(rad,xx,yy);
#endif            
#if order>=12
                    gvec(78:90) = abfs12(rad,xx,yy);
#endif    
                    gvec(1:nsizeG) = gvec(1:nsizeG)*ff1                    


                    !! grad and Laplacian of signal at particular wavenumber qq
                    qq = dble(i1)*pi/(h(i)/hovdx)/dble(nk)
                    hchecksumL = hchecksumL + (0.5d0 - 0.5d0*cos(x*qq)*cos(y*qq))*dot_product(bvecL,gvec)
                    hchecksumX = hchecksumX + cos(y*qq)*sin(x*qq)*dot_product(bvecX,gvec)
                    hchecksumY = hchecksumY + cos(x*qq)*sin(y*qq)*dot_product(bvecY,gvec)
                            
                 end do                         
           
                 !! Normalise with qq or qq**2
                 hchecksumL = hchecksumL/(qq**2)
                 hchecksumX = hchecksumX/qq
                 hchecksumY = hchecksumY/qq

                 !! Modify hchecksum to break out of h-reduction loop if required
                 if(hchecksumL.gt.amp_tol.or.hchecksumL.lt.zero) then
write(6,*) i,i1,"stopping because of L",ii,hchecksum,res_tol/hh 
                    hchecksum = two*res_tol/hh              
                 end if
                 if(abs(hchecksumX).gt.amp_tol) then
                    hchecksum = two*res_tol/hh
write(6,*) i,i1,"stopping because of X",ii                 
                 end if
                 if(abs(hchecksumY).gt.amp_tol) then
                    hchecksum = two*res_tol/hh
write(6,*) i,i1,"stopping because of Y",ii
                 end if
                 if(isnan(hchecksum)) then
                    hchecksum = two*res_tol/hh
write(6,*) i,i1,"stopping because of NaN",ii                    
                 end if
              
              end do
           end if

           !! Limit is half original h
           if(ii.gt.log(0.6)/log(reduction_factor)) then
              hchecksum = two*res_tol/hh !! Limit total reduction to XX%.
write(6,*) i,i1,"stopping because of max reduction limit",ii
           end if


!hchecksum = 2.0*res_tol/hh
           !! Check the h-reduction criteria
    
!#if order==4
!           hchecksum = two*res_tol/hh
!#endif    
           if(hchecksum.ge.res_tol/hh)then   !! breakout of do-while
              reduce_h=.false.
!write(6,*) i,"stopping due to residual",ii,hh,hchecksum,res_tol/hh 
           else  !! continue reducing h, copy tmp neighbour lists to main neighbour lists, and set h(i)
              h(i) = hh
              ii = ii + 1         !! Counter for number of times reduced
              ij_count(i)=0
              ij_link(i,k)=0
              do k=1,full_j_count_i
                 j=full_j_link_i(k)
                 ij_count(i) = ij_count(i) + 1
                 ij_link(i,ij_count(i)) = j
              end do        
           end if
                
        end do
        !! Temporary: store size in neighbourcountreal(i)
        neighbourcountreal(i) = dble(ij_count(i))
     end do
!     !$OMP END PARALLEL DO    
     deallocate(amat,mmat,bvecL,bvecX,bvecY,gvec,xvec)      
          
     
     !! The remainder of this subroutine is just for analysis...
     qq = zero
     !$OMP PARALLEL DO REDUCTION(+:qq)
     do i=1,npfb
        qq = qq + neighbourcountreal(i)**2
     end do
     !$OMP END PARALLEL DO
     qq = sqrt(qq/npfb)   
     
     !! Output to screen
     write(6,*) "ij_count mean,min:",qq,floor(minval(neighbourcountreal(1:npfb)))
               
     !! Deallocation
     deallocate(full_j_link_i,neighbourcountreal)
#endif
     return
  end subroutine adapt_stencils     
!! ------------------------------------------------------------------------------------------------
!! FUNCTIONS TO CALCULATE THE TAYLOR MONOMIALS
!! ------------------------------------------------------------------------------------------------
  function monomials12(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(13) :: cxvec
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,x10,y10,x11,y11,x12,y12
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y;x6=x5*x;y6=y5*y;x7=x6*x;y7=y6*y;x8=x7*x;y8=y7*y
     x9=x8*x;y9=y8*y;x10=x9*x;y10=y9*y;x11=x10*x;y11=y10*y;x12=x11*x;y12=y11*y
     cxvec(1) = (1.0/479001600.0)*x12;cxvec(2)=(12.0/479001600.0)*x11*y;cxvec(3)=(66.0/479001600.0)*x10*y2
     cxvec(4)=(220.0/479001600.0)*x9*y3;cxvec(5)=(495.0/479001600.0)*x8*y4;cxvec(6) = (792.0/479001600.0)*x7*y5
     cxvec(7)=(924.0/479001600.0)*x6*y6;cxvec(8)=(792.0/479001600.0)*x5*y7;cxvec(9)=(495.0/479001600.0)*x4*y8
     cxvec(10)=(220.0/479001600.0)*x3*y9;cxvec(11)=(66.0/479001600.0)*x2*y10;cxvec(12)=(12.0/479001600.0)*x*y11
     cxvec(13)=(1.0/479001600.0)*y12
  end function monomials12
  function monomials11(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(12) :: cxvec
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,x10,y10,x11,y11
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y;x6=x5*x;y6=y5*y;x7=x6*x;y7=y6*y;x8=x7*x;y8=y7*y
     x9=x8*x;y9=y8*y;x10=x9*x;y10=y9*y;x11=x10*x;y11=y10*y
     cxvec(1) = (1.0/39916800.0)*x11;cxvec(2)=(11.0/39916800.0)*x10*y;cxvec(3)=(55.0/39916800.0)*x9*y2
     cxvec(4)=(165.0/39916800.0)*x8*y3;cxvec(5)=(330.0/39916800.0)*x7*y4;cxvec(6) = (462.0/39916800.0)*x6*y5
     cxvec(7)=(462.0/39916800.0)*x5*y6;cxvec(8)=(330.0/39916800.0)*x4*y7;cxvec(9)=(165.0/39916800.0)*x3*y8
     cxvec(10)=(55.0/39916800.0)*x2*y9;cxvec(11)=(11.0/39916800.0)*x*y10;cxvec(12)=(1.0/39916800.0)*y11
  end function monomials11
  function monomials10(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(11) :: cxvec
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,x10,y10
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y;x6=x5*x;y6=y5*y;x7=x6*x;y7=y6*y;x8=x7*x;y8=y7*y
     x9=x8*x;y9=y8*y;x10=x9*x;y10=y9*y
     cxvec(1) = (1.0/3628800.0)*x10;cxvec(2)=(1.0/362880.0)*x9*y;cxvec(3)=(1.0/80640.0)*x8*y2
     cxvec(4)=(1.0/30240.0)*x7*y3;cxvec(5)=(1.0/17280.0)*x6*y4;cxvec(6) = (1.0/14400.0)*x5*y5
     cxvec(7)=(1.0/17280.0)*x4*y6;cxvec(8)=(1.0/30240.0)*x3*y7;cxvec(9)=(1.0/80640.0)*x2*y8
     cxvec(10) = (1.0/362880.0)*x*y9;cxvec(11) = (1.0/3628800.0)*y10
  end function monomials10
  function monomials9(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(10) :: cxvec
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y;x6=x5*x;y6=y5*y;x7=x6*x;y7=y6*y;x8=x7*x;y8=y7*y
     x9=x8*x;y9=y8*y
     cxvec(1) = (1.0/362880.0)*x9;cxvec(2)=(1.0/40320.0)*x8*y;cxvec(3)=(1.0/10080.0)*x7*y2
     cxvec(4)=(1.0/4320.0)*x6*y3;cxvec(5)=(1.0/2880.0)*x5*y4
     cxvec(6) = (1.0/2880.0)*x4*y5;cxvec(7)=(1.0/4320.0)*x3*y6
     cxvec(8)=(1.0/10080.0)*x2*y7;cxvec(9)=(1.0/40320.0)*x*y8;cxvec(10) = (1.0/362880.0)*y9
  end function monomials9
  function monomials8(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(9) :: cxvec
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y;x6=x5*x;y6=y5*y;x7=x6*x;y7=y6*y;x8=x7*x;y8=y7*y         
     cxvec(1) = (1.0/40320.0)*x8;cxvec(2)=(1.0/5040.0)*x7*y;cxvec(3)=(1.0/1440.0)*x6*y2;cxvec(4)=(1.0/720.0)*x5*y3;
     cxvec(5) = (1.0/576.0)*x4*y4;cxvec(6)=(1.0/720.0)*x3*y5;cxvec(7)=(1.0/1440.0)*x2*y6;cxvec(8)=(1.0/5040.0)*x*y7
     cxvec(9) = (1.0/40320.0)*y8
  end function monomials8
  function monomials7(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(8) :: cxvec
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y;x6=x5*x;y6=y5*y;x7=x6*x;y7=y6*y         
     cxvec(1) = (1.0/5040.0)*x7;cxvec(2)=(1.0/720.0)*x6*y;cxvec(3)=(1.0/240.0)*x5*y2;cxvec(4)=(1.0/144.0)*x4*y3;
     cxvec(5) = (1.0/144.0)*x3*y4;cxvec(6)=(1.0/240.0)*x2*y5;cxvec(7)=(1.0/720.0)*x*y6;cxvec(8)=(1.0/5040.0)*y7
  end function monomials7
  function monomials6(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(7) :: cxvec
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y;x6=x5*x;y6=y5*y         
     cxvec(1) = (1.0/720.0)*x6;cxvec(2)=(1.0/120.0)*x5*y;cxvec(3)=(1.0/48.0)*x4*y2;cxvec(4)=(1.0/36.0)*x3*y3;
     cxvec(5) = (1.0/48.0)*x2*y4;cxvec(6)=(1.0/120.0)*x*y5;cxvec(7)=(1.0/720.0)*y6
  end function monomials6
  function monomials5(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(6) :: cxvec
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y         
     cxvec(1) = (1.0/120.0)*x5;cxvec(2)=(1.0/24.0)*x4*y;cxvec(3)=(1.0/12.0)*x3*y2
     cxvec(4) = (1.0/12.0)*x2*y3;cxvec(5)=(1.0/24.0)*x*y4;cxvec(6)=(1.0/120.0)*y5
  end function monomials5
  function monomials4(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(5) :: cxvec
     real(rkind) :: x2,y2,x3,y3,x4,y4
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y
     cxvec(1) = (1.0/24.0)*x4;cxvec(2)=(1.0/6.0)*x3*y;
     cxvec(3) = 0.25d0*x2*y2;cxvec(4)=(1.0/6.0)*x*y3;cxvec(5)=(1.0/24.0)*y4
  end function monomials4
  function monomials3(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(4) :: cxvec
     real(rkind) :: x2,y2,x3,y3
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y 
     cxvec(1) = (1.0/6.0)*x3;cxvec(2) = 0.5d0*x2*y;cxvec(3) = 0.5d0*x*y2;cxvec(4) = (1.0/6.0)*y3;
  end function monomials3
  function monomials2(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(5) :: cxvec
     real(rkind) :: x2,y2
     x2=x*x;y2=y*y
     cxvec(1) = x;cxvec(2)=y;cxvec(3)=0.5d0*x2;cxvec(4)=x*y;cxvec(5)=0.5d0*y2
  end function monomials2
!! ------------------------------------------------------------------------------------------------
!! ------------------------------------------------------------------------------------------------
!! ------------------------------------------------------------------------------------------------
#if ABF==1
!! ------------------------------------------------------------------------------------------------
!! ABFs generated from partial derivatives of an RBF
!! The "original" LABFM
!! ------------------------------------------------------------------------------------------------
  function abfs8(rad,x,y) result(ggvec)   !! EIGHT
     real(rkind),intent(in) :: rad,x,y
     real(rkind),dimension(9) :: ggvec
     real(rkind) :: rad3,rad2,r15
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8
     rad2 = max(rad*rad,eta2);rad3=max(rad*rad*rad,eta3);r15 = rad3*rad3*rad3*rad3*rad3
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y;x6=x5*x;y6=y5*y;x7=x6*x;y7=y6*y
     x8=x7*x;y8=y7*y
     ggvec(1) = (20160.0*x6*y2 - 75600.0*x4*y4+37800.0*x2*y6-1575.0*y8)/r15        
     ggvec(2) = (11025.0*x*y7-66150.0*x3*y5+52920.0*x5*y3-5040.0*x7*y)/r15
     ggvec(3) = (720.0*x8-24840.0*x6*y2+74250.0*x4*y4-33975.0*x2*y6+1350.0*y8)/r15
     ggvec(4) = (61425.0*x3*y5-9450.0*x*y7-56700*x5*y3-7560.0*x7*y)/r15
     ggvec(5) = (29700.0*(x2*y6+x6*y2)-1080.0*(x8+y8)-73575.0*x4*y4)/r15
     ggvec(6) = (61425.0*x5*y3-9450.0*x7*y-56700*x3*y5-7560.0*x*y7)/r15
     ggvec(7) = (720.0*y8-24840.0*x2*y6+74250.0*x4*y4-33975.0*x6*y2+1350.0*x8)/r15
     ggvec(8) = (11025.0*x7*y-66150.0*x5*y3+52920.0*x3*y5-5040.0*x*y7)/r15
     ggvec(9) = (20160.0*x2*y6 - 75600.0*x4*y4+37800.0*x6*y2-1575.0*x8)/r15 
  end function abfs8
  function abfs7(rad,x,y) result(ggvec)     !! SEVEN
     real(rkind),intent(in) :: rad,x,y
     real(rkind),dimension(8) :: ggvec
     real(rkind) :: rad3,rad2,r13
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7
     rad2 = max(rad*rad,eta2);rad3=max(rad*rad*rad,eta3);r13 = rad3*rad3*rad3*rad2*rad2
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y;x6=x5*x;y6=y5*y;x7=x6*x;y7=y6*y
     ggvec(1) = (6300.0*x3*y4 - 2520.0*x5*y2 - 1575.0*x*y6)/r13               
     ggvec(2) = (720.0*x6*y - 225.0*y7 + 4050.0*x2*y5 - 5400.0*x4*y3)/r13
     ggvec(3) = (1350.0*x*y6 - 120.0*x7 + 3000.0*x5*y2 + 5925.0*x3*y4)/r13
     ggvec(4) = (180.0*y7 - 3555.0*x2*y5 + 5580.0*x4*y3 - 1080.0*x6*y)/r13
     ggvec(5) = (180.0*x7 - 3555.0*x5*y2 + 5580.0*x3*y4 - 1080.0*x*y6)/r13
     ggvec(6) = (1350.0*x6*y - 120.0*y7 + 3000.0*x2*y5 + 5925.0*x4*y3)/r13
     ggvec(7) = (720.0*x*y6 - 225.0*x7 + 4050.0*x5*y2 - 5400.0*x3*y4)/r13
     ggvec(8) = (6300.0*x4*y3 - 2520.0*x2*y5 - 1575.0*x6*y)/r13 
  end function abfs7
  function abfs6(rad,x,y) result(ggvec)   !! SIX
     real(rkind),intent(in) :: rad,x,y
     real(rkind),dimension(7) :: ggvec
     real(rkind) :: rad3,rad2,r11
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6
     rad2 = max(rad*rad,eta2);rad3=max(rad*rad*rad,eta3);r11 = rad3*rad3*rad3*rad2
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y;x6=x5*x;y6=y5*y
     ggvec(1) = (360.0*x4*y2 - 540.0*x2*y4 + 45.0*y6)/r11                  
     ggvec(2) = (600.0*x3*y3 - 225.0*x*y5 - 120.0*x5*y)/r11
     ggvec(3) = (477.0*x2*y4 - 408.0*x4*y2 - 36.0*y6 + 24.0*x6)/r11
     ggvec(4) = (180.0*x*y5 - 585.0*x3*y3 + 180.0*x5*y)/r11
     ggvec(5) = (477.0*x4*y2 - 408.0*x2*y4 - 36.0*x6 + 24.0*y6)/r11
     ggvec(6) = (600.0*x3*y3 - 225.0*x5*y - 120.0*x*y5)/r11
     ggvec(7) = (360.0*x2*y4 - 540.0*x4*y2 + 45.0*x6)/r11
  end function abfs6
  function abfs5(rad,x,y) result(ggvec)      !! FIVE
     real(rkind),intent(in) :: rad,x,y
     real(rkind),dimension(6) :: ggvec
     real(rkind) :: rad3,r9
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5
     rad3=max(rad*rad*rad,eta3);r9=rad3*rad3*rad3
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y
     ggvec(1) = (45.0*x*y4 - 60.0*x3*y2)/r9                    
     ggvec(2) = (9.0*y5 - 72.0*x2*y3 + 24.0*x4*y)/r9
     ggvec(3) = (63.0*x3*y2 - 36.0*x*y4 - 6.0*x5)/r9
     ggvec(4) = (63.0*x2*y3 - 36.0*x4*y - 6.0*y5)/r9
     ggvec(5) = (9.0*x5 - 72.0*x3*y2 + 24.0*x*y4)/r9
     ggvec(6) = (45.0*x4*y - 60.0*x2*y3)/r9
  end function abfs5
  function abfs4(rad,x,y) result(ggvec)     !! FOUR
     real(rkind),intent(in) :: rad,x,y
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: rad3,rad2,r7
     real(rkind) :: x2,y2,x3,y3,x4,y4
     rad2 = max(rad*rad,eta2);rad3=max(rad*rad*rad,eta3);r7=rad2*rad2*rad3
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y
     ggvec(1) = (12.0*x2*y2-3.0*y4)/r7                
     ggvec(2) = (9.0*x*y3-6.0*x3*y)/r7
     ggvec(3) = (2.0*x4-11.0*x2*y2+2.0*y4)/r7
     ggvec(4) = (9.0*x3*y-6.0*x*y3)/r7
     ggvec(5) = (12.0*x2*y2-3.0*x4)/r7
  end function abfs4
  function abfs3(rad,x,y) result(ggvec)     !! THREE
     real(rkind),intent(in) :: rad,x,y
     real(rkind),dimension(4) :: ggvec
     real(rkind) :: rad3,rad2,r5
     real(rkind) :: x2,y2,x3,y3
     rad2 = max(rad*rad,eta2);rad3=max(rad*rad*rad,eta3);r5=rad3*rad2
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y
     ggvec(1) = -3.0*y2*x/r5                 
     ggvec(2) = (2.0*x2*y - y3)/r5
     ggvec(3) = (2.0*x*y2 - x3)/r5
     ggvec(4) = -3.0*x2*y/r5
  end function abfs3
  function abfs2(rad,x,y) result(ggvec)     !! TWO AND ONE
     real(rkind),intent(in) :: rad,x,y
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: rad3,rad2
     real(rkind) :: x2,y2
     rad2 = max(rad*rad,eta2);rad3=max(rad*rad*rad,eta3)
     x2=x*x;y2=y*y 
     ggvec(1) = x/max(rad,eta2) 
     ggvec(2) = y/max(rad,eta2)
     ggvec(3) = y2/rad3    
     ggvec(4) = -x*y/rad3
     ggvec(5) = x2/rad3
  end function abfs2
!! ------------------------------------------------------------------------------------------------
#elif ABF==2
!! HERMITE ABFs: Bivariate Hermite polynomials (probabalistic kind)
!! Formula based on Area, Dimitrov & Godoy (2015), J. Math. Anal. Appl. 421(1):830-841
!! with a=c=1, b=0.
!! Generated with: H_{p,q}(x,y) = H_{p}(x/sqrt(2))*H_{q}(y/sqrt(2))*2^((p+q)/2)
!!
!! NB sqrt2 and oosqrt are set in kind_parameters
!! NB ff1 multiplies the Hermite polynomial by an RBF to improve accuracy
!! ------------------------------------------------------------------------------------------------
  function abfs12(dummy,x,y) result(ggvec)         !! TWELVE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(13) :: ggvec
     write(6,*) "WARNING, ONLY CODED TO 10th ORDER. STOPPING"
     stop
              
  end function abfs12
  function abfs11(dummy,x,y) result(ggvec)         !! ELEVEN
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(12) :: ggvec
     write(6,*) "WARNING, ONLY CODED TO 10th ORDER. STOPPING"
     stop
  end function abfs11
  function abfs10(dummy,x,y) result(ggvec)         !! TEN
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(11) :: ggvec
     ggvec(1) = Hermite10(x)
     ggvec(2) = Hermite9(x)*Hermite1(y)
     ggvec(3) = Hermite8(x)*Hermite2(y)
     ggvec(4) = Hermite7(x)*Hermite3(y)
     ggvec(5) = Hermite6(x)*Hermite4(y)
     ggvec(6) = Hermite5(x)*Hermite5(y)
     ggvec(7) = Hermite4(x)*Hermite6(y)
     ggvec(8) = Hermite3(x)*Hermite7(y)
     ggvec(9) = Hermite2(x)*Hermite8(y)
     ggvec(10)= Hermite1(x)*Hermite9(y)
     ggvec(11)= Hermite10(y)
  end function abfs10
  function abfs9(dummy,x,y) result(ggvec)       !! NINE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(10) :: ggvec
     ggvec(1) = Hermite9(x)
     ggvec(2) = Hermite8(x)*Hermite1(y)
     ggvec(3) = Hermite7(x)*Hermite2(y)
     ggvec(4) = Hermite6(x)*Hermite3(y)
     ggvec(5) = Hermite5(x)*Hermite4(y)
     ggvec(6) = Hermite4(x)*Hermite5(y)
     ggvec(7) = Hermite3(x)*Hermite6(y)
     ggvec(8) = Hermite2(x)*Hermite7(y)
     ggvec(9) = Hermite1(x)*Hermite8(y)
     ggvec(10)= Hermite9(y)
  end function abfs9
  function abfs8(dummy,x,y) result(ggvec)            !! EIGHT
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(9) :: ggvec
     ggvec(1) = Hermite8(x)
     ggvec(2) = Hermite7(x)*Hermite1(y)
     ggvec(3) = Hermite6(x)*Hermite2(y)
     ggvec(4) = Hermite5(x)*Hermite3(y)
     ggvec(5) = Hermite4(x)*Hermite4(y)
     ggvec(6) = Hermite3(x)*Hermite5(y)
     ggvec(7) = Hermite2(x)*Hermite6(y)
     ggvec(8) = Hermite1(x)*Hermite7(y)
     ggvec(9) = Hermite8(y)
  end function abfs8
  function abfs7(dummy,x,y) result(ggvec)         !! SEVEN
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(8) :: ggvec
     ggvec(1) = Hermite7(x)
     ggvec(2) = Hermite6(x)*Hermite1(y)
     ggvec(3) = Hermite5(x)*Hermite2(y)
     ggvec(4) = Hermite4(x)*Hermite3(y)
     ggvec(5) = Hermite3(x)*Hermite4(y)
     ggvec(6) = Hermite2(x)*Hermite5(y)
     ggvec(7) = Hermite1(x)*Hermite6(y)
     ggvec(8) = Hermite7(y)
  end function abfs7
  function abfs6(dummy,x,y) result(ggvec)        !! SIX
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(7) :: ggvec
     ggvec(1) = Hermite6(x)
     ggvec(2) = Hermite5(x)*Hermite1(y)
     ggvec(3) = Hermite4(x)*Hermite2(y)
     ggvec(4) = Hermite3(x)*Hermite3(y)
     ggvec(5) = Hermite2(x)*Hermite4(y)
     ggvec(6) = Hermite1(x)*Hermite5(y)
     ggvec(7) = Hermite6(y)
  end function abfs6
  function abfs5(dummy,x,y) result(ggvec)      !! FIVE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(6) :: ggvec
     ggvec(1) = Hermite5(x)
     ggvec(2) = Hermite4(x)*Hermite1(y)
     ggvec(3) = Hermite3(x)*Hermite2(y)
     ggvec(4) = Hermite2(x)*Hermite3(y)
     ggvec(5) = Hermite1(x)*Hermite4(y)
     ggvec(6) = Hermite5(y)
  end function abfs5
  function abfs4(dummy,x,y) result(ggvec)      !! FOUR
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(5) :: ggvec
     ggvec(1) = Hermite4(x)
     ggvec(2) = Hermite3(x)*Hermite1(y)
     ggvec(3) = Hermite2(x)*Hermite2(y)
     ggvec(4) = Hermite1(x)*Hermite3(y)
     ggvec(5) = Hermite4(y)
  end function abfs4
  function abfs3(dummy,x,y) result(ggvec)     !!! THREE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(4) :: ggvec
     ggvec(1) = Hermite3(x)
     ggvec(2) = Hermite2(x)*Hermite1(y)
     ggvec(3) = Hermite1(x)*Hermite2(y)
     ggvec(4) = Hermite3(y)
  end function abfs3
  function abfs2(dummy,x,y) result(ggvec)     !! TWO AND ONE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(5) :: ggvec
     ggvec(1) = Hermite1(x)
     ggvec(2) = Hermite1(y)
     ggvec(3) = Hermite2(x)
     ggvec(4) = Hermite1(x)*Hermite1(y)
     ggvec(5) = Hermite2(y)
  end function abfs2
!! ------------------------------------------------------------------------------------------------
#elif ABF==3
!! LEGENDRE ABFs: Legendre polynomials
!! NB ff1 multiplies the Legendre polynomial by an RBF to improve accuracy
!! ------------------------------------------------------------------------------------------------
  function abfs10(dummy,x,y) result(ggvec)      !! TEN
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(11) :: ggvec
     ggvec(1) = Legendre10(x)
     ggvec(2) = Legendre9(x)*Legendre1(y)
     ggvec(3) = Legendre8(x)*Legendre2(y)
     ggvec(4) = Legendre7(x)*Legendre3(y)
     ggvec(5) = Legendre6(x)*Legendre4(y)
     ggvec(6) = Legendre5(x)*Legendre5(y)
     ggvec(7) = Legendre4(x)*Legendre6(y)
     ggvec(8) = Legendre3(x)*Legendre7(y)
     ggvec(9) = Legendre2(x)*Legendre8(y)
     ggvec(10)= Legendre1(x)*Legendre9(y)
     ggvec(11)= Legendre10(y)
  end function abfs10
  function abfs9(dummy,x,y) result(ggvec)      !! NINE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(10) :: ggvec
     ggvec(1) = Legendre9(x)
     ggvec(2) = Legendre8(x)*Legendre1(y)
     ggvec(3) = Legendre7(x)*Legendre2(y)
     ggvec(4) = Legendre6(x)*Legendre3(y)
     ggvec(5) = Legendre5(x)*Legendre4(y)
     ggvec(6) = Legendre4(x)*Legendre5(y)
     ggvec(7) = Legendre3(x)*Legendre6(y)
     ggvec(8) = Legendre2(x)*Legendre7(y)
     ggvec(9) = Legendre1(x)*Legendre8(y)
     ggvec(10)= Legendre9(y)
  end function abfs9
  function abfs8(dummy,x,y) result(ggvec)      !! EIGHT
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(9) :: ggvec
     ggvec(1) = Legendre8(x)
     ggvec(2) = Legendre7(x)*Legendre1(y)
     ggvec(3) = Legendre6(x)*Legendre2(y)
     ggvec(4) = Legendre5(x)*Legendre3(y)
     ggvec(5) = Legendre4(x)*Legendre4(y)
     ggvec(6) = Legendre3(x)*Legendre5(y)
     ggvec(7) = Legendre2(x)*Legendre6(y)
     ggvec(8) = Legendre1(x)*Legendre7(y)
     ggvec(9) = Legendre8(y)
  end function abfs8
  function abfs7(dummy,x,y) result(ggvec)      !! SEVEN
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(8) :: ggvec
     ggvec(1) = Legendre7(x)
     ggvec(2) = Legendre6(x)*Legendre1(y)
     ggvec(3) = Legendre5(x)*Legendre2(y)
     ggvec(4) = Legendre4(x)*Legendre3(y)
     ggvec(5) = Legendre3(x)*Legendre4(y)
     ggvec(6) = Legendre2(x)*Legendre5(y)
     ggvec(7) = Legendre1(x)*Legendre6(y)
     ggvec(8) = Legendre7(y)
  end function abfs7
  function abfs6(dummy,x,y) result(ggvec)      !! SIX
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(7) :: ggvec
     ggvec(1) = Legendre6(x)
     ggvec(2) = Legendre5(x)*Legendre1(y)
     ggvec(3) = Legendre4(x)*Legendre2(y)
     ggvec(4) = Legendre3(x)*Legendre3(y)
     ggvec(5) = Legendre2(x)*Legendre4(y)
     ggvec(6) = Legendre1(x)*Legendre5(y)
     ggvec(7) = Legendre6(y)
  end function abfs6
  function abfs5(dummy,x,y) result(ggvec)      !! FIVE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(6) :: ggvec
     ggvec(1) = Legendre5(x)
     ggvec(2) = Legendre4(x)*Legendre1(y)
     ggvec(3) = Legendre3(x)*Legendre2(y)
     ggvec(4) = Legendre2(x)*Legendre3(y)
     ggvec(5) = Legendre1(x)*Legendre4(y)
     ggvec(6) = Legendre5(y)
  end function abfs5
  function abfs4(dummy,x,y) result(ggvec)      !! FOUR
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(5) :: ggvec
     ggvec(1) = Legendre4(x)
     ggvec(2) = Legendre3(x)*Legendre1(y)
     ggvec(3) = Legendre2(x)*Legendre2(y)
     ggvec(4) = Legendre1(x)*Legendre3(y)
     ggvec(5) = Legendre4(y)
  end function abfs4
  function abfs3(dummy,x,y) result(ggvec)     !!! THREE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(4) :: ggvec
     ggvec(1) = Legendre3(x)
     ggvec(2) = Legendre2(x)*Legendre1(y)
     ggvec(3) = Legendre1(x)*Legendre2(y)
     ggvec(4) = Legendre3(y)
  end function abfs3
  function abfs2(dummy,x,y) result(ggvec)     !! TWO AND ONE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(5) :: ggvec
     ggvec(1) = Legendre1(x)
     ggvec(2) = Legendre1(y)
     ggvec(3) = Legendre2(x)
     ggvec(4) = Legendre1(x)*Legendre1(y)
     ggvec(5) = Legendre2(y)
  end function abfs2
!! ------------------------------------------------------------------------------------------------
#elif ABF==4  
!! LAGUERRE ABFs: Laguerre polynomials
!! NB ff1 multiplies the Laguerre polynomial by an RBF to improve accuracy
!! ------------------------------------------------------------------------------------------------
  function abfs7(dummy,x,y) result(ggvec)       
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(8) :: ggvec
     write(6,*) "WARNING, Laguerre ABFs only implemented up to 6th order so far. Stopping."
     stop
  end function abfs7  
  function abfs6(dummy,x,y) result(ggvec)        !! SIX
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(7) :: ggvec
     real(rkind) :: xx,yy
     xx = oosqrt2*x;yy=oosqrt2*y
     ggvec(1) = Laguerre6(xx)*0.125d0
     ggvec(2) = Laguerre5(xx)*Laguerre1(yy)*0.125d0
     ggvec(3) = Laguerre4(xx)*Laguerre2(yy)*0.125d0
     ggvec(4) = Laguerre3(xx)*Laguerre3(yy)*0.125d0
     ggvec(5) = Laguerre2(xx)*Laguerre4(yy)*0.125d0
     ggvec(6) = Laguerre1(xx)*Laguerre5(yy)*0.125d0
     ggvec(7) = Laguerre6(yy)*0.125d0
  end function abfs6
  function abfs5(dummy,x,y) result(ggvec)      !! FIVE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(6) :: ggvec
     real(rkind) :: xx,yy
     xx = oosqrt2*x;yy=oosqrt2*y
     ggvec(1) = Laguerre5(xx)*oosqrt2*0.25d0
     ggvec(2) = Laguerre4(xx)*Laguerre1(yy)*oosqrt2*0.25d0
     ggvec(3) = Laguerre3(xx)*Laguerre2(yy)*oosqrt2*0.25d0
     ggvec(4) = Laguerre2(xx)*Laguerre3(yy)*oosqrt2*0.25d0
     ggvec(5) = Laguerre1(xx)*Laguerre4(yy)*oosqrt2*0.25d0
     ggvec(6) = Laguerre5(yy)*oosqrt2*0.25d0
  end function abfs5
  function abfs4(dummy,x,y) result(ggvec)      !! FOUR
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: xx,yy
     xx = oosqrt2*x;yy=oosqrt2*y
     ggvec(1) = Laguerre4(xx)*0.25d0
     ggvec(2) = Laguerre3(xx)*Laguerre1(yy)*0.25d0
     ggvec(3) = Laguerre2(xx)*Laguerre2(yy)*0.25d0
     ggvec(4) = Laguerre1(xx)*Laguerre3(yy)*0.25d0
     ggvec(5) = Laguerre4(yy)*0.25d0
  end function abfs4
  function abfs3(dummy,x,y) result(ggvec)     !!! THREE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(4) :: ggvec
     real(rkind) :: xx,yy
     xx = oosqrt2*x;yy=oosqrt2*y
     ggvec(1) = Laguerre3(xx)*oosqrt2*0.5d0
     ggvec(2) = Laguerre2(xx)*Laguerre1(yy)*oosqrt2*0.5d0
     ggvec(3) = Laguerre1(xx)*Laguerre2(yy)*oosqrt2*0.5d0
     ggvec(4) = Laguerre3(yy)*oosqrt2*0.5d0
  end function abfs3
  function abfs2(dummy,x,y) result(ggvec)     !! TWO AND ONE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: xx,yy
     xx = x*oosqrt2;yy=y*oosqrt2
     ggvec(1) = Laguerre1(xx)*oosqrt2
     ggvec(2) = Laguerre1(yy)*oosqrt2
     ggvec(3) = Laguerre2(xx)*0.5d0
     ggvec(4) = Laguerre1(xx)*Laguerre1(yy)*0.5d0
     ggvec(5) = Laguerre2(yy)*0.5d0
  end function abfs2
!! ------------------------------------------------------------------------------------------------  
#elif ABF==5  
  function abfs7(dummy,x,y) result(ggvec)
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(8) :: ggvec
     write(6,*) "WARNING, Taylor ABFs only implemented up to 6th order so far. Stopping."
     stop     
  end function abfs7
  function abfs6(dummy,x,y) result(ggvec)
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(7) :: ggvec
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y;x6=x5*x;y6=y5*y         
     ggvec(1) = ff1*(1.0/720.0)*x6
     ggvec(2) = ff1*(1.0/120.0)*x5*y
     ggvec(3) = ff1*(1.0/48.0)*x4*y2
     ggvec(4) = ff1*(1.0/36.0)*x3*y3
     ggvec(5) = ff1*(1.0/48.0)*x2*y4
     ggvec(6) = ff1*(1.0/120.0)*x*y5
     ggvec(7) = ff1*(1.0/720.0)*y6
  end function abfs6
  function abfs5(dummy,x,y) result(ggvec)
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(6) :: ggvec
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y;x5=x4*x;y5=y4*y         
     ggvec(1) = ff1*(1.0/120.0)*x5
     ggvec(2) = ff1*(1.0/24.0)*x4*y
     ggvec(3) = ff1*(1.0/12.0)*x3*y2
     ggvec(4) = ff1*(1.0/12.0)*x2*y3
     ggvec(5) = ff1*(1.0/24.0)*x*y4
     ggvec(6) = ff1*(1.0/120.0)*y5
  end function abfs5
  function abfs4(dummy,x,y) result(ggvec)
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: x2,y2,x3,y3,x4,y4
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y
     ggvec(1) = ff1*(1.0/24.0)*x4
     ggvec(2) = ff1*(1.0/6.0)*x3*y
     ggvec(3) = ff1*0.25d0*x2*y2
     ggvec(4) = ff1*(1.0/6.0)*x*y3
     ggvec(5) = ff1*(1.0/24.0)*y4
  end function abfs4
  function abfs3(dummy,x,y) result(ggvec)
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(4) :: ggvec
     real(rkind) :: x2,y2,x3,y3
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y 
     ggvec(1) = ff1*(1.0/6.0)*x3
     ggvec(2) = ff1*0.5d0*x2*y
     ggvec(3) = ff1*0.5d0*x*y2
     ggvec(4) = ff1*(1.0/6.0)*y3
  end function abfs3
  function abfs2(dummy,x,y) result(ggvec)
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: x2,y2
     x2=x*x;y2=y*y
     ggvec(1) = ff1*x
     ggvec(2) = ff1*y
     ggvec(3) = ff1*0.5d0*x2
     ggvec(4) = ff1*x*y
     ggvec(5) = ff1*0.5d0*y2
  end function abfs2
#elif ABF==6
!! Henry Broadley Fourier Basis
!! ------------------------------------------------------------------------------------------------
  function abfs10(dummy,x,y) result(ggvec)         !! TEN
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(11) :: ggvec
     ggvec(1) = cos(5.0*x)
     ggvec(2) = sin(5.0*x)*sin(y)
     ggvec(3) = cos(4.0*x)*cos(y)
     ggvec(4) = sin(4.0*x)*sin(2.0*y)
     ggvec(5) = cos(3.0*x)*cos(2.0*y)
     ggvec(6) = sin(3.0*x)*sin(3.0*y)
     ggvec(7) = cos(2.0*x)*cos(3.0*y)
     ggvec(8) = sin(2.0*x)*sin(4.0*y)
     ggvec(9) = cos(x)*cos(4.0*y)
     ggvec(10)= sin(x)*sin(5.0*y)
     ggvec(11)= cos(5.0*y)
  end function abfs10
  function abfs9(dummy,x,y) result(ggvec)         !! NINE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(10) :: ggvec
     ggvec(1) = sin(5.0*x)
     ggvec(2) = cos(4.0*x)*sin(y)
     ggvec(3) = sin(4.0*x)*cos(y)
     ggvec(4) = cos(3.0*x)*sin(2.0*y)
     ggvec(5) = sin(3.0*x)*cos(2.0*y)
     ggvec(6) = cos(2.0*x)*sin(3.0*y)
     ggvec(7) = sin(2.0*x)*cos(3.0*y)
     ggvec(8) = cos(x)*sin(4.0*y)
     ggvec(9) = sin(x)*cos(4.0*y)
     ggvec(10)= sin(5.0*y)
  end function abfs9
  function abfs8(dummy,x,y) result(ggvec)         !! EIGHT
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(9) :: ggvec
     ggvec(1) = cos(4.0*x)
     ggvec(2) = sin(4.0*x)*sin(y)
     ggvec(3) = cos(3.0*x)*cos(1.0*y)
     ggvec(4) = sin(3.0*x)*sin(2.0*y)
     ggvec(5) = cos(2.0*x)*cos(2.0*y)
     ggvec(6) = sin(2.0*x)*sin(3.0*y)
     ggvec(7) = cos(1.0*x)*cos(3.0*y)
     ggvec(8) = sin(x)*sin(4.0*y)
     ggvec(9) = cos(4.0*y)
  end function abfs8
  function abfs7(dummy,x,y) result(ggvec)         !! SEVEN
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(8) :: ggvec
     ggvec(1) = sin(4.0*x)
     ggvec(2) = cos(3.0*x)*sin(y)
     ggvec(3) = sin(3.0*x)*cos(1.0*y)
     ggvec(4) = cos(2.0*x)*sin(2.0*y)
     ggvec(5) = sin(2.0*x)*cos(2.0*y)
     ggvec(6) = cos(1.0*x)*sin(3.0*y)
     ggvec(7) = sin(x)*cos(3.0*y)
     ggvec(8) = sin(4.0*y)
  end function abfs7
  function abfs6(dummy,x,y) result(ggvec)        !! SIX
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(7) :: ggvec
     ggvec(1) = cos(3.0*x)
     ggvec(2) = sin(3.0*x)*sin(y)
     ggvec(3) = cos(2.0*x)*cos(1.0*y)
     ggvec(4) = sin(2.0*x)*sin(2.0*y)
     ggvec(5) = cos(1.0*x)*cos(2.0*y)
     ggvec(6) = sin(x)*sin(3.0*y)
     ggvec(7) = cos(3.0*y)
  end function abfs6
  function abfs5(dummy,x,y) result(ggvec)      !! FIVE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(6) :: ggvec
     ggvec(1) = sin(3.0*x)
     ggvec(2) = cos(2.0*x)*sin(y)
     ggvec(3) = sin(2.0*x)*cos(1.0*y)
     ggvec(4) = cos(1.0*x)*sin(2.0*y)
     ggvec(5) = sin(x)*cos(2.0*y)
     ggvec(6) = sin(3.0*y)
  end function abfs5
  function abfs4(dummy,x,y) result(ggvec)      !! FOUR
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(5) :: ggvec
     ggvec(1) = cos(2.0*x)
     ggvec(2) = sin(2.0*x)*sin(y)
     ggvec(3) = cos(x)*cos(y)
     ggvec(4) = sin(x)*sin(2.0*y)
     ggvec(5) = cos(2.0*y)         
  end function abfs4
  function abfs3(dummy,x,y) result(ggvec)     !!! THREE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(4) :: ggvec
     ggvec(1) = sin(2.0*x)
     ggvec(2) = cos(x)*sin(y)
     ggvec(3) = sin(x)*cos(y)
     ggvec(4) = sin(2.0*y)
  end function abfs3
  function abfs2(dummy,x,y) result(ggvec)     !! TWO AND ONE
     real(rkind),intent(in) :: x,y,dummy
     real(rkind),dimension(5) :: ggvec
     ggvec(1) = sin(x)
     ggvec(2) = sin(y)
     ggvec(3) = cos(x)
     ggvec(4) = sin(x)*sin(y)
     ggvec(5) = cos(y)
  end function abfs2  
#endif
!! ------------------------------------------------------------------------------------------------
!! Below this line are functions which return univariate polynomials from which some ABFs above
!! are constructed
!! ------------------------------------------------------------------------------------------------
  !! Univariate Hermite Polynomials (Physicists kind)
  function Hermite1(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = z
  end function Hermite1
  function Hermite2(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = z*z - 1.0d0
  end function Hermite2
  function Hermite3(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = z*z*z - 3.0d0*z
  end function Hermite3
  function Hermite4(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = z*z*z*z - 6.0d0*z*z + 3.0d0
  end function Hermite4
  function Hermite5(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z3
     real(rkind) :: Hres
     z2=z*z;z3=z*z2
     Hres = z2*z3 - 10.0d0*z3 + 15.0*z
  end function Hermite5
  function Hermite6(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z3
     real(rkind) :: Hres
     z2=z*z;z3=z*z2
     Hres = z3*z3 - 15.0d0*z2*z2 + 45.0d0*z2 - 15.0d0
  end function Hermite6
  function Hermite7(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z3,z4
     real(rkind) :: Hres
     z2=z*z;z3=z*z2;z4=z3*z
     Hres = z4*z3 - 21.0d0*z3*z2 + 105.0d0*z3 - 105.0d0*z
  end function Hermite7
  function Hermite8(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z4
     real(rkind) :: Hres
     z2=z*z;z4=z2*z2
     Hres = z4*z4 - 28.0d0*z4*z2 + 210.0d0*z4 - 420.0d0*z2 + 105.0d0
  end function Hermite8
  function Hermite9(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z3,z4
     real(rkind) :: Hres
     z2=z*z;z3=z*z2;z4=z3*z
     Hres = z4*z3*z2 - 36.0d0*z4*z3 + 378.0d0*z3*z2 - 1260.0d0*z3 + 945.0d0*z
  end function Hermite9
  function Hermite10(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z4
     real(rkind) :: Hres
     z2=z*z;z4=z2*z2
     Hres = z4*z4*z2 - 45.0d0*z4*z4 + 630.0d0*z4*z2 - 3150.0d0*z4 &
            + 4725.0d0*z2 - 945.0d0
  end function Hermite10
!! ------------------------------------------------------------------------------------------------
!! Univariate Legendre polynomials
!! ------------------------------------------------------------------------------------------------
  function Legendre1(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = z
  end function Legendre1
  function Legendre2(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = 0.5d0*(3.0d0*z*z-1.0d0)
  end function Legendre2
  function Legendre3(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = 0.5d0*(5.0d0*z*z*z-3.0d0*z)
  end function Legendre3
  function Legendre4(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2
     z2=z*z
     Lres = 0.125d0*(35.0d0*z2*z2 - 30.0d0*z2 + 3.0d0)
  end function Legendre4
  function Legendre5(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z3
     z2=z*z;z3=z2*z
     Lres = 0.125d0*(63.0d0*z2*z3 - 70.0d0*z3 + 15.0d0*z)
  end function Legendre5
  function Legendre6(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z4
     z2=z*z;z4=z2*z2
     Lres = 0.0625d0*(231.0d0*z2*z4 - 316.0d0*z4 + 105.0d0*z2 - 5.0d0)
  end function Legendre6
  function Legendre7(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z3,z4
     z2=z*z;z3=z2*z;z4=z2*z2
     Lres = 0.0625d0*(429.0d0*z3*z4 - 693.0d0*z2*z3 + 315.0d0*z3 - 35.0d0*z)
  end function Legendre7
  function Legendre8(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z4
     z2=z*z;z4=z2*z2
     Lres = 7.8125d-3*(6435.0d0*z4*z4 - 12012.0d0*z4*z2 + 6930.0d0*z4 - 1260.0d0*z + 35.0d0)
  end function Legendre8
  function Legendre9(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z4,z3
     z2=z*z;z4=z2*z2;z3=z2*z
     Lres = 7.8125d-3*(12155.0d0*z3*z2*z4 - 25740.0d0*z4*z3 + 18018.0d0*z3*z2 &
          - 4620.0d0*z3 + 315.0d0*z)
  end function Legendre9
  function Legendre10(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z4
     z2=z*z;z4=z2*z2
     Lres = 3.90625d-3*(46189.0d0*z4*z4*z2 - 109395.0d0*z4*z4 + 90090.0d0*z4*z2 &
          - 30030.0d0*z4 + 3465.0d0*z2 - 63.0d0)
  end function Legendre10
!! ------------------------------------------------------------------------------------------------
!! Univariate Laguerre polynomials
!! ------------------------------------------------------------------------------------------------
  function Laguerre1(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = 1.0d0 - z
  end function Laguerre1
  function Laguerre2(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = 0.5d0*(z*z-4.0d0*z+2.0d0)
  end function Laguerre2
  function Laguerre3(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = (-z*z*z+9.0d0*z*z-18.0d0*z+6.0d0)/6.0d0
  end function Laguerre3
  function Laguerre4(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2
     z2=z*z
     Lres = (z2*z2 - 16.0d0*z2*z + 72.0d0*z2 - 96.0d0*z + 24.0d0)/24.0d0
  end function Laguerre4
  function Laguerre5(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z3
     z2=z*z;z3=z2*z
     Lres = (-z3*z2 +25.0d0*z2*z2 - 200.0d0*z3 + 600.0d0*z2 - 600.0d0*z + 120.0d0)/120.0d0
  end function Laguerre5
  function Laguerre6(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z4
     z2=z*z;z4=z2*z2
     Lres = (z4*z2 - 36.0d0*z4*z + 450.0d0*z4 - 2400.0d0*z2*z + 5400.0d0*z2 - 4320.0d0*z &
             + 720.0d0)/720.0d0
  end function Laguerre6
!! ------------------------------------------------------------------------------------------------
end module moments
