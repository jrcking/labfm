module nodes
!! This module contains routines to create node distributions, do iterative shifting
!! find neighbours, and output field variable files.
!! The three main distribution creation routines to use are:
!!    a) create_particles_banalytic - noisy cartesian distribution, with full support boundary
!!       and boundary prescribed with analytic function.
!!    b) create_particles_bperiodic - as a), but periodic boundaries. Periodic boundaries done via 
!!       parent-child linking (i.e. the inefficient way), as it's conceptually simpler when messing around
!!    c) create_particles_bperiodic_varh - as b), but with varying resolution.
  use kind_parameters
  use common_parameter
  use common_2d
  use analytic_functions
  use omp_lib
  implicit none
  
  integer(ikind),dimension(:),allocatable :: dontshift

!! There are npfb particles in the domain
!! 1 to nb_n are Neumann
!! nb_n + 1 to nb are Dirichlet
!! 1 to nb are "solid boundary"
!! nb+1 to npfb are "internal"
!! npfb + 1 to npfb + nb_n  are for satisfaction of Neumann bcs
!! npfb + nb_n + 1 to np are for other BCs (e.g. periodic, analytic...)

contains
  subroutine create_particles_banalytic
     !! Create particles in a Cartesian arrangement on a square domain (side length xmax-xmin)
     !! nx+1 particles per side length, including particles ON domain boundary.
     !! Create ceiling(sup_size/dx) layers of BC particles surrounding domain.
     !! Initialise variable u to an analytic function
     !! irelation is unallocated, to indicate some analytic function is needed for boundaries
     integer(ikind) :: i,tmp_i,tmp_j,j,imp,nss
     real(rkind) :: ns
     real(rkind),dimension(:,:),allocatable :: rptmp
     real(rkind),dimension(:),allocatable :: htmp

     time = 0.0d0
     dx = (xmax - xmin)/dble(nx);dv=dx*dx   ! set particle spacing based on particles per domain side length
     h0 = hovdx*dx;sup_size = ss*h0;h2=h0*h0;h3=h2*h0;ss2=sup_size*sup_size  !h, support 
     eta = 1.0d-8*h0;eta2 = eta*eta;eta3=eta*eta2
     hmin = h0

     nss = ceiling(sup_size/dx) + 1
     tmp_i = (nx+2*nss + 1)**2  !! np should equal this, if h=2dx. make the 11 bigger if bigger support... or 2*nss +1...
     npfb = (nx+1)**2    
     nb=0;nb_n=0 
     allocate(rp(tmp_i,dims))
     allocate(h(tmp_i));h=h0

     !! Create "fluid" particles
     call random_seed()
     !$OMP PARALLEL DO PRIVATE(tmp_i,tmp_j) 
     do i=1,npfb     !! Fluid particles in a square, nx+1 by nx+1,
           tmp_i=(i-1)/(nx+1)+1
           tmp_j=i - (tmp_i-1)*(nx+1)        
           call random_number(ns);ns = (ns - 0.5)*dx*tmp_noise
           rp(i,1) = xmin + dble(tmp_i-1)*dx + ns
           call random_number(ns);ns = (ns - 0.5)*dx*tmp_noise
           rp(i,2) = xmin + dble(tmp_j-1)*dx + ns
     end do
     !$OMP END PARALLEL DO

     !! Create boundary particles
     imp = npfb
     if(.true.)then
     do i=1,ceiling(sup_size/dx) + 1 - 1    ! Boundary particles  (-2 gives enough support..)
        do j=1,nx + 1 + ceiling(sup_size/dx) !+ 1
           imp = imp + 1;rp(imp,1) = xmin - dble(i)*dx;rp(imp,2) = xmin + dble(j-1)*dx
           imp = imp + 1;rp(imp,2) = xmax + dble(i)*dx;rp(imp,1) = xmin + dble(j-1)*dx
           imp = imp + 1;rp(imp,1) = xmax + dble(i)*dx;rp(imp,2) = xmax - dble(j-1)*dx
           imp = imp + 1;rp(imp,2) = xmin - dble(i)*dx;rp(imp,1) = xmax - dble(j-1)*dx
        end do
     end do
     end if
     np = imp 
    
     !! Set u(i) to a given differentiable function - this is set in sphtools
     allocate(u(np)) 
     allocate(v(np));allocate(w(np))
     !$OMP PARALLEL DO 
     do i=1,np  
        u(i) = ftn(rp(i,1),rp(i,2))
     end do
     !$OMP END PARALLEL DO
     

open(unit=30,file='fort.30')
do i=1,np
write(30,*) rp(i,:)
end do
flush(30);close(30)
     return
  end subroutine create_particles_banalytic
!! ------------------------------------------------------------------------------------------------
  subroutine create_particles_nosupport
     integer(ikind) :: i,tmp_i,tmp_j,j,imp,nss
     real(rkind) :: nsx,nsy,dr,rmin,rmax,tn,dbound

     time = 0.0d0
     dx = (xmax - xmin)/dble(nx);dv=dx*dx   ! set particle spacing based on particles per domain side length
     h0 = hovdx*dx;sup_size = ss*h0;h2=h0*h0;h3=h2*h0;ss2=sup_size*sup_size  !h, support 
     eta = 1.0d-8*h0;eta2 = eta*eta;eta3=eta*eta2
     hmin = h0

     nss = ceiling(sup_size/dx) + 1
     tmp_i = (nx*2 + 1)**2  !! np should equal this, if h=2dx. make the 11 bigger if bigger support... or 2*nss +1...
     npfb = (nx+1)**2    
     nb=0;nb_n=0 
     allocate(rp(tmp_i,dims))
     allocate(h(tmp_i));h=h0    
      

     !! Create Dirichlet particles
     dr = dx
     nss = floor((xmax-xmin)/dr)
     do i=1,nss
        nb=nb+1;rp(nb,1)=xmin+dble(i-1)*dr;rp(nb,2)=xmin
        nb=nb+1;rp(nb,1)=xmax-dble(i-1)*dr;rp(nb,2)=xmax
        nb=nb+1;rp(nb,2)=xmin+dble(i-1)*dr;rp(nb,1)=xmax
        nb=nb+1;rp(nb,2)=xmax-dble(i-1)*dr;rp(nb,1)=xmin
     end do
     npfb=nb
     nb_n=0

     !! Create a series of rings (squares) of nodes of decreasing size to fill the domain.
     !! tmp_i counts the rings
     call random_seed()
     rmin=xmin+dr;rmax=xmax-dr;tmp_i=0
     do while(rmin.le.rmax-dr)  !! Loop over rings of nodes...
        dbound = rmin-xmin
        if(dbound.le.h0) then
           dr = dx;tn = 0.0d0
        else
           dr = dx;tn = tmp_noise
        end if
        nss = nint((rmax-rmin)/dr)
        do i=1,nss
           call random_number(nsx);nsx=(nsx-0.5)*dr*tn;call random_number(nsy);nsy=(nsy-0.5)*dr*tn
           npfb=npfb+1;rp(npfb,1)=rmin+dble(i-1)*dr+nsx;rp(npfb,2)=rmin+nsy
           call random_number(nsx);nsx=(nsx-0.5)*dr*tn;call random_number(nsy);nsy=(nsy-0.5)*dr*tn
           npfb=npfb+1;rp(npfb,1)=rmax-dble(i-1)*dr+nsx;rp(npfb,2)=rmax+nsy
           call random_number(nsx);nsx=(nsx-0.5)*dr*tn;call random_number(nsy);nsy=(nsy-0.5)*dr*tn
           npfb=npfb+1;rp(npfb,2)=rmin+dble(i-1)*dr+nsy;rp(npfb,1)=rmax+nsx
           call random_number(nsx);nsx=(nsx-0.5)*dr*tn;call random_number(nsy);nsy=(nsy-0.5)*dr*tn
           npfb=npfb+1;rp(npfb,2)=rmax-dble(i-1)*dr+nsy;rp(npfb,1)=rmin+nsx
        end do
        rmin = rmin+dr;rmax=rmax-dr
     end do
     npfb=npfb+1;rp(npfb,:)=xmin + 0.5d0*(xmax-xmin)  !! Final node in centre...

     np = npfb
    
     !! Set u(i) to a given differentiable function - this is set in sphtools
     allocate(u(np)) 
     allocate(v(np));allocate(w(np))
     !$OMP PARALLEL DO 
     do i=1,np  
        u(i) = ftn(rp(i,1),rp(i,2))
     end do
     !$OMP END PARALLEL DO
         

open(unit=30,file='fort.30')
do i=1,np
write(30,*) rp(i,:)
end do
flush(30);close(30)

     return
  end subroutine create_particles_nosupport
!! ------------------------------------------------------------------------------------------------
  subroutine create_particles_bperiodic
     !! Create particles in a Cartesian arrangement on a square domain (side length xmax-xmin)
     !! nx+1 particles per side length, from dx/2,3dx/2,5dx/2 etc.
     !! Create ceiling(sup_size/dx) layers of BC particles surrounding domain.
     !! Initialise variable u to an analytic function
     !! BC have irelation(i) = j, with j mirror particle in domain, proving periodic BCs
     integer(ikind) :: i,tmp_i,tmp_j,j,imp,nss,nmirror
     real(rkind) :: ns
     integer(ikind) :: aspect_r

     aspect_r = 1
     ymax=ymax*dble(aspect_r);ymin=ymin*dble(aspect_r)

     time = 0.0d0
     dx = (xmax - xmin)/dble(nx+1);dv=dx*dx   ! set particle spacing based on particles per domain side length
     h0 = hovdx*dx;sup_size = ss*h0;h2=h0*h0;h3=h2*h0;ss2=sup_size*sup_size  !h, support 
     eta = 1.0d-8*h0;eta2 = eta*eta;eta3=eta*eta2
     hmin = h0   
     
     xbcond = 1;ybcond=1  

     nss = ceiling(sup_size/dx) + 1
     tmp_i = aspect_r*(nx+8*nss + 1)**2  !! np should equal this, if h=2dx. make the 11 bigger.. or 2*nss +1...
     npfb = aspect_r*(nx+1)**2     
     nb=0;nb_n=0
     allocate(rp(tmp_i,dims))
     allocate(h(tmp_i));h=h0

     !! Create "fluid" particles
     call random_seed()
     !$OMP PARALLEL DO PRIVATE(tmp_i,tmp_j,ns) 
     do i=1,npfb     !! Fluid particles in a square, nx+1 by nx+1,
           tmp_i=(i-1)/(nx+1)+1
           tmp_j=i - (tmp_i-1)*(nx+1)
           call random_number(ns);ns = (ns - 0.5)*dx*tmp_noise
           rp(i,2) = ymin + 0.5*dx + dble(tmp_i-1)*dx + ns
           call random_number(ns);ns = (ns - 0.5)*dx*tmp_noise
           rp(i,1) = xmin + 0.5*dx + dble(tmp_j-1)*dx + ns
     end do
     !$OMP END PARALLEL DO


     call create_mirror_particles            

     !write(6,*) minval(rp(1:np,1)),maxval(rp(1:np,1))
     !write(6,*) minval(rp(1:np,2)),maxval(rp(1:np,2))

     call iteratively_shift(100)
     call create_mirror_particles    
    
     !! Set u(i) to a given differentiable function - this is set in sphtools
     allocate(u(tmp_i)) 
     allocate(v(tmp_i));allocate(w(tmp_i),ro(tmp_i),Yspec(tmp_i),E(tmp_i))
     !$OMP PARALLEL DO 
     do i=1,npfb  
        u(i) = ftn(rp(i,1),rp(i,2))
     end do
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO PRIVATE(j)
     do i=npfb+1,np
        j=irelation(i)
        u(i) = u(j)
     end do
     !$OMP END PARALLEL DO

open(unit=30,file='fort.30')
do i=1,np
write(30,*) rp(i,:)
end do
flush(30);close(30)

     return
  end subroutine create_particles_bperiodic
!! ------------------------------------------------------------------------------------------------  
  subroutine create_particles_RT
     !! Create particles for Rayleigh-Taylor flow. L-R is periodic. U-D is symmetric
     integer(ikind) :: i,tmp_i,tmp_j,j,imp,nss,nmirror
     real(rkind) :: ns

     time = 0.0d0
     dx = (xmax - xmin)/dble(nx+1);dv=dx*dx   ! set particle spacing based on particles per domain side length
     h0 = hovdx*dx;sup_size = ss*h0;h2=h0*h0;h3=h2*h0;ss2=sup_size*sup_size  !h, support 
     eta = 1.0d-8*h0;eta2 = eta*eta;eta3=eta*eta2
     hmin = h0     
     
     xbcond =1;ybcond=2

     nss = ceiling(sup_size/dx) + 1
     tmp_i = (nx+8*nss + 1)**2  !! np should equal this, if h=2dx. make the 11 bigger if bigger support... or 2*nss +1...
     npfb = (nx+1)**2     
     nb=0;nb_n=0
     allocate(rp(tmp_i,dims))
     allocate(h(tmp_i));h=h0

     !! Create "fluid" particles
     call random_seed()
     !$OMP PARALLEL DO PRIVATE(tmp_i,tmp_j,ns) 
     do i=1,npfb     !! Fluid particles in a square, nx+1 by nx+1,
           tmp_i=(i-1)/(nx+1)+1
           tmp_j=i - (tmp_i-1)*(nx+1)
           call random_number(ns);ns = (ns - 0.5)*dx*tmp_noise
           rp(i,1) = xmin + 0.5*dx + dble(tmp_i-1)*dx + ns
           call random_number(ns);ns = (ns - 0.5)*dx*tmp_noise
           rp(i,2) = xmin + 0.5*dx + dble(tmp_j-1)*dx + ns
     end do
     !$OMP END PARALLEL DO

     call create_mirror_particles
     
     call iteratively_shift(100)
     call create_mirror_particles    
    
     !! Set u(i) to a given differentiable function - this is set in sphtools
     allocate(u(tmp_i)) 
     allocate(v(tmp_i));allocate(w(tmp_i),ro(tmp_i),Yspec(tmp_i),E(tmp_i))
     !$OMP PARALLEL DO 
     do i=1,npfb  
        u(i) = ftn(rp(i,1),rp(i,2))
     end do
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO PRIVATE(j)
     do i=npfb+1,np
        j=irelation(i)
        u(i) = u(j)
     end do
     !$OMP END PARALLEL DO

!open(unit=30,file='fort.30')
!do i=1,np
!write(30,*) rp(i,:)
!end do
!flush(30);close(30)
     return
  end subroutine create_particles_RT 
!! ------------------------------------------------------------------------------------------------
  subroutine create_particles_bperiodic_varh
     !! Create particles in a Cartesian arrangement on a square domain (side length xmax-xmin)
     !! nx+1 particles per side length, including particles ON domain boundary.
     !! Create ceiling(sup_size/dx) layers of BC particles surrounding domain.
     !! Initialise variable u to an analytic function
     !! irelation is unallocated, to indicate some analytic function is needed for boundaries
     integer(ikind) :: i,tmp_i,tmp_j,j,imp,nss
     real(rkind) :: ns,dx0,dy,x,y
     real(rkind),dimension(:,:),allocatable :: rptmp
     real(rkind),dimension(:),allocatable :: htmp

     nplink = 16*hovdx*hovdx*2 !! Make bigger for var res...
     time = 0.0d0
     dx = (xmax - xmin)/dble(nx);dv=dx*dx   ! set particle spacing based on particles per domain side length
     h0 = hovdx*dx;sup_size = ss*h0;h2=h0*h0;h3=h2*h0;ss2=sup_size*sup_size  !h, support 
     eta = 1.0d-8*h0;eta2 = eta*eta;eta3=eta*eta2
     hmin = h0/2.0d0

     nss = ceiling(sup_size/dx) + 1
     tmp_i = (4*nx+2*nss + 1)**2  !! np should equal this, if h=2dx. make the 11 bigger if bigger support... or 2*nss +1...
     nb=0;nb_n=0 
     allocate(rp(tmp_i,dims))
     allocate(h(tmp_i));h=h0

     !! Create "fluid" particles
     call random_seed()
     i=0;dx0=dx
     x=xmin+0.5d0*dx0
     do while (x.le.xmax)!-0.2d0*dx0)
        if(x.gt.0.15.and.x.le.0.35)then
           dx = dx0*(1.0d0+((h0/hmin)-1.0d0)*sin(pi*(x-0.25d0)/0.2)**2)/(h0/hmin)
        else if(x.gt.0.65.and.x.le.0.85)then
           dx = dx0*(1.0d0+((h0/hmin)-1.0d0)*sin(pi*(x-0.75d0)/0.2)**2)/(h0/hmin)
        else
           dx = dx0
        end if
        y=ymin+0.5d0*dx
        do while (y.le.ymax)!-0.2d0*dx)
           i=i+1;
           call random_number(ns);ns = (ns - 0.5)*dx*tmp_noise;rp(i,1)=x+ns
           call random_number(ns);ns = (ns - 0.5)*dx*tmp_noise;rp(i,2)=y+ns
           h(i) = dx*hovdx
           y=y+dx
        end do
        x=x+dx
     end do
     npfb = i;dx=dx0
     
     call create_mirror_particles         
     call iteratively_shift(50)
     call create_mirror_particles
     
         
     allocate(rptmp(np,dims),htmp(np))
     rptmp(1:np,:)=rp(1:np,:)
     htmp(1:np)=h(1:np)
!     deallocate(rp,h)
!     allocate(rp(np),hp
!     deallocate(rptmp,htmp)
    
     !! Set u(i) to a given differentiable function - this is set in sphtools
     allocate(u(np)) 
     allocate(v(np));allocate(w(np))
     !$OMP PARALLEL DO 
     do i=1,np  
        u(i) = ftn(rp(i,1),rp(i,2))
     end do
     !$OMP END PARALLEL DO

open(unit=30,file='fort.30')
do i=1,np
write(30,*) rp(i,:)
end do
flush(30);close(30)
     return
  end subroutine create_particles_bperiodic_varh
!! ------------------------------------------------------------------------------------------------
  subroutine output_uv(n_out)
     !! Little subroutine to write out field variables. 
     integer(ikind),intent(in) :: n_out
     integer(ikind) :: i,j,k
     character(70) :: fname

     !! set the name of the file...
     if( n_out .lt. 10 ) then 
        write(fname,'(A16,I1)') './data_out/uv/uv', n_out
     else if( n_out .lt. 100 ) then 
        write(fname,'(A16,I2)') './data_out/uv/uv', n_out
     else if( n_out .lt. 1000 ) then
        write(fname,'(A16,I3)') './data_out/uv/uv', n_out
     else
        write(fname,'(A16,I4)') './data_out/uv/uv', n_out
     end if

     open(unit = 20,file=fname)  
!     write(20,*) nx+1,dx,xmin,xmin,0      !! Write parameters


     do i=1,npfb
        write(20,*) rp(i,1),rp(i,2),h(i)/h0,u(i),v(i),ro(i),Yspec(i),E(i)
     end do
     flush(20)
     close(20)

     write(13,*) time,npfb,n_out,0.1
     flush(13)
     return
  end subroutine output_uv
!! ------------------------------------------------------------------------------------------------
  subroutine iteratively_shift(kk)
     use neighbours
     !! Subroutine to create a nice shifted distribution...
     integer(ikind),intent(in) :: kk
     integer(ikind) :: ll,i,j,k
     real(rkind),dimension(dims) :: rij,gradw,dr_tmp
     real(rkind),dimension(:,:),allocatable :: dr
     real(rkind) :: rad,tmp,qq
     real(rkind) :: qkd_mag,ns,drmag
    
     allocate(dr(npfb,dims))
  
     !! Low order shifting loop
     do ll=1,kk
        call create_mirror_particles       
        call find_neighbours
                

        !! Find shifting vector...
        !$OMP PARALLEL DO PRIVATE(k,j,rij,rad,qq,gradw,dr_tmp)
        do i=nb+1,npfb
           qkd_mag = 1.0d-1*h0!(i)
           dr_tmp = 0.0d0
           do k=1,ij_count(i)
              j=ij_link(i,k)
              rij = rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=hovdx*rad/h(i)

              gradw(:) = qkd_mag*(0.5d0*qq - 1.0d0)*rij(:)/max(rad,epsilon(rad))
              if(qq.gt.2.0) gradw(:) = 0.0d0
              dr_tmp = dr_tmp + gradw(:)
           end do
           dr(i,:) = dr_tmp
           
           rad = sqrt(dot_product(dr_tmp,dr_tmp))
           if(rad.gt.0.1d0*h0) then
              dr(i,:) = dr_tmp*0.1d0*h0/rad
           end if
           
        end do
        !$OMP END PARALLEL DO
        
        !! Move particles...
        !$OMP PARALLEL DO
        do i=nb+1,npfb
           if(allocated(dontshift))then
              if(dontshift(i).eq.0)then
                 rp(i,:) = rp(i,:) - dr(i,:)
              end if
           else
              rp(i,:) = rp(i,:) - dr(i,:)              
           end if
        end do
        !$OMP END PARALLEL DO
!        write(6,*) maxval(dr(1:npfb-nb,2))/dx

        deallocate(ij_link,ij_count)
        if(allocated(ij_count4)) deallocate(ij_count4)
     end do

     deallocate(dr)
     return
  end subroutine iteratively_shift
!! ------------------------------------------------------------------------------------------------  
  subroutine create_mirror_particles
     !! Subroutine loops over all nodes, and creates mirrors for those
     !! near periodic or symmetric domain limits, for a square domain.
     !! -- NB:
     !! -----:  irelation(j)=i where i is the parent-node of node j
    real(rkind),dimension(2) :: rcorn
    real(rkind) :: cdist
    integer(ikind) :: i,j,imp,k
    integer(ikind) :: nmirror,nmirror_esti
      
    nmirror_esti = npfb*4  ! Estimate for max number of mirrors
    if(allocated(irelation)) deallocate(irelation,vrelation)
    allocate(irelation(npfb+1:npfb+nmirror_esti)) 
    allocate(vrelation(npfb+1:npfb+nmirror_esti))     
    imp = 0     
         

    do i=1,npfb
       
       !! LEFT AND RIGHT BOUNDARIES
       if(rp(i,1).le.xmin+ss*h(i))then ! Close to left bound
          if(xbcond.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k) = 1
             rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2)=rp(i,2)
          else if(xbcond.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k) = 2
             rp(k,1) = 2.0d0*xmin - rp(i,1);rp(k,2)=rp(i,2)
          end if
       end if   
       
       if(rp(i,1).ge.xmax-ss*h(i))then ! Close to right bound
          if(xbcond.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k) = 1
             rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2)=rp(i,2)
          else if(xbcond.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k) = 2
             rp(k,1) = 2.0d0*xmax - rp(i,1);rp(k,2)=rp(i,2)
          end if
       end if 
       
       !! UPPER AND LOWER BOUNDARIES
       if(rp(i,2).le.ymin+ss*h(i))then ! Close to lower bound
          if(ybcond.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k) = 1
             rp(k,1) = rp(i,1);rp(k,2)=rp(i,2) + ymax - ymin
          else if(ybcond.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k) = 3
             rp(k,1) = rp(i,1);rp(k,2)= 2.0d0*ymin - rp(i,2)
          end if
       end if   
       
       if(rp(i,2).ge.ymax-ss*h(i))then ! Close to upper bound
          if(ybcond.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k) = 1
             rp(k,1) = rp(i,1);rp(k,2)=rp(i,2) - ymax + ymin
          else if(ybcond.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k) = 3
             rp(k,1) = rp(i,1);rp(k,2)= 2.0d0*ymax - rp(i,2)
          end if
       end if                
       !! CORNER BOUNDARIES
       rcorn = (/xmin,ymin/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h(i))then  !! Close to lower left corner
          if(xbcond.ne.0.and.ybcond.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond.eq.1.and.ybcond.eq.1)then
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = rp(i,2) + ymax - ymin
                vrelation(k)=1                
             else if(xbcond.eq.2.and.ybcond.eq.1)then
                rp(k,1) = 2.0d0*xmin - rp(i,1);rp(k,2) = rp(i,2) + ymax - ymin
                vrelation(k)=2                
             else if(xbcond.eq.1.and.ybcond.eq.2)then
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = 2.0d0*ymin - rp(i,2)
                vrelation(k)=3                
             else if(xbcond.eq.2.and.ybcond.eq.2)then
                rp(k,1) = 2.0d0*xmin - rp(i,1);rp(k,2) = 2.0d0*ymin - rp(i,2)
                vrelation(k)=4                
             end if
          end if
       end if
       
       rcorn = (/xmax,ymin/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h(i))then  !! close to lower right corner
          if(xbcond.ne.0.and.ybcond.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond.eq.1.and.ybcond.eq.1)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = rp(i,2) + ymax - ymin
                vrelation(k)=1                
             else if(xbcond.eq.2.and.ybcond.eq.1)then
                rp(k,1) = 2.0d0*xmax - rp(i,1);rp(k,2) = rp(i,2) + ymax - ymin
                vrelation(k)=2                
             else if(xbcond.eq.1.and.ybcond.eq.2)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = 2.0d0*ymin - rp(i,2)
                vrelation(k)=3                
             else if(xbcond.eq.2.and.ybcond.eq.2)then
                rp(k,1) = 2.0d0*xmax - rp(i,1);rp(k,2) = 2.0d0*ymin - rp(i,2)
                vrelation(k)=4                
             end if
          end if
       end if
       
       rcorn = (/xmin,ymax/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h(i))then  !! close to upper left corner
          if(xbcond.ne.0.and.ybcond.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond.eq.1.and.ybcond.eq.1)then
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = rp(i,2) - ymax + ymin
                vrelation(k)=1                
             else if(xbcond.eq.2.and.ybcond.eq.1)then
                rp(k,1) = 2.0d0*xmin - rp(i,1);rp(k,2) = rp(i,2) - ymax + ymin
                vrelation(k)=2                
             else if(xbcond.eq.1.and.ybcond.eq.2)then
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = 2.0d0*ymax - rp(i,2)
                vrelation(k)=3                
             else if(xbcond.eq.2.and.ybcond.eq.2)then
                rp(k,1) = 2.0d0*xmin - rp(i,1);rp(k,2) = 2.0d0*ymax - rp(i,2)
                vrelation(k)=4                
             end if
          end if
       end if
       
       rcorn = (/xmax,ymax/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h(i))then  !! Close to upper right corner
          if(xbcond.ne.0.and.ybcond.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond.eq.1.and.ybcond.eq.1)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = rp(i,2) - ymax + ymin
                vrelation(k)=1                
             else if(xbcond.eq.2.and.ybcond.eq.1)then
                rp(k,1) = 2.0d0*xmax - rp(i,1);rp(k,2) = rp(i,2) - ymax + ymin
                vrelation(k)=2                
             else if(xbcond.eq.1.and.ybcond.eq.2)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = 2.0d0*ymax - rp(i,2)
                vrelation(k)=3                
             else if(xbcond.eq.2.and.ybcond.eq.2)then
                rp(k,1) = 2.0d0*xmax - rp(i,1);rp(k,2) = 2.0d0*ymax - rp(i,2)
                vrelation(k)=4                
             end if
          end if
       end if       
    end do       
              
    nmirror = imp
    np = npfb + nmirror      
   
    return
  end subroutine create_mirror_particles
!! ------------------------------------------------------------------------------------------------      
  subroutine reapply_mirror_bcs
     integer(ikind) :: i,j
    
    
     !! Update properties in the boundary particles
     !$OMP PARALLEL DO PRIVATE(i)
     do j=npfb+1,np
        i = irelation(j)
        ro(j) = ro(i)
        if(vrelation(j).eq.1)then
           u(j) = u(i)
           v(j) = v(i) 
        else if(vrelation(j).eq.2)then
           u(j) = -u(i)
           v(j) = v(i)        
        else if(vrelation(j).eq.3)then
           u(j) = u(i)
           v(j) = -v(i) 
        else if(vrelation(j).eq.4)then
           u(j) = -u(i)
           v(j) = -v(i) 
        end if   
        
        Yspec(j)=Yspec(i)
        if(allocated(E)) E(j)=E(i)
     end do
     !$OMP END PARALLEL DO     


     return
  end subroutine reapply_mirror_bcs  
!! ------------------------------------------------------------------------------------------------  
!! ---------------------------------
  function vardx(dx0,dr) result(tmp)   !! This function sets how dx changes with radius...
    real(rkind),intent(in) :: dx0,dr
    real(rkind) :: tmp
    tmp = dx0*hmin/h0 + dx0*dr**0.5
  end function vardx
!! ---------------------------------
!! ------------------------------------------------------------------------------------------------  
  subroutine create_particles_disc
     integer(ikind) :: i,tmp_i,tmp_j,j,imp,nring
     real(rkind) :: nsx,nsy,dr,rmin,rmax,tn,dbound,thta,dx0,ns

     time = 0.0d0
     xmax = 1.6d0;xmin = -1.6d0;ymax=xmax;ymin=xmin
     dx = 2.0*pi/4.0/nx  ! set particle spacing based 4*nx particles round disc
     h0 = hovdx*dx;sup_size = 2.0d0*ss*h0;h2=h0*h0;h3=h2*h0;ss2=sup_size*sup_size  !h, support 
     eta = 1.0d-8*h0;eta2 = eta*eta;eta3=eta*eta2
     hmin = h0/10.0d0
     dx0 = dx
     
     tmp_i = ceiling(2.0*pi*dble(nx)**2)
     nb=0;nb_n=0 
     allocate(rp(tmp_i,dims))
     allocate(h(tmp_i));h=h0   
     allocate(dontshift(tmp_i));dontshift=0 
      

     !! Create Dirichlet particles outer circle
     do i=1,4*nx+1
        thta = dble((i-1)*dx)/1.0d0
        nb=nb+1;rp(nb,1)=cos(thta);rp(nb,2)=sin(thta)
     end do
     npfb=nb
     !! Internal nodes
     call random_seed()
     dr = 1.0d0 - dx0
     nring = 0
     dx = dx0
     do while(dr.gt.0.5d0*dx0*hmin/h0)
        j = floor(dr*2.0*pi/dx)
        nring=nring+1
        dx = vardx(dx0,dr)        
        do i=1,j
           thta = dble(i-1)*2*pi/dble(j)
           npfb=npfb+1;
           if(nring.le.5)then
              rp(npfb,1)=dr*cos(thta)
              rp(npfb,2)=dr*sin(thta) 
              dontshift(npfb)=1    
           else          
              call random_number(ns);ns = (ns - 0.5)*dx*tmp_noise;rp(npfb,1)=dr*cos(thta)+ns
              call random_number(ns);ns = (ns - 0.5)*dx*tmp_noise;rp(npfb,2)=dr*sin(thta)+ns
           end if
           h(npfb)=hovdx*dx
        end do     
        if(nring.eq.3) nb=npfb !! Make first 3 rings also bound
        dr = dr - dx        
     end do
     np=npfb

     !! Ghost nodes for shifting...
     do j=1,5
     do i=1,4*nx+1
        thta = dble((i-1)*dx0)/1.0d0
        dr = 1.0d0 + dble(j)*dx0
        np=np+1;rp(np,1)=dr*cos(thta);rp(np,2)=dr*sin(thta)
     end do  
     end do
     
open(unit=30,file='fort.30')
do i=1,np
write(30,*) rp(i,:)
end do
flush(30);close(30)        

     call iteratively_shift(10)
    
     !! Set u(i) to a given differentiable function - this is set in sphtools
     allocate(u(np)) 
     allocate(v(np));allocate(w(np))
     !$OMP PARALLEL DO 
     do i=1,np  
        u(i) = ftn(rp(i,1),rp(i,2))
     end do
     !$OMP END PARALLEL DO
         

open(unit=31,file='fort.31')
do i=1,np
write(31,*) rp(i,:)
end do
flush(31);close(31)

     deallocate(dontshift)
     return
  end subroutine create_particles_disc  
!! ------------------------------------------------------------------------------------------------
end module nodes
