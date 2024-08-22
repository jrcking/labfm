module neighbours
  use kind_parameters
  use common_parameter
  use common_2d
  use sphtools
  use omp_lib
  implicit none

  private
  public :: find_neighbours

  integer(ikind), dimension(:), allocatable ::ist, ip, nc
  integer(ikind), dimension(:), allocatable :: cellpart,ic_count
  integer(ikind) :: ncx,ncy,nct
  real(rkind) :: xmaxt,xmint,ymaxt,ymint,unoss
  integer(ikind),dimension(:,:),allocatable :: ic_link
  
 
contains
!! ------------------------------------------------------------------------------------------------
  subroutine find_neighbours
    !! This subroutine does some allocation, then calls routines to put nodes in
    !! cells, reorder them, and look through cells to build neighbour lists..
    integer nct_temp,n_count
    real(rkind) :: q,qlim
    integer(ikind) :: i,j,k
     
    !! Build cells and create cell linked-lists     
    if(.not.allocated(ic_count))then
       call cell_generator
       call neighbour_cells
    end if

    ! allocate linking lists
    allocate(ij_count(npfb))
    allocate(ij_link(npfb,nplink))
    n_count=np
    
    ! allocate cell lists
    nct_temp=nct+ncx+ncy+3
    allocate(ist(nct_temp+1))
    allocate(nc(nct_temp))
    allocate(ip(n_count))
    nc=0

    ! put nodes in different cells
    call divide(1,np) 
    call neighboring_list_parallel
    
    !! Order neighbours by distance
    call order_neighbours
    
    !! Make ij_count4 (smaller stencils)
!    allocate(ij_count4(npfb))
!    do i=1,npfb
!       qlim = (1.7d0*h(i))**two
!       q = zero
!       k=0
!       do while(q.le.qlim)
!          k=k+1
!          j=ij_link(i,k)
!          q=(rp(i,1)-rp(j,1))**two + (rp(i,2)-rp(j,2))**two
!       end do
!       ij_count4(i) = k
!!       write(6,*) ij_count(i),k
!    end do


    deallocate(ist)
    deallocate(nc)
    deallocate(ip)
    deallocate(cellpart)

  end subroutine find_neighbours
!! ------------------------------------------------------------------------------------------------
  subroutine cell_generator
 
     !! cell size and max/min domain size
     !! NB. we use large cells (2.0d0*sup_size) to account for nodes shifted from stencil centre
     unoss = 1.0d0/(sup_size)
     xmaxt = xmax + sup_size;xmint = xmin - sup_size
     ymaxt = ymax + sup_size;ymint = ymin - sup_size
     
     !! Number of cells in each direction
     ncx = int((xmaxt-xmint)*unoss) + 1     
     ncy = int((ymaxt-ymint)*unoss) + 1
     
     nct=ncx*ncy ! Number of cells in a XY sheet

  end subroutine cell_generator
!! ------------------------------------------------------------------------------------------------
  subroutine neighbour_cells
    use common_2d
    integer(ikind) :: ic,icx,icy,jc,nsheet
   
    nsheet = ncx*ncy
   
    !! Allocate the cell-link-lists
    allocate(ic_count(nct),ic_link(nct,9))
 
    !! Loop over all cells
    do ic = 1,nct
       icx = mod(ic,ncx)      !! Position in row of cells
       icy = (ic/ncx) +1   !! Which row...

       ! This cell
       ic_count(ic)=1;ic_link(ic,1)=ic
       
       if(icx+1.le.ncx)then
          !! East cell
          jc = ic+1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          if(icy+1.le.ncy)then
             !! NE cell
             jc = ic + ncx + 1    
             ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if          
          if(icy-1.ge.1)then
             !! SE cell
             jc = ic - ncx + 1
             ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if
       end if       
       if(icx-1.ge.1)then
          !! West cell
          jc = ic - 1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          if(icy+1.le.ncy)then
             !! NW cell
             jc = ic + ncx - 1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if          
          if(icy-1.ge.1)then
             !! SW cell
             jc = ic - ncx - 1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if         
       end if       
       if(icy+1.le.ncy) then
          !! North cell
          jc = ic + ncx
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
       end if      
       if(icy-1.ge.1)then
          !! South cell
          jc = ic - ncx
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
       end if          
    end do
   
    return
  end subroutine neighbour_cells 
!! ------------------------------------------------------------------------------------------------
  subroutine divide(nini,nend) 
    !! This subroutine divides the nodes amongst the cells, and creates
    !! a list of particles ordered by their cell
    integer(ikind), intent(in) :: nini,nend 
    real(rkind) :: ddx,ddy
    integer(ikind) :: nc_ii_max, k, icell, jcell, ii

    nc=0_ikind  ! nc(i,mkind)=0

    nc_ii_max = 0
    allocate(cellpart(nend-nini+1))

    !$omp parallel do private(ddx,ddy,icell,jcell,ii) reduction(+:nc) reduction(max:nc_ii_max)
    do k=nini, nend    ! loop over specified range of particles
       ddx = rp(k,1) - xmint    ! find x,y of particle if origin at BL corner of domain
       ddy = rp(k,2) - ymint
  
       icell = int( ddx *unoss) + 1       !! Cell label in x
       jcell = int( ddy *unoss) + 1       !! Cell label in y
       ii = icell + (jcell - 1)*ncx        !! ii is linear cell position in the matrix of cells

       nc(ii) = nc(ii)+1     ! add 1 to the number of particles in this cell
       !! If the number of particles in a cell is greater than
       !! some upper limit then stop
       nc_ii_max = nc(ii)  ! ensure nc_ii_max is max(nc(ii) for all ii
       cellpart(k)=ii             !! the index of the cell in which particle i is
    end do
    !$omp end parallel do

    !! loop over all cells
    ist(1)=1_ikind
    do ii=1, nct
       ist(ii+1)=ist(ii)+nc(ii)   ! ist(ii) is the starting index for each cell
       nc(ii)=0_ikind             ! erase the nc in each cell here 
    end do

    !! index look-up: ip(j)=k where j is the cell-ordered index and k is the original index
    do k=nini,nend
       ip(ist(cellpart(k))+nc(cellpart(k)))=k   
       nc(cellpart(k))=nc(cellpart(k))+1
    end do 

  end subroutine divide
!! ------------------------------------------------------------------------------------------------
  subroutine  neighboring_list_parallel
    !! This subroutine loops through every node i, finds the cell ic containing
    !! node i, then goes through every cell jc adjacent to ic to check for neighbours  
    integer(ikind) :: i,ic
    integer(ikind) :: jc,kc
    real(rkind),dimension(dims) :: rpi

    ij_count=0
    ij_link=0 !neighbour list arrays
        
    !! Loop over all particles (in parallel)
    !$omp parallel do private(ic,jc)
    do i=1,npfb 
       !! Which cell are we in?
       ic = cellpart(i)         

       rpi(:) = rp(i,:)  !! search locus is shifted to node il
       
       !! Loop over all cells neighbouring cell ic
       do kc=1,ic_count(ic)
          jc = ic_link(ic,kc)
          
          !! Add any neighbours in cell jc to neighbour list of particle i
          call neighbour_node_cell(i,rpi,jc)
   
       end do
    end do
    !$omp end parallel do
  end subroutine neighboring_list_parallel
!! ------------------------------------------------------------------------------------------------
  subroutine neighbour_node_cell(ii,ri,jc)
    !! looks through nodes in cell jc and checks whether they are neighbours of 
    !! node ii which has position (vector) ri.
    integer(ikind),intent(in) :: ii,jc !! jc is cell
    real(rkind),dimension(:),intent(in) :: ri
    real(rkind) :: rr2tmp,rad,stencilsize2
    integer(ikind) :: j,jj,is,ie
    real(rkind),dimension(dims) :: rij
    
    if(nc(jc).ne.0)then !! if the cell isn't empty
       stencilsize2=(h(ii)*ss)**2
       is =ist(jc);ie=ist(jc+1)-1  !! Start and end indices of particles in cell jc
       do jj=is,ie       !! Loop over all particles in cell jc
          j=ip(jj)       !! j is regular index, jj is cell-ordered index

          rij(:) = ri(:)-rp(j,:)
          rr2tmp=dot_product(rij,rij)  !! Distance squared
          if(rr2tmp .le. stencilsize2)then              
             ij_count(ii)=ij_count(ii)+1;              !! Increment count
             ij_link(ii,ij_count(ii))=j                !! add to list  
          endif          
       end do
    end if
  
    return
  end subroutine neighbour_node_cell
!! ------------------------------------------------------------------------------------------------  
  subroutine order_neighbours
    !! This routine sorts the neighbour list in order of increasing distance.
    integer(ikind) :: i,j,k
    real(rkind),dimension(dims) :: rij
    real(rkind),dimension(:),allocatable :: rij2
  
    !! Allocate space for rij2
    allocate(rij2(nplink))
  
    !! Loop over nodes
    !$omp parallel do private(j,k,rij,rij2)
    do i=1,npfb
    
       !! First evaluate square of distance
       rij2 = 0.0d0
       do k=1,ij_count(i)
          j=ij_link(i,k)
          rij = rp(i,:) - rp(j,:)
          rij2(k) = dot_product(rij,rij)
!          rij2(k) = dble(j)
       end do
       
       !! Now sort by value of rij2
       call quicksort_neighbours(i,rij2,1,ij_count(i))       
    
    
    end do
    !$omp end parallel do
    
    deallocate(rij2)
  

    return
  end subroutine order_neighbours
!! ------------------------------------------------------------------------------------------------
  recursive subroutine quicksort_neighbours(ii,a, first, last)
    integer(ikind),intent(in) :: ii  !! Index of node whose neighbours we're sorting
    real(rkind) ::  a(*), x, t
    integer(ikind),intent(in) :: first, last  !! Start and end indices of sort
    integer(ikind) i, j,it

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
       do while (a(i) < x)
          i=i+1
       end do
       do while (x < a(j))
          j=j-1
       end do
       if (i >= j) exit
       t = a(i);  a(i) = a(j);  a(j) = t
       it = ij_link(ii,i); ij_link(ii,i) = ij_link(ii,j);ij_link(ii,j) = it
       i=i+1
       j=j-1
    end do
    if (first.lt.i-1) call quicksort_neighbours(ii,a, first, i-1)
    if (j+1.lt.last)  call quicksort_neighbours(ii,a, j+1, last)
    
  end subroutine quicksort_neighbours
end module neighbours
