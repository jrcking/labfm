program labfm
  use kind_parameters
  use common_parameter
  use common_2d
  use omp_lib
  use nodes
  use neighbours
  use moments
  use basic_convergence_studies
  use burgers_equation
  use ns_equations
  use filtering
  implicit none

  integer(ikind) :: k,kk

  call initial_setup

  !! Loop over a range of resolutions
  nx = 20!! 1/2 the initial resolution
  do k=1,10
       !! Create the particles and give initial values
     nx = nx*2  !! Increase the resolution by a factor of 2 each time...

!     call create_particles_banalytic
     call create_particles_bperiodic
!     call create_particles_bperiodic_varh

     !! Build the neighbour lists
     call find_neighbours
     
     !! Adapt the stencils
!     call adapt_stencils
     
     !! Restrict stencil to the first XX neighbours
!     ij_count(:)=70

     !! Calculate all the interparticle weights and any moments we might need
     !! This is the key part of LABFM
     call calc_interparticle_weights
     call filter_coefficients 

     !! Call subroutine to do whatever test we choose...
!     call gradient_convergence_test
!     call laplacian_convergence_test
!     call freq_response_test             !! Set nx=~80 and comment out nx=nx*2 above
     call filter_test
!     call filter_some_noise
!     call vortex_resolve_test
!     call stability_test
!     call solve_burgers_equation
!     call solve_ns_equations

     call output_uv(k)

     !! Deallocate particle properties and neighbour lists
     deallocate(rp,u);if(allocated(v)) deallocate(v);if(allocated(w)) deallocate(w)
     deallocate(h)
     deallocate(ro,Yspec)
     deallocate(ij_count,ij_link)
     if(allocated(irelation)) deallocate(irelation)
     if(allocated(ij_count4)) deallocate(ij_count4)
     if(allocated(ibtype)) deallocate(ibtype)
     if(allocated(ij_w_grad)) deallocate(ij_w_grad,ij_w_lap,ij_w_hyp,ij_w_hyp2)
     if(allocated(filter_coeff)) deallocate(filter_coeff)
     if(allocated(hqw)) deallocate(hqw)

  end do

  close(1)
  stop
end program labfm
!! ------------------------------------------------------------------------------------------------
subroutine initial_setup  
  use kind_parameters
  use common_parameter
  use common_2d

  !! Domain size
  xmin = -0.5d0;xmax = 0.5d0!*2.0d0*pi
  ymin=xmin;ymax=xmax
  !! Time begins at zero
  time = 0.0d0
  
  !! Characteristic length scale...
  lambda = (xmax - xmin)*1.0d0

  !! Particles per smoothing length (depends on order)
#if order==2
  hovdx = 1.2d0  
#elif order==3
  hovdx = 1.4d0
#elif order==4
  hovdx = 1.4d0
#elif order==5
  hovdx = 2.0d0
#elif order==6
  hovdx = 2.2d0
#elif order==7
  hovdx = 2.4d0
#elif order==8
  hovdx = 2.4d0  !2.7d0 for variable resolution
#elif order==9
  hovdx = 2.9d0
#elif order==10
  hovdx = 3.1d0
#elif order==11
  hovdx = 3.4d0
#elif order==12
  hovdx = 3.7d0
#else
  write(6,*) "Compiled with no order specified. Stopping"
  Stop
#endif

  
  !! For asymmetric stencils these three lines are modified
  hovdx_av=hovdx  
  hovdx_max = hovdx
  hovdx_min = hovdx!2.0d0
  
  !! Stencil size 
  ss = 2.0
  nplink = 4.0*ss*ss*hovdx*hovdx  !! nplink allows for square stencil with side length 2*ss

  !! Level of noise in node distribution
  tmp_noise = 0.5d0!0.5d0 
  
  !! A Reynolds number?? (also parameter for some Poisson stuff)
  Re=1.0d1

  !! Powers of pi and wavelength  
  pi2 = pi*pi;pi3=pi*pi2;pi4=pi3*pi
  l2 = lambda*lambda;l3=l2*lambda;l4=l3*lambda

  !! Open a files for outputting
  open(unit=1,file='./data_out/L2norm')
  open(unit=10,file='./data_out/pde_L2norm')
  open(unit=11,file='./data_out/l2_time')
  open(unit=13,file='./data_out/uv/time_out')
  open(unit=14,file='./data_out/DM_eigens')
  open(unit=15,file='./data_out/filtered_stats')
  
end subroutine initial_setup
