module common_2d

  use iso_c_binding
  use common_parameter
  use kind_parameters      

  implicit none

  !! Particle properties
  real(rkind), dimension(:), allocatable, target :: u,v,w,h,filter_coeff,ro,Yspec,E
  real(rkind), dimension(:,:), allocatable, target :: rp,rshift
  real(rkind),dimension(:,:),allocatable :: grad
  real(rkind),dimension(:),allocatable :: lap,div
  real(rkind),dimension(:),allocatable :: dudt
  real(rkind), bind(c) :: dx,dv

  !! flag to indicate the correction type, pde type...
  integer(ikind) :: i_correction_type, i_eqn
  integer(ikind) :: i_uv_burgers !! tell analytic functions which sign needed for diff components of u,v in burgers

  !! TIME
  real(rkind) :: time,dt,time_end
  integer(ikind) :: itime

  !! Neighbours
  integer(ikind),dimension(:),allocatable :: ij_count,ij_count4
  integer(ikind),dimension(:,:),allocatable :: ij_link

  !! Weightings - more memory, computationally cheaper...
  real(rkind),dimension(:,:,:),allocatable :: ij_w_grad
  real(rkind),dimension(:,:),allocatable :: ij_w_lap,ij_w_hyp,ij_w_hyp2
  real(rkind),dimension(:),allocatable :: hqw     !! Hopkins quadrature weights

  !! Parents and boundaries...
  integer(ikind),dimension(:),allocatable :: irelation,ibtype,vrelation
  integer(ikind) :: xbcond,ybcond

  integer :: n_threads   !! Parallelisation!!

  integer(ikind), bind(c) :: np,npfb,npfb_esti,nx,nb,nb_n,nplink  !! Particle numbers

  real(rkind),dimension(3) :: L2norm,L2norm_y   !! L2 norms of error for outputting

  real(rkind), bind(c) :: xmin,xmax,ymin,ymax,tmp_noise

  real(rkind), bind(c) :: h0,sup_size,h3,hovdx,ss,h2,ss2,hovdx_av,hmin
  real(rkind) :: pi2,pi3,pi4,hovdx_max,hovdx_min

  !! Parameters - wavelengths, diffusivities, non-dim numbers etc...
  real(rkind) :: kl,lambda,l2,l3,l4,kappa,Re,BF

  !! Big sparse matrix
  real(rkind), dimension(:), allocatable :: lhs_mat
  integer(ikind), dimension(:), allocatable :: ija

  real(rkind), bind(c) :: eta,eta2,eta3
        
end module common_2d
