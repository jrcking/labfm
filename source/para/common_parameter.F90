module common_parameter
  use kind_parameters
  implicit none 

  integer(ikind) ,parameter :: npar=140000
  !! NOTE: npmax should be npar*nplink
  integer(ikind) ,parameter :: dims = 2

  !! integer,parameter::nmax=500000

  real(rkind), parameter :: pi=3.141592653589793238462643383279502884197d0
  real(rkind), parameter :: sqrt2=dsqrt(2.0d0)
  real(rkind), parameter :: oosqrt2=1.0d0/dsqrt(2.0d0)
!!  real(rkind), parameter :: pi=4.0_rkind * DATAN(1.0_rkind)
  real(rkind), parameter :: zero=0.0d0 
  real(rkind), parameter :: one=1.0d0
  real(rkind), parameter :: two=2.0d0
  real(rkind), parameter :: three=3.0d0
  real(rkind), parameter :: four=4.0d0 
  real(rkind), parameter :: half = one/two     
  
  !! Runge Kutta coefficients ---------------------------------------------------------------------
  !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
  real(rkind),parameter :: rk3_4s_2r_a21 = 11847461282814.0d0/36547543011857.0d0
  real(rkind),parameter :: rk3_4s_2r_a32 = 3943225443063.0d0/7078155732230.0d0
  real(rkind),parameter :: rk3_4s_2r_a43 = -346793006927.0d0/4029903576067.0d0
  real(rkind),parameter :: rk3_4s_2r_b1 = 1017324711453.0d0/9774461848756.0d0
  real(rkind),parameter :: rk3_4s_2r_b2 = 8237718856693.0d0/13685301971492.0d0
  real(rkind),parameter :: rk3_4s_2r_b3 = 57731312506979.0d0/19404895981398.0d0
  real(rkind),parameter :: rk3_4s_2r_b4 = -101169746363290.0d0/37734290219643.0d0
  real(rkind),parameter :: rk3_4s_2r_bh1 = 15763415370699.0d0/46270243929542.0d0
  real(rkind),parameter :: rk3_4s_2r_bh2 = 514528521746.0d0/5659431552419.0d0
  real(rkind),parameter :: rk3_4s_2r_bh3 = 27030193851939.0d0/9429696342944.0d0
  real(rkind),parameter :: rk3_4s_2r_bh4 = -69544964788955.0d0/30262026368149.0d0
  real(rkind),dimension(3),parameter :: rk3_4s_2r_a=(/rk3_4s_2r_a21,rk3_4s_2r_a32,rk3_4s_2r_a43/)
  real(rkind),dimension(4),parameter :: rk3_4s_2r_b=(/rk3_4s_2r_b1,rk3_4s_2r_b2,rk3_4s_2r_b3,rk3_4s_2r_b4/)
  real(rkind),dimension(4),parameter :: rk3_4s_2r_bh=(/rk3_4s_2r_bh1,rk3_4s_2r_bh2,rk3_4s_2r_bh3,rk3_4s_2r_bh4/)
  real(rkind),dimension(4),parameter :: rk3_4s_2r_bmbh = rk3_4s_2r_b - rk3_4s_2r_bh
  real(rkind),parameter :: rk3_4s_2r_c1=zero
  real(rkind),parameter :: rk3_4s_2r_c2=rk3_4s_2r_a21
  real(rkind),parameter :: rk3_4s_2r_c3=rk3_4s_2r_a32 + rk3_4s_2r_b1
  real(rkind),parameter :: rk3_4s_2r_c4=one
  real(rkind),dimension(4),parameter :: rk3_4s_2r_c=(/rk3_4s_2r_c1,rk3_4s_2r_c2,rk3_4s_2r_c3,rk3_4s_2r_c4/)   
  

end module common_parameter
