module analytic_functions
  use kind_parameters
  use common_2d
  implicit none

contains

!! ftype 1 - exponentially decaying sinsin (vorticity in TG flow)
!! ftype 2 - some 6th order polynomial in (x-a) (y-b) - offset as George suggested.
!! ftype 3 - Burgers travelling wave initial conditions and solution
!! ftype 4 - Polynomial peak on unit domain w/ Homogeneous Dirichlet conds for Poisson
!! ftype 5 - first X terms of square wave fourier series
!! ftype 6 - Re-entrant corner problem for Laplace eqn...
!! ftype 7 - exponential peak problem for Poisson 
!! ftype 8 - Oscillatory for Helmholtz...
#define ftype 2
!! ------------------------------------------------------------------------------------------------
#if ftype==1
  function ftn(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    tmp = sin(2.0*pi*x/lambda)*sin(2.0*pi*y/lambda)*exp(-8.0*kappa*pi2*time/l2)
  end function ftn
  function dfdx(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    tmp = (2.0*pi/lambda)*cos(2.0*pi*x/lambda)*sin(2.0*pi*y/lambda)*exp(-8.0*kappa*pi2*time/l2)
  end function dfdx
  function dfdy(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    tmp = (2.0*pi/lambda)*sin(2.0*pi*x/lambda)*cos(2.0*pi*y/lambda)*exp(-8.0*kappa*pi2*time/l2)
  end function dfdy
  function d2fdx2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    tmp = -(4.0*pi2/(l2))*sin(2.0*pi*x/lambda)*sin(2.0*pi*y/lambda)*exp(-8.0*kappa*pi2*time/l2)
  end function d2fdx2
  function d2fdy2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    tmp = -(4.0*pi2/(l2))*sin(2.0*pi*x/lambda)*sin(2.0*pi*y/lambda)*exp(-8.0*kappa*pi2*time/l2)
  end function d2fdy2
!! ------------------------------------------------------------------------------------------------
#elif ftype==2
  !!! 16th order polynomial in x and y
  function ftn(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5
    x1=x-0.1453*lambda;y1=y-0.16401*lambda
    x2=x1*x1;x3=x2*x1;x4=x3*x1;x5=x4*x1;y2=y1*y1;y3=y2*y1;y4=y3*y1;y5=y4*y1;
    tmp = 1.0 + x1 + x2 + x3 + x4 + x5 + x1*x5 + y1 + y2 + y3 + y4 + y5 + y1*y5 + y4*x4 + y4*y4*x4*x4
  end function ftn
  function dfdx(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5
    x1=x-0.1453*lambda;y1=y-0.16401*lambda
    x2=x1*x1;x3=x2*x1;x4=x3*x1;x5=x4*x1;y2=y1*y1;y3=y2*y1;y4=y3*y1;y5=y4*y1;
    tmp = 1.0 + 2.0*x1 + 3.0*x2 + 4.0*x3 + 5.0*x4 + 6.0*x5 + 4.0*y4*x3 + 8.0*y4*y4*x4*x3

  end function dfdx
  function dfdy(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5
    x1=x-0.1453*lambda;y1=y-0.16401*lambda
    x2=x1*x1;x3=x2*x1;x4=x3*x1;x5=x4*x1;y2=y1*y1;y3=y2*y1;y4=y3*y1;y5=y4*y1;
    tmp = 1.0 + 2.0*y1 + 3.0*y2 + 4.0*y3 + 5.0*y4 + 6.0*y5 + 4.0*y3*x4 + 8.0*y4*y3*x4*x4

  end function dfdy
  function d2fdx2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5
    x1=x-0.1453*lambda;y1=y-0.16401*lambda
    x2=x1*x1;x3=x2*x1;x4=x3*x1;x5=x4*x1;y2=y1*y1;y3=y2*y1;y4=y3*y1;y5=y4*y1;
    tmp = 2.0 + 6.0*x1 + 12.0*x2 + 20.0*x3 + 30.0*x4 + 12.0*y4*x2 + 56.0*y4*y4*x3*x3

  end function d2fdx2
  function d2fdy2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5
    x1=x-0.1453*lambda;y1=y-0.16401*lambda
    x2=x1*x1;x3=x2*x1;x4=x3*x1;x5=x4*x1;y2=y1*y1;y3=y2*y1;y4=y3*y1;y5=y4*y1;
    tmp = 2.0 + 6.0*y1 + 12.0*y2 + 20.0*y3 + 30.0*y4 + 12.0*y2*x4 + 56.0*y3*y3*x4*x4

  end function d2fdy2
!! ------------------------------------------------------------------------------------------------
#elif ftype==3
  !! Solution to Burgers equation
  function ftn(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,alph
    alph = Re*(-time-4.0*x+4.0*y)/32.0
    tmp = 0.25/(1.0 + exp(alph))
  end function ftn
  function dfdx(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,alph
    alph = Re*(-time-4.0*x+4.0*y)/32.0
    tmp = (1.0/32.0)*Re*exp(alph)/(1.0 + exp(alph))**2.0
  end function dfdx
  function dfdy(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,alph
    alph = Re*(-time-4.0*x+4.0*y)/32.0
    tmp = -(1.0/32.0)*Re*exp(alph)/(1.0 + exp(alph))**2.0
  end function dfdy
  function d2fdx2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,alph
    alph = Re*(-time-4.0*x+4.0*y)/32.0
    tmp = (Re*Re/256.0)*(exp(2.0*alph)-exp(alph))/(1.0 + exp(alph))**3.0
    if(i_uv_burgers.eq.1) tmp = -tmp
  end function d2fdx2
  function d2fdy2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,alph
    alph = Re*(-time-4.0*x+4.0*y)/32.0
    tmp = (Re*Re/256.0)*(exp(2.0*alph)-exp(alph))/(1.0 + exp(alph))**3.0
    if(i_uv_burgers.eq.1) tmp = -tmp
  end function d2fdy2
!! ------------------------------------------------------------------------------------------------
#elif ftype==4
  !! Polynomial peak for Poisson problem  
  function ftn(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    tmp = (2.0d0**(4.0d0*Re))*(x**Re)*((1.0d0-x)**Re)*(y**Re)*((1.0d0-y)**Re)
  end function ftn
  function dfdx(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,tmp1,tmp2
    tmp1 = (2.0d0**(4.0d0*Re))*(y**Re)*((1.0d0-y)**Re)
    tmp2 = (1.0d0-2.0d0*x)*Re*(x-x*x)**(Re-1.0d0)
    tmp = tmp1*tmp2
  end function dfdx
  function dfdy(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,tmp1,tmp2
    tmp1 = (2.0d0**(4.0d0*Re))*(x**Re)*((1.0d0-x)**Re)
    tmp2 = (1.0d0-2.0d0*y)*Re*(y-y*y)**(Re-1.0d0)
    tmp = tmp1*tmp2  
  end function dfdy
  function d2fdx2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,tmp1,tmp2
    tmp1 = (2.0d0**(4.0d0*Re))*(y**Re)*((1.0d0-y)**Re)
    tmp2 = (Re-1.0d0)*Re*((1.0d0-2.0d0*x)**2.0)*(x-x*x)**(Re-2.0d0) - 2.0*Re*(x-x*x)**(Re-1.0d0)    
    tmp = tmp2*tmp1
  end function d2fdx2
  function d2fdy2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,tmp2,tmp1
    tmp1 = (2.0d0**(4.0d0*Re))*(x**Re)*((1.0d0-x)**Re)
    tmp2 = (Re-1.0d0)*Re*((1.0d0-2.0d0*y)**2.0)*(y-y*y)**(Re-2.0d0) - 2.0*Re*(y-y*y)**(Re-1.0d0)    
    tmp = tmp2*tmp1
  end function d2fdy2
!! ------------------------------------------------------------------------------------------------
#elif ftype==5
  function ftn(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    real(rkind) :: xarg
    xarg = 2.0*pi*(x/lambda - 0.25d0)
    tmp = sin(2.0*pi*y/lambda)*(4.0/pi)*(sin(xarg)+sin(3.0*xarg)/3.0 + sin(5.0*xarg)/5.0 + sin(7.0*xarg)/7.0 + sin(9.0*xarg)/9.0 &
          + sin(11.0*xarg)/11.0 + sin(13.0*xarg)/13.0 + sin(15.0*xarg)/15.0)! + sin(17.0*xarg)/17.0 + sin(19.0*xarg)/19.0)
  end function ftn
  function dfdx(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    real(rkind) :: xarg
    xarg = 2.0*pi*(x/lambda - 0.25d0)
    tmp = (2.0*pi/lambda)*sin(2.0*pi*y/lambda)*(4.0/pi)* &
          (cos(xarg)+cos(3.0*xarg) + cos(5.0*xarg) + cos(7.0*xarg) + cos(9.0*xarg) + cos(11.0*xarg) + cos(13.0*xarg) + &
           cos(15.0*xarg))! + cos(17.0*xarg) + cos(19.0*xarg) )
  end function dfdx
  function dfdy(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    real(rkind) :: xarg
    xarg = 2.0*pi*(x/lambda - 0.25d0)
    tmp = (2.0*pi/lambda)*cos(2.0*pi*y/lambda)* &
          (4.0/pi)*(sin(xarg)+sin(3.0*xarg)/3.0 + sin(5.0*xarg)/5.0 + sin(7.0*xarg)/7.0 + sin(9.0*xarg)/9.0 &
          + sin(11.0*xarg)/11.0 + sin(13.0*xarg)/13.0 + sin(15.0*xarg)/15.0)! + sin(17.0*xarg)/17.0 + sin(19.0*xarg)/19.0)
  end function dfdy
  function d2fdx2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    real(rkind) :: xarg
    xarg = 2.0*pi*(x/lambda - 0.25d0)
    tmp = -(4.0*pi2/l2)*sin(2.0*pi*y/lambda)*(4.0/pi)* &
          (sin(xarg) + 3.0*sin(3.0*xarg) + 5.0*sin(5.0*xarg) + 7.0*sin(7.0*xarg) + 9.0*sin(9.0*xarg) + 11.0*sin(11.0*xarg) + &
           13.0*sin(13.0*xarg) + 15.0*sin(15.0*xarg))! + 17.0*sin(17.0*xarg) + 19.0*sin(19.0*xarg))
  end function d2fdx2
  function d2fdy2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    real(rkind) :: xarg
    xarg = 2.0*pi*(x/lambda - 0.25d0)
    tmp = -(4.0*pi2/l2)*sin(2.0*pi*y/lambda)* &
          (4.0/pi)*(sin(xarg)+sin(3.0*xarg)/3.0 + sin(5.0*xarg)/5.0 + sin(7.0*xarg)/7.0 + sin(9.0*xarg)/9.0 &
          + sin(11.0*xarg)/11.0 + sin(13.0*xarg)/13.0 + sin(15.0*xarg)/15.0)! + sin(17.0*xarg)/17.0 + sin(19.0*xarg)/19.0)
  end function d2fdy2
!! ------------------------------------------------------------------------------------------------
#elif ftype==6
  !! Re-entrant corner Poisson problem
  function ftn(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,th
    r = sqrt(x*x+y*y)
    if(x.ge.0.0d0)then
       if(y.ge.0.0d0)then
          th = asin(y/r)
       else
          th = 2.0d0*pi-asin(-y/r)
       endif
    else if(x.lt.0.0d0)then
       th = pi - asin(y/r)
    end if    
    if(x.eq.0.0d0) then
       if(y.le.0.0d0) th=3.0d0*pi/2.0d0
       if(y.gt.0.0d0) th=0.5d0*pi
    end if
    tmp = sin(Re*th)*r**(Re)
  end function ftn
  function dfdx(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,th
    tmp = 0.0
  end function dfdx
  function dfdy(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,th
    tmp = 0.0
  end function dfdy
  function d2fdx2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,th
    tmp = 0.0d0
  end function d2fdx2
  function d2fdy2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,th
    tmp = 0.0d0
  end function d2fdy2 
!! ------------------------------------------------------------------------------------------------
#elif ftype==7
  !! Exponential peak problem
  function ftn(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,rr,xc,yc
    xc = 0.0d0;yc=0.0d0;rr=(x-xc)**2 + (y-yc)**2
    tmp = exp(-Re*rr)
  end function ftn
  function dfdx(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,th
    tmp = 0.0
  end function dfdx
  function dfdy(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,th
    tmp = 0.0
  end function dfdy
  function d2fdx2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,th
    tmp = 0.0d0
  end function d2fdx2
  function d2fdy2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,rr,xc,yc
    xc = 0.0d0;yc=0.0d0;rr=(x-xc)**2 + (y-yc)**2
    tmp = exp(-Re*rr)*(4.0d0*Re*Re*rr - 4.0d0*Re)
  end function d2fdy2 
!! ------------------------------------------------------------------------------------------------
#elif ftype==8
  !! 2D Gaussian
  function ftn(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,a,b
    r=sqrt(x*x+y*y);a=10.0d0;b=2.0d0
    tmp = exp(-a*r**b)
!    tmp = 0.5*(erf(-a*x)+1.0d0)
  end function ftn
  function dfdx(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,a,b
    r=sqrt(x*x+y*y);a=10.0d0;b=2.0d0
    tmp = -a*b*x*((x*x+y*y)**(0.5d0*b-1.0d0))*exp(-a*r**b)
!    tmp = -a*exp(-a*a*x*x)/sqrt(pi)
  end function dfdx
  function dfdy(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,a,b
    r=sqrt(x*x+y*y);a=10.0d0;b=2.0d0
    tmp = -a*b*y*((x*x+y*y)**(0.5d0*b-1.0d0))*exp(-a*r**b)
  end function dfdy
  function d2fdx2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,a,b
    r=sqrt(x*x+y*y);a=10.0d0;b=2.0d0
    tmp = a*b*(a*b*x*x*(r**(2.0d0*b))+((1.0d0-b)*x*x-y*y)*(r**b))*exp(-a*r**b)/(r**4.0d0)
  end function d2fdx2
  function d2fdy2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,r,a,b
    r=sqrt(x*x+y*y);a=10.0d0;b=2.0d0
    tmp = a*b*(a*b*y*y*(r**(2.0d0*b))+((1.0d0-b)*y*y-x*x)*(r**b))*exp(-a*r**b)/(r**4.0d0)  
  end function d2fdy2 
#elif ftype==9
!!---------------------------------------------------------------------------------------------------------
  !!! 7th order polynomial in x only
  function ftn(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,x1,x2,x3,x4,x5
    x1=x-0.1453*lambda
    x2=x1*x1;x3=x2*x1;x4=x3*x1;x5=x4*x1
    tmp = 1.0 + x1 + x2 + x3 + x4 + x5 + x1*x5 + x5*x2
  end function ftn
  function dfdx(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,x1,x2,x3,x4,x5
    x1=x-0.1453*lambda
    x2=x1*x1;x3=x2*x1;x4=x3*x1;x5=x4*x1
    tmp = 1.0 + 2.0*x1 + 3.0*x2 + 4.0*x3 + 5.0*x4 + 6.0*x5 + 7.0*x5*x1
  end function dfdx
  function dfdy(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    tmp = 0.0d0
  end function dfdy
  function d2fdx2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp,x1,x2,x3,x4,x5
    x1=x-0.1453*lambda
    x2=x1*x1;x3=x2*x1;x4=x3*x1;x5=x4*x1
    tmp = 2.0 + 6.0*x1 + 12.0*x2 + 20.0*x3 + 30.0*x4 + 42.0*x5
  end function d2fdx2
  function d2fdy2(x,y) result(tmp)
    real(rkind),intent(in) :: x,y
    real(rkind) :: tmp
    tmp = 0.0d0
  end function d2fdy2  
#endif
!! ------------------------------------------------------------------------------------------------
end module analytic_functions
