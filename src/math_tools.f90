! ------------------------------------------------------- !
!                                                         !    
! This is part of the SPA (soil-plant-atmosphere) model.  !
!                                                         !    
!  Copyright (C) 2010--2012 SPA contributors.             !
!  See docs/Copyright_Detail for copying conditions.      !
!                                                         !
! Model version:                                          !
! ------------------------------------------------------- !


module math_tools

  !! This module contains general maths routines necessary for !!
  !! calculating derivatives and bisections.  See individual   !!
  !! routines for details.                                     !!

  use gv_scale_declarations, only: nmax

  implicit none

  ! Set the default state of variable/procedures to be private.  This
  ! limits their use to within this module, and avoids accidental overlap
  ! by use of the same name in different modules. 
  private
  ! Once set to private, we then have to explicitly declare anything we
  ! want to 'share' with the rest of SPA.
  ! Procedures..
  public :: ode_int, zbrent
  ! Variables..
  public :: dxsav, kmax

  ! variables shared across SPA..
  integer :: kmax   !
  real    :: dxsav  !

  ! variables private to this module..
  integer,parameter :: kmaxx  = 200,  & ! descriptions
                       maxstp = 10000   !   would
  real,parameter    :: tiny   = 1.e-30  !        nice!

contains
  !
  !----------------------------------------------------------------------
  !
  ! public procedures first.
  !
  !----------------------------------------------------------------------
  !
  subroutine ode_int( called_from , ystart , nvar , x1 , x2 , eps , h1 , hmin , nok , nbad , derivs )

    ! This is an integrator for ordinary differential equations. We use !
    ! the RUNGE_KUTTA to track the dynamic behaviour of various model   !
    ! state variables.. uch a leaf water potential, and water stored on !
    ! leaf surfaces. The RUNGE_KUTTA finds the time step that ensures   !
    ! dynamics are smooth and that the relevant feedbacks are properly  !
    ! incorporated. For a full description see Press et al. (1986).     !

    use gv_scale_declarations, only: nmax
    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
    integer,intent(in)          :: nvar 
    real,intent(in)             :: h1, hmin, x1, x2
    real,intent(inout)          :: ystart(nvar), eps
    integer,intent(out)         :: nbad, nok

    ! Interfaces are the correct way to pass procedures as arguments.
    ! (note 'derivs' is actually one of canopy_water_store,soil_water_store,lwp_diff_eqn)
!    external :: derivs
    interface
       subroutine derivs(time,y,dydt)
         use gv_scale_declarations, only: nmax
         real,intent(in)    :: time
         real,intent(in)    :: y(nmax)
         real,intent(out)   :: dydt(nmax)
       end subroutine derivs
    end interface

    ! local variables..
    integer :: i, kount, nstp
    real    :: h, hdid, hnext, x, xsav ,         &
               dydx(nmax), y(nmax), yscal(nmax), &
               xp(kmaxx), yp(nmax,kmaxx) 

    ! calculations..
    x = x1
    h = sign( h1 , x2-x1 )
    nok = 0
    nbad = 0
    kount = 0
    do i = 1 , nvar
       y(i) = ystart(i)
    enddo
    if ( kmax .gt. 0 ) xsav = x - 2. * dxsav
    do nstp = 1 , MAXSTP
       call derivs( x , y , dydx )
       do i = 1 , nvar
          yscal(i) = abs(y(i))+abs(h*dydx(i))+TINY
       enddo
       if ( kmax .gt. 0 ) then
          if ( abs( x - xsav ) .gt. abs( dxsav ) ) then
             if ( kount .lt. kmax - 1 ) then
                kount = kount + 1
                xp(kount) = x
                do i = 1 , nvar
                   yp(i,kount) = y(i)
                enddo
                xsav = x
             endif
          endif
       endif
       if ( (x+h-x2) * (x+h-x1) .gt. 0. ) h = x2 - x

       call runge_kutta_q_step( trim(called_from)//":ode_int" , y , dydx , nvar , x , h , eps , yscal , hdid , hnext , derivs )
       if (hdid.eq.h) then
          nok = nok+1
       else
          nbad = nbad+1
       endif
       if ( (x-x2) * (x2-x1) .ge. 0. ) then
          do i = 1 , nvar
             ystart(i) = y(i)
          enddo
          if ( kmax .ne. 0 ) then
             kount = kount + 1
             xp(kount) = x
             do i = 1 , nvar
                yp(i,kount) = y(i)
             enddo
          endif
          return
       endif
       if ( abs(hnext) .lt. hmin ) then
         write(message,*) "stepsize smaller than permitted minimum in ode_int",new_line('x'),&
                          "ode_int was called by: ",trim(called_from)
         call write_log( trim(message) , msg_warning , __FILE__ , __LINE__ )
       endif
       h  =  hnext
    enddo

  end subroutine ode_int
  !
  !----------------------------------------------------------------------
  !
  real function zbrent( called_from , func , x1 , x2 , tol )

    ! This is a bisection routine. When ZBRENT is called, we provide a    !
    !  reference to a particular function and also two values which bound !
    !  the arguments for the function of interest. ZBRENT finds a root of !
    !  the function (i.e. the point where the function equals zero), that !
    !  lies between the two bounds.                                       !
    ! For a full description see Press et al. (1986).                     !

    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
    real,intent(in)             :: tol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
       real function func( xval )
         real,intent(in) :: xval
       end function func
    end interface

!   external :: derivs
!    interface
!       subroutine derivs(time,y,dydt)
!         use gv_scale_declarations, only: nmax
!         real,intent(in)    :: time
!         real,intent(in)    :: y(nmax)
!         real,intent(out)   :: dydt(nmax)
!       end subroutine derivs
!    end interface


    ! local variables..
    integer            :: iter
    integer,parameter  :: ITMAX = 30
    real               :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    real,parameter     :: EPS = 3.e-8

    ! calculations...
    a  = x1
    b  = x2
    fa = func( a )
    fb = func( b )

    ! Check that we haven't (by fluke) already started with the root..
    if ( fa .eq. 0. ) then
      zbrent = a
      return
    elseif ( fb .eq. 0. ) then
      zbrent = b
      return
    endif

    ! Ensure the supplied x-values give y-values that lie either
    ! side of the root and if not flag an error message...
    if ( sign(1.0,fa) .eq. sign(1.0,fb) ) then
       fa = func( a )
       fb = func( b )
       write(message,*)"Supplied values must bracket the root of the function.",new_line('x'),  &
          "     ","You supplied x1:",x1,new_line('x'),                     &
          "     "," and x2:",x2,new_line('x'),                             &
          "     "," which give function values of fa :",fa,new_line('x'),  &
          "     "," and fb:",fb," .",new_line('x'),                        &
          " zbrent was called by: ",trim(called_from)
       call write_log( trim(message) , msg_error , __FILE__ , __LINE__ )
       fa = func( a )
       fb = func( b )
    endif
    c = b
    fc = fb
    do iter = 1 , ITMAX

       ! If the new value (f(c)) doesn't bracket
       ! the root with f(b) then adjust it.. 
       if ( sign(1.0,fb) .eq. sign(1.0,fc) ) then
          c  = a
          fc = fa
          d  = b - a
          e  = d
       endif
       if ( abs(fc) .lt. abs(fb) ) then
          a  = b
          b  = c
          c  = a
          fa = fb
          fb = fc
          fc = fa
       endif
       tol1 = 2. * EPS * abs(b) + 0.5 * tol
       xm   = .5 * ( c - b )
       if ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0. ) ) then
          zbrent = b
          return
       endif
       if ( ( abs(e) .ge. tol1 ) .and. ( abs(fa) .gt. abs(fb) ) ) then
          s = fb / fa
          if ( a .eq. c ) then
             p = 2. * xm * s
             q = 1. - s
          else
             q = fa / fc
             r = fb / fc
             p = s * ( 2. * xm * q * ( q - r ) - ( b - a ) * ( r - 1. ) )
             q = ( q - 1. ) * ( r - 1. ) * ( s - 1. )
          endif
          if ( p .gt. 0. ) q = -q
          p = abs( p )
          if ( (2.*p) .lt. min( 3.*xm*q-abs(tol1*q) , abs(e*q) ) ) then
             e = d
             d = p / q
          else
             d = xm
             e = d
          endif
       else
          d = xm
          e = d
       endif
       a  = b
       fa = fb
       if ( abs(d) .gt. tol1 ) then
          b = b + d
       else
          b = b + sign( tol1 , xm )
       endif
       fb = func(b)
    enddo
    write(message,*) "zbrent has exceeded maximum iterations",new_line('x'),&
                     "zbrent was called by: ",trim(called_from)
    call write_log( message , msg_warning , __FILE__ , __LINE__ )
    zbrent = b

  end function zbrent
  !
  !----------------------------------------------------------------------
  !
  ! procedures below are private, ie their use is limited to this module
  !
  !----------------------------------------------------------------------
  !
  real function bilinear_interp( data , corners , location )

    ! Bilinearly interpolates to find value at location based  !
    ! on values in data occuring at corners.  Assumes that     !
    ! data  starting at bottom left and !
    ! moving clockwise,ie 2nd is top left, 3rd is top right.   !

    implicit none

    ! arguments..
    real,dimension(4),  intent(in) :: data     ! the values of the 4 nearest grid pts
    real,dimension(2,2),intent(in) :: corners  ! (1,*) = x1,x2, (2,*) = y1,y2
    real,dimension(2),  intent(in) :: location ! (x,y)

    ! local variables..,
    real :: xfrac, yfrac

    ! check that there is not almost-zero x distance...
    if ( ( corners(1,2) - corners(1,1) ) .lt. tiny ) then
      xfrac = 1
    else
      xfrac = ( location(1) - corners(1,1) ) / ( corners(1,2) - corners(1,1) )
    endif

    ! check that there is not almost-zero y distance...
    if ( ( corners(2,2) - corners(2,1) ) .lt. tiny ) then
      yfrac = 1
    else
      yfrac = ( location(2) - corners(2,1) ) / ( corners(2,2) - corners(2,1) )
    endif

    ! calculate the solution...
    bilinear_interp =     yfrac     *     xfrac     * data(3)  &
                    +     yfrac     * ( 1 - xfrac ) * data(2)  &
                    + ( 1 - yfrac ) *     xfrac     * data(4)  &
                    + ( 1 - yfrac ) * ( 1 - xfrac ) * data(1)

  end function bilinear_interp
  !
  !----------------------------------------------------------------------
  !
  subroutine runge_kutta_check( y , dydx , n , x , h , yout , yerr , derivs )

    ! > subroutine summary? < !

    use gv_scale_declarations, only: nmax

    implicit none

    ! arguments..
    integer,intent(in) :: n
    real,intent(in)    :: h,x,dydx(n),y(n)
    real,intent(out)   :: yerr(n),yout(n)

    ! Interfaces are the correct way to pass procedures as arguments.
    ! (note 'derivs' is actually one of canopy_water_store,soil_water_store,lwp_diff_eqn)
    external :: derivs
    interface
       subroutine derivs( time , y , dydt )
         use gv_scale_declarations, only: nmax
         real,intent(in)   :: y(nmax)
         real,intent(in)   :: time
         real,intent(out)  :: dydt(nmax)
       end subroutine derivs
    end interface

    ! local variables..
    integer        :: i
    real           :: ak2(nmax),ak3(nmax),ak4(nmax),ak5(nmax),ak6(nmax),ytemp(nmax)
    real,parameter :: A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,B32=9./40., &
         B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,B53=-70./27.,B54=35./27.,      &
         B61=1631./55296.,B62=175./512.,B63=575./13824.,B64=44275./110592.,         &
         B65=253./4096.,C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.,        &
         DC1=C1-2825./27648.,DC3=C3-18575./48384.,DC4=C4-13525./55296.,             &
         DC5=-277./14336.,DC6=C6-.25

    ! calculations...
    do i = 1 , n
       ytemp(i) = y(i) + B21 * h * dydx(i)
    enddo
    call derivs( x + A2 * h , ytemp , ak2 )
    do i = 1 , n
       ytemp(i) = y(i) + h * ( B31 * dydx(i) + B32 * ak2(i) )
    enddo
    call derivs( x + A3 * h , ytemp , ak3 )
    do i = 1 , n
       ytemp(i) = y(i) + h * ( B41 * dydx(i) + B42 * ak2(i) + B43 * ak3(i) )
    enddo
    call derivs( x + A4 * h , ytemp , ak4 )
    do i = 1 , n
       ytemp(i) = y(i) + h * ( B51 * dydx(i) + B52 * ak2(i) + B53 * ak3(i) + B54 * ak4(i) )
    enddo
    call derivs( x + A5 * h , ytemp , ak5 )
    do i = 1 , n
       ytemp(i) = y(i) + h * ( B61 * dydx(i) + B62 * ak2(i) + B63 * ak3(i) + B64 * ak4(i) + B65 * ak5(i) )
    enddo
    call derivs( x + A6 * h , ytemp , ak6 )
    do i = 1 , n
       yout(i) = y(i) + h * ( C1 * dydx(i) + C3 * ak3(i) + C4 * ak4(i) + C6 * ak6(i) )
       yerr(i) = h * ( DC1 * dydx(i) + DC3 * ak3(i) + DC4 * ak4(i) + DC5 * ak5(i) + DC6 * ak6(i) )
    enddo

  end subroutine runge_kutta_check
  !
  !----------------------------------------------------------------------
  !
  subroutine runge_kutta_q_step( called_from , y , dydx , n , x , htry , eps , yscal , hdid , hnext , derivs )

    ! > subroutine summary? < !

    use gv_scale_declarations, only: nmax
    use log_tools

    implicit none

    ! arguments..
    character(len=*),intent(in) :: called_from    ! name of procedure calling (used to pass through for errors)
    integer,intent(in)          :: n
    real,intent(in)             :: eps, htry, dydx(n), yscal(n)
    real,intent(inout)          :: x
    real,intent(out)            :: hdid, hnext, y(n)

    ! Interfaces are the correct way to pass procedures as arguments.
    ! (note 'derivs' is actually one of canopy_water_store,soil_water_store,lwp_diff_eqn)
    interface
       subroutine derivs( time , y , dydt )
         use gv_scale_declarations, only: nmax
         real,intent(in)  :: y(nmax)
         real,intent(in)  :: time
         real,intent(out) :: dydt(nmax)
       end subroutine derivs
    end interface

    ! local variables..
    integer        :: i
    real           :: errmax, h, htemp, xnew, yerr(nmax), ytemp(nmax)
    real,parameter :: ERRCON = 1.89e-4, & !
                      PGROW  = -.2,     & !
                      PSHRNK = -.25,    & !
                      SAFETY = 0.9        !

    ! calculations...
    h = htry
1   call runge_kutta_check( y , dydx , n , x , h , ytemp , yerr , derivs )
    errmax = 0.
    do  i = 1 , n
       errmax = max( errmax , abs( yerr(i) / yscal(i) ) )
    enddo
    errmax = errmax / eps
    if ( errmax .gt. 1. ) then
       htemp = SAFETY * h * ( errmax**PSHRNK )
       h = sign( max( abs(htemp) , 0.1*abs(h) ) , h )
       xnew = x + h
       if ( xnew .eq. x ) then
         write(message,*) "stepsize underflow in runge_kutta_q_step",new_line('x'),&
                          "runge_kutta_q_step called from: ",trim(called_from)
         call write_log( message , msg_warning , __FILE__ , __LINE__ )
       endif
       goto 1
    else
       if ( errmax .gt. ERRCON ) then
          hnext = SAFETY * h * ( errmax**PGROW )
       else
          hnext = 5. * h
       endif
       hdid = h
       x = x + h
       do i = 1 , n
          y(i) = ytemp(i)
       enddo

    endif

  end subroutine runge_kutta_q_step
  !
  !----------------------------------------------------------------------
  !
end module math_tools
!
!------------------------------------------------------------------------
!
