! compile: gfortran quaternion.f90 quatTest.f90 -o quatTest
program quatTest
use quaternion
implicit none
real(4), dimension(4) :: a, au, q, p, pp, q2, qt, qi, at, qa
!real(4), dimension(4) :: qi, e
real(4) :: ang
! for cart2sph
real(4), dimension(3) :: c, s

!call quatrotate( (/1.0,1.0,1.0,1.0/), (/.99,0.0,0.13,0.0/), pp )
p = (/0.0,1.0,0.0,0.0/) !pp
a = p
write(0,*) "X-axis: ",p(2:4)
! rotate around X-axis 90 degrees
call quatrotor( (/0.0,1.0,0.0,0.0/), 1.5708, q )
q2 = q; q2(3:4) = -q(3:4)
qt = q
! "rotate" the X-axis vector with this rotor
call quatrotate( p, q, pp )
call quatrotate( a, q2, au )
at = q2
write(0,*) "90-X: ",pp(2:4), au(2:4)
! rotate around Y-axis -90 degrees
call quatrotor( (/0.0,0.0,1.0,0.0/), -1.5708, q )
call quatmult( at, q, q2 )
at = q2
q2 = q; q2(3:4) = -q(3:4)
! "rotate" the X-axis vector with this rotor
call quatrotate( pp, q, p )
call quatrotate( au, q2, a )
write(0,*) "-90-Y: ",p(2:4), a(2:4)
! is it the same for a total rotor?
call quatmult( qt, q, q2 )
qt = q2
p = (/0.0,1.0,0.0,0.0/) !pp
call quatrotate( p, qt, pp )
call quatrotate( p, at, au )
call quatinv( qt, qi )
write(0,*) " rotors:", qt, at, qi
call quatrotate( p, qi, a )
write(0,*) " total quat: ", pp(2:4), au(2:4), a(2:4)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(0,*)
p = (/0.0,1.0,0.0,0.0/) !pp
write(0,*) "X-axis: ",p(2:4)
! rotate around Y-axis -90 degrees
call quatrotor( (/0.0,0.0,0.0,1.0/), -1.5708, q )
q2 = q; q2(3:4) = -q(3:4)
! "rotate" the X-axis vector with this rotor
call quatrotate( p, q, pp )
call quatrotate( p, q2, au )
call quatinv( q, qi )
write(0,*) " rotors:", q, q2, qi
call quatrotate( p, qi, a )
write(0,*) "-90-Z: ", pp(2:4), au(2:4), a(2:4)

!call quatrotor( (/0.0,1.0,1.0,1.0/), 0.0, pp )
!write(0,*) "rotor unit:",pp
!!call quatmult( (/.99,0.0,0.13,0.0/), (/.99,0.0,0.13,0.0/), p )
!call quatmult( pp, (/.99,0.0,0.13,0.0/), p )
!call quatrotate( (/1.0,1.0,1.0,1.0/), p, pp )
!write(0,*) pp

write(6,*) "specify a point in 3D space subject to rotation"
write(6,'(A)',advance="no") "(x,y,z) = "
read(5,*) p(2), p(3), p(4)
p(1) = 0.0

write(6,*) "define rotation axis vector"
write(6,'(A)',advance="no") "(x,y,z) = "
read(5,*) a(2), a(3), a(4)
a(1) = 0.0

write(6,*) "define rotation angle (degrees)"
write(6,'(A)',advance="no") "angle = "
read(5,*) ang

ang = ang*3.1415926535/180.0 !convert to radians
call quatunit( a, au ) !make this a unit vector

call quatrotor( a, ang, q )
write(6,*) "using rotor quaternion:", q
call rotoraxis( q, qa )
write(6,*) "recovered angle axis notation:",qa
!!!! Built in !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
call quatrotate( p, q, pp )
!!!! routines !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!call quatinv( q, qi )
!write(6,*) "using inverse rotor quaternion:", qi
!call quatmult( q, p, e )
!write(6,*) "q*p=", e
!call quatmult( e, qi, pp )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
write(6,*) "new position after rotation:", pp(2:4)

write(6,*) "orientations give what spherical angle from atan2()?"
write(6,*) "(x,y,z) (r,the,phi)"
c = (/ 1.0, 0.001, 0.001 /)
call cart2sph(c, s)
write(6,*) nint(c), s
c = (/ -1.0, -0.001, -0.001 /)
call cart2sph(c, s)
write(6,*) nint(c), s
c = (/ -1.0, 0.001, -0.001 /)
call cart2sph(c, s)
write(6,*) nint(c), s
c = (/ 0.001, 1.0, 0.001 /)
call cart2sph(c, s)
write(6,*) nint(c), s
c = (/ -0.001, -1.0, -0.001 /)
call cart2sph(c, s)
write(6,*) nint(c), s
c = (/ 0.001, 0.001, 1.0 /)
call cart2sph(c, s)
write(6,*) nint(c), s
c = (/ -0.001, -0.001, -1.0 /)
call cart2sph(c, s)
write(6,*) nint(c), s

write(6,*) "theta is angle down from Z axis"
write(6,*) "X axis is phi=0.0 "
write(6,*) "-X axis is phi=-pi or pi "
write(6,*) " coordinates in local space craft, make sense as"
write(6,*) "   +X=forward, -X backwards"
write(6,*) "   +Z theta=0 UP, -Z theta=pi DOWN"
write(6,*) "   +Y phi>0 LEFT port, -Y phi<0 RIGHT starboard"
write(6,*)
write(6,*) " Phi could be offset to start at 0 and go to 2pi"
write(6,*) " Human field of view is 210 deg Horizontal and 150 vertical. 1.4 aspectratio."
write(6,*) " games commonly use 100 degrees"

end program quattest
