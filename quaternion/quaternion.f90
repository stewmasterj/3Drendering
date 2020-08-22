! vim:fdm=marker
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
module quaternion

! quaternions are 4 element arrays [s, i, j, k]
! vectors are 3 element arrays [x, y, z]

! quatmult(a,b,c) multiply a by b to result in c
! quatinv(a,b)    invert a to produce b
! quatnorm(a,b)   norm of a is b
! quatunit(a,b)   b is the unitvector of a
! quatexp(a,b)    exponential of a is b
! quatlog(a,b)    natural logarithm of a is b
! quatpow(a,b,c)  a raised to the b power is c
! quatdotnorm(a,b,c)  a dot b is c scalar

! geodesicdist(a,b,c) compute geodesic distance from a to b
! quatrotate(a,b,c) rotate pure vector a by unit rotor b, giving c
! quat3rotate(a,b,c) rotate pure vector a(3) by unit rotor b, giving c(3)
! quatrotor(a,b,c) create a rotor from axis a and angle b (radians) to give rotor c.
! abrotor(a,b,c)  put in two vectors a,b to get rotor that rotates a to b.
! rotor2matrix(q,R) create rotation matrix from rotor
! angs2rotor(a,b) Euler angles to rotor
! rotoraxis(q,a)  return angle,axis notation, a, from rotor, q.

! cart2sph(a,b)   cartesian vector a converted to spherical vector b
! sph2cart(a,b)   spherical vector a converted to cartesian vector b
! cross_product(a,b,c) 3D cross product of vector a and b, returns c.

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine quatmult( a, b, c ) !{{{
implicit none
real(4), dimension(4), intent(in)  :: a, b
real(4), dimension(4), intent(out) :: c !c is output
! ab = c
!(a+bi+cj+dk)(e+fi+gj+hk)=
!(ae−bf−cg−dh) + (af+be+ch−dg)i + (ag−bh+ce+df)j + (ah+bg−cf+de)k 
c(1) = (a(1)*b(1)-a(2)*b(2)-a(3)*b(3)-a(4)*b(4)) 
c(2) = (a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3)) !i 
c(3) = (a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2)) !j 
c(4) = (a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1)) !k 
end subroutine quatmult !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine quatinv( a, b ) !{{{
implicit none
real(4), dimension(4), intent(in)  :: a
real(4), dimension(4), intent(out) :: b !b is output
real(4) :: d
!(a+bi+cj+dk)^{-1}=1/(aa+bb+cc+dd)(a-bi-cj-dk)
d = a(1)*a(1)+a(2)*a(2)+a(3)*a(3)+a(4)*a(4)
! b is the conugate of a i.e.:   b = a*
b(1) =  a(1)
b(2) = -a(2)
b(3) = -a(3)
b(4) = -a(4)
b = b/d
end subroutine quatinv !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine quatnorm( a, b ) !{{{
implicit none
real(4), dimension(4), intent(in)  :: a
real(4), intent(out) :: b !b is output
!|q| = sqrt{qq*} = sqrt{q*q} = sqrt{aa+bb+cc+dd}
b = sqrt( a(1)*a(1)+a(2)*a(2)+a(3)*a(3)+a(4)*a(4) )
end subroutine quatnorm !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine quatunit( a, b ) !{{{
implicit none
real(4), dimension(4), intent(in)  :: a
real(4), dimension(4), intent(out) :: b !b is output
real(4) :: d
!call quatnorm(a, d)
d = sqrt( a(1)*a(1)+a(2)*a(2)+a(3)*a(3)+a(4)*a(4) )
b = a/d
end subroutine quatunit !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine quatexp( a, b ) !{{{
implicit none
real(4), dimension(4), intent(in)  :: a
real(4), dimension(4)  :: aa
real(4), dimension(4), intent(out) :: b !b is output
real(4) :: d
!exp(q) = sum{n=0,infty}{q^{n}/n!} = exp{a}(cos|v|+v/|v|sin|v|)
aa = a
aa(1) = 0.0
!call quatnorm(aa, d)
d = sqrt( a(1)*a(1)+a(2)*a(2)+a(3)*a(3)+a(4)*a(4) )
b = aa/d*sin(d)
b(1) = cos(d)
b = exp(a(1))*b
end subroutine quatexp !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine quatlog( a, b ) !{{{
implicit none
real(4), dimension(4), intent(in)  :: a
real(4), dimension(4), intent(out) :: b !b is output
real(4) :: qn, vn, d
! ln(q) = ln(|q|) + v/|v|acos(a/|q|)
!call quatnorm(aa, d)
d = a(2)*a(2)+a(3)*a(3)+a(4)*a(4)
qn = sqrt( a(1)*a(1) + d )
vn = sqrt( d )
b(2:4) = a(2:4)/vn*acos(a(1)/qn)
b(1) = log(qn)
end subroutine quatlog !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine quatpow( a, b, c ) !{{{
implicit none
real(4), dimension(4), intent(in)  :: a
real(4), dimension(4)  :: qn, ql
real(4), intent(in) :: b
real(4), dimension(4), intent(out) :: c !c is output
!https://www.mathworks.com/help/aerotbx/ug/quatpower.html
!qp=quatpower(q,pow) calculates q to the power of pow for a normalized quaternion, q.
!qp = quatpower(quatnormalize([1.0 0 1.0 0]),2)
!qt=exp(t⋅log(q)),with t∈R.

call quatunit(a,qn) !make sure input is normalized (unitvector)
call quatlog(qn,ql)
call quatexp(b*ql,c)

end subroutine quatpow !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine quatdotnorm( a, b, c ) !{{{
implicit none
real(4), dimension(4), intent(in)  :: a, b
real(4), intent(out) :: c !c is output
real(4) :: x, y, z
!cosθ=s1s2+x1x2+y1y2+z1z2/|q1||q2|
x = sum(a*b)
y = sqrt( a(1)*a(1)+a(2)*a(2)+a(3)*a(3)+a(4)*a(4) )
z = sqrt( b(1)*b(1)+b(2)*b(2)+b(3)*b(3)+b(4)*b(4) )
c = x/(y*z)
end subroutine quatdotnorm !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine geodesicdist( a, b, c ) !{{{
implicit none
real(4), dimension(4), intent(in)  :: a, b
real(4), intent(out) :: c !c is output
c = acos(2.0*sum(a*b)**2 -1.0 )
end subroutine geodesicdist !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine quatrotate( a, b, c ) !{{{
implicit none
real(4), dimension(4), intent(in)  :: a, b
real(4), dimension(4)  :: e, bi
real(4), dimension(4), intent(out) :: c !c is output
! b must be a unit vector rotor 
! a must have no scalar component a(1)=0
! c = a' = bab^-1
call quatinv( b, bi )
! c = b*a*bi = e*bi
call quatmult( b, a, e )
call quatmult( e, bi, c )
end subroutine quatrotate !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine quat3rotate( a3, b, c3 ) !{{{
implicit none
real(4), dimension(3), intent(in)  :: a3
real(4), dimension(4), intent(in)  :: b
real(4), dimension(4)  :: a, c, e, bi
real(4), dimension(3), intent(out) :: c3 !c is output
! b must be a unit vector rotor 
! a must have no scalar component a(1)=0
! c = a' = bab^-1
call quatinv( b, bi )
! c = b*a*bi = e*bi
a(1) = 0.0; a(2:4) = a3(1:3)
call quatmult( b, a, e )
call quatmult( e, bi, c )
c3(1:3) = c(2:4)
end subroutine quat3rotate !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine quatrotor( a, b, c ) !{{{
implicit none
real(4), dimension(4), intent(in)  :: a
real(4), intent(in) :: b
real(4) :: n
real(4), dimension(4), intent(out) :: c !c is output
! a must be a unit vector defining the axis of rotation
! b is the angle about the axis in radians
!q=[cos(θ/2),sin(θ/2)v^]
c(1) = cos(0.5*b)
n = sin(0.5*b)
c(2:4) = n*a(2:4)
end subroutine quatrotor !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine ab2rotor( a, b, c ) !{{{
implicit none
real(4), dimension(3), intent(in) :: a, b
real(4), dimension(4), intent(out) :: c
real(4), dimension(4) :: q
real(4), dimension(3) :: ax
real(4) :: ang, da, db

da = SQRT(dot_product(a,a))
db = SQRT(dot_product(b,b))
ang = ACOS(dot_product(a,b)/(da*db))
!if (ang.gt.3.141593) ang = ang-3.1415926535

call cross_product(a/da,b/db,ax) ! get the rotation axis
da = SQRT(dot_product(ax,ax))
if (abs(da).lt.1.e-10) then
 c(2:4) = a; c(1) = 0.0
 return
endif

q=0.0
q(2:4) = ax/da !set rotation axis
call quatrotor( q, -ang, c )

end subroutine ab2rotor !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine rotor2matrix( q, R ) !{{{
implicit none
real(4), dimension(4), intent(in) :: q
real(4), dimension(3,3), intent(out) :: R
real(4) :: s
call quatnorm( q, s )
s = 1.0/(s*s)
R(1,:) = (/ 1.0 - 2.0*s*(q(3)*q(3) +q(4)*q(4)), 2.0*s*(q(2)*q(3)-q(4)*q(1)), 2.0*s*(q(2)*q(4)+q(3)*q(1)) /)
R(2,:) = (/ 2.0*s*(q(2)*q(3)+q(4)*q(1)), 1.0-2.0*s*(q(2)*q(2)+q(4)*q(4)), 2.0*s*(q(3)*q(4)-q(2)*q(1)) /)
R(3,:) = (/ 2.0*s*(q(2)*q(4)-q(3)*q(1)), 2.0*s*(q(3)*q(4)+q(2)*q(1)), 1.0-2.0*s*(q(2)*q(2)+q(3)*q(3)) /)
end subroutine rotor2matrix !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine angs2rotor( a, q ) !{{{
implicit none
real(4), dimension(3), intent(in) :: a !Euler angles,  yaw (Z), pitch (Y), roll (X)
real(4), dimension(4), intent(out) :: q ! quaternion rotor
real(4) :: cy, sy, cp, sp, cr, sr
! from : https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
cy = cos(a(3) * 0.5)
sy = sin(a(3) * 0.5)
cp = cos(a(2) * 0.5)
sp = sin(a(2) * 0.5)
cr = cos(a(1) * 0.5)
sr = sin(a(1) * 0.5)

q(1) = cy * cp * cr + sy * sp * sr
q(2) = cy * cp * sr - sy * sp * cr
q(3) = sy * cp * sr + cy * sp * cr
q(4) = sy * cp * cr - cy * sp * sr
end subroutine angs2rotor !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine rotoraxis( q, a ) !{{{
implicit none
real(4), dimension(4), intent(in) :: q
real(4), dimension(4), intent(out) :: a
real(4) :: n
a=0.0
a(1) = acos(q(1))
n = sin(a(1))
if (abs(a(1)).lt.1.e-24) then
 n=1.0 !zero angle, so don't fuck with the axis
 a(2:4) = 1.0 !just give some axis, derp, zero angle so doesn't matter
else
  a(2:4) = q(2:4)/n !axis of rotation
endif
a(1) = a(1)*2.0   !angle of rotation in radian
! wiki just has this: https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
!n = sqrt(sum(q(2:4)**2))
!a(2:4) = q(2:4)/n
!a(1) = 2.0*atan2(n,q(1))
end subroutine rotoraxis !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine cart2sph( a, b ) !{{{
implicit none
real(4), dimension(3), intent(in)  :: a
real(4), dimension(3), intent(out) :: b
! a is (x, y, z)
! b is (r, theta, phi)
b(1) = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
b(2) = atan2(sqrt(a(1)*a(1)+a(2)*a(2)), a(3))
b(3) = atan2(a(2), a(1))
end subroutine cart2sph !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine sph2cart( a, b ) !{{{
implicit none
real(4), dimension(3), intent(in)  :: a
real(4), dimension(3), intent(out) :: b
! a is (r, theta, phi)
! b is (x, y, z)
!b(1) = a(1)*sin(a(3))*cos(a(2))
!b(2) = a(1)*sin(a(3))*sin(a(3))
!b(3) = a(1)*cos(a(3))
b(1) = a(1)*sin(a(2))*cos(a(3))
b(2) = a(1)*sin(a(2))*sin(a(3))
b(3) = a(1)*cos(a(2))
end subroutine sph2cart !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine cross_product(a,b,c) !{{{
implicit none
real(4), dimension(3), intent(in)  :: a, b
real(4), dimension(3), intent(out) :: c
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
end subroutine !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
end module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
