! Draw some stuff to a spherical background
! Ross J. Stewart
! September 21, 2019
! vim:fdm=marker
! compile: gfortran ../../../quaternion.f90 drawBG.f90 -o drawBG 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
module subobjects !{{{
! this type should help you to formulate the subcomponents
! they will all be merged at the end when writing the file
type subobj
 integer :: nn, nt
 real(4), dimension(:,:), allocatable :: n !nodes
 integer, dimension(:,:), allocatable :: c !node color
 integer, dimension(:,:), allocatable :: t !triangles
end type subobj
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define your subobjects here:

type(subobj) :: px,nx, py,ny, pz,nz
type(subobj) :: rlz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
subroutine soinit( so, nn, nt ) !{{{
implicit none
type(subobj), intent(inout) :: so
integer, intent(in) :: nn, nt

so%nn = nn
so%nt = nt
allocate( so%n(3,nn), so%c(3,nn), so%t(3,nt) )

end subroutine soinit !}}}
end module subobjects !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
program drawBG
use subobjects
use quaternion
implicit none
real(4), dimension(3) :: x, c
real(4), dimension(4) :: qr
real(4) :: PI, r, a, o, th, ph
integer :: i, j, is, ni

PI = ACOS(-1.0)
a = PI*0.01 ! half letter angle
o = PI*0.001 ! half offset angle for letter trianlges, like serif lengths
r = 100.0 ! distance to points on sphere, i.e. radius


!objects are in local cartesian coordinates, we'll draw this from a spherical
! coordinate system (r, theta, phi)
! define Giant X labels {{{
call soinit( px, 9, 4 ) ! 9 odes, 4 triangles for X shape
th = 0.5*PI
ph = 0.0
x = (/ r, th, ph /) ! center coordinate
call sph2cart( x, px%n(:,1) )

x = (/ r, th+a-o, ph-a-o /) ! top left left
call sph2cart( x, px%n(:,2) )
x = (/ r, th+a+o, ph-a+o /) ! top left top
call sph2cart( x, px%n(:,3) )

x = (/ r, th+a-o, ph+a+o /) ! top right right
call sph2cart( x, px%n(:,4) )
x = (/ r, th+a+o, ph+a-o /) ! top right top
call sph2cart( x, px%n(:,5) )

x = (/ r, th-a+o, ph-a-o /) ! bottom left left
call sph2cart( x, px%n(:,6) )
x = (/ r, th-a-o, ph-a+o /) ! bottom left bottom
call sph2cart( x, px%n(:,7) )

x = (/ r, th-a+o, ph+a+o /) ! bottom right right
call sph2cart( x, px%n(:,8) )
x = (/ r, th-a-o, ph+a-o /) ! bottom right bottom
call sph2cart( x, px%n(:,9) )

! define the triangles
px%t(:,1) = (/ 1, 2, 3 /)
px%t(:,2) = (/ 1, 4, 5 /)
px%t(:,3) = (/ 1, 6, 7 /)
px%t(:,4) = (/ 1, 8, 9 /)

c = (/ 0, 120, 0 /) ! green
do i = 1, 9
  px%c(:,i) = c ! set the colors for all nodes the same
enddo

call soinit( nx, 9, 4 )
nx = px ! set the opposite
c = (/ 120, 0, 0 /) ! red
do i = 1, 9
  nx%n(1,i) = -nx%n(1,i) !negate the x component
  nx%c(:,i) = c ! set the colors for all nodes the same
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!}}}

! define Giant Y labels {{{
call soinit( py, 7, 3 ) ! 7 nodes, 3 triangles for Y shape
th = PI*0.5    !theta
ph = PI*0.5 !phi
x = (/ r, th, ph /) ! center coordinate
call sph2cart( x, py%n(:,1) )

x = (/ r, th+a-o, ph-a-o /) ! top left left
call sph2cart( x, py%n(:,2) )
x = (/ r, th+a+o, ph-a+o /) ! top left top
call sph2cart( x, py%n(:,3) )

x = (/ r, th+a-o, ph+a+o /) ! top right right
call sph2cart( x, py%n(:,4) )
x = (/ r, th+a+o, ph+a-o /) ! top right top
call sph2cart( x, py%n(:,5) )

x = (/ r, th-a*1.414, ph-o /) ! bottom left left
call sph2cart( x, py%n(:,6) )
x = (/ r, th-a*1.414, ph+o /) ! bottom left bottom
call sph2cart( x, py%n(:,7) )

! define the triangles
py%t(:,1) = (/ 1, 2, 3 /)
py%t(:,2) = (/ 1, 4, 5 /)
py%t(:,3) = (/ 1, 6, 7 /)

c = (/ 0, 120, 0 /) ! green
do i = 1, 7
  py%c(:,i) = c ! set the colors for all nodes the same
enddo

call soinit( ny, 7, 3 )
ny = py ! set the opposite
c = (/ 120, 0, 0 /) ! red
do i = 1, 7
  ny%n(2,i) = -ny%n(2,i) !negate the y component
  ny%c(:,i) = c ! set the colors for all nodes the same
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!}}}

! define Giant Z labels {{{
call soinit( pz, 7, 3 ) ! 7 nodes, 3 triangles for Y shape
th = 0.0  !theta
ph = 0.0     !phi
x = (/ r, th-o, ph+0.25*PI /) ! center coordinate
call sph2cart( x, pz%n(:,1) )

x = (/ r, th+a, ph-0.25*PI /) ! top left 
call sph2cart( x, pz%n(:,2) )
x = (/ r, th+a-2*o, ph /) ! top middle
call sph2cart( x, pz%n(:,3) )
x = (/ r, th+a, ph+0.25*PI /) ! top right
call sph2cart( x, pz%n(:,4) )

x = (/ r, th-a, ph-0.25*PI /) ! bottom left 
call sph2cart( x, pz%n(:,5) )
x = (/ r, th-a+2*o, ph /) ! bottom middle
call sph2cart( x, pz%n(:,6) )
x = (/ r, th-a, ph+0.25*PI /) ! bottom right
call sph2cart( x, pz%n(:,7) )

! define the triangles
pz%t(:,1) = (/ 2, 3, 4 /) !top
pz%t(:,2) = (/ 2, 1, 5 /) !middle
pz%t(:,3) = (/ 5, 6, 7 /) !bottom

c = (/ 0, 120, 0 /) ! green
do i = 1, 7
  pz%c(:,i) = c ! set the colors for all nodes the same
enddo

call soinit( nz, 7, 3 )
nz = pz ! set the opposite
c = (/ 120, 0, 0 /) ! red
do i = 1, 7
  nz%n(3,i) = -nz%n(3,i) !negate the y component
  nz%c(:,i) = c ! set the colors for all nodes the same
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!}}}

!define some rules for measurement
call soinit( rlz, 120, 40 ) !{{{
is = 0
ni = 0
do i = 1, 7 !theta
 th = real(i)*PI*0.125
 do j = 1, 4 !phi
  ph = real(j*2-5)*PI*0.25
  ni = ni + 1
  is = is + 1
  x = (/ r, th, ph /)  !center
  call sph2cart( x, rlz%n(:,is) )
  c = 60
  if (rlz%n(1,is).gt.0.0) c(1) = 120
  if (rlz%n(2,is).gt.0.0) c(2) = 120
  if (rlz%n(3,is).gt.0.0) c(3) = 120
  rlz%c(:,is) = c
  is = is + 1
  x = (/ r, th, ph-a /) !left
  call sph2cart( x, rlz%n(:,is) )
  rlz%c(:,is) = c
  is = is + 1
  x = (/ r, th, ph+a /) !right
  call sph2cart( x, rlz%n(:,is) )
  rlz%c(:,is) = c
  ! define triangle
  rlz%t(:,ni) = (/ is-1, is-2, is /)
 enddo
enddo
th = 0.5*PI
do j = 1, 15
 if (j.eq.4.or.j.eq.8.or.j.eq.12) cycle
 ni = ni + 1
 is = is + 1
 ph = real(j-8)*PI*0.125
 x = (/ r, th, ph /)  !center
 call sph2cart( x, rlz%n(:,is) )
 c = 60; c(3) = 0
 if (rlz%n(1,is).gt.0.0) c(1) = 120
 if (rlz%n(2,is).gt.0.0) c(2) = 120
 rlz%c(:,is) = c
 is = is + 1
 x = (/ r, th-a, ph /) !left
 call sph2cart( x, rlz%n(:,is) )
 rlz%c(:,is) = c
 is = is + 1
 x = (/ r, th+a, ph /) !right
 call sph2cart( x, rlz%n(:,is) )
 rlz%c(:,is) = c
 ! define triangle
 rlz%t(:,ni) = (/ is-1, is-2, is /)
enddo !}}}
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write all subobjects with correct offset IDs
is = 0

write(6,'(A)') "name drawBG"
write(6,'(A)') "mode wire"
write(6,'(A)') "offset 0.0 0.0 0.0 #point offset to global coordinates"
write(6,'(A)') "velocity 0.0 0.0 0.0 #initial velocity of object"
write(6,'(A)') "orient 0.0 1.0 1.0 1.0 #angle axis notation, converted to rotor quaternion"
write(6,'(A)') "spin 0.0 0.0 0.0 1.0 #angle axis notation, angle per centisecond"
write(6,'(A)') "# point list with local coordinates and RGB for each vertex"
is = px%nn + py%nn + pz%nn
is = is*2.0 + rlz%nn
write(6,'(A,x,i4)') "points",is !{{{
! write the px nodes
is = 0 ! node offset
do i = 1, px%nn
 write(6,'(i4,x,3f12.3,x,4i4)') i+is, px%n(:,i), px%c(:,i), 0
enddo
! write the nx nodes
is = px%nn ! node offset
do i = 1, nx%nn
 write(6,'(i4,x,3f12.3,x,4i4)') i+is, nx%n(:,i), nx%c(:,i), 0
enddo
! write the py nodes
is = is + nx%nn ! node offset
do i = 1, py%nn
 write(6,'(i4,x,3f12.3,x,4i4)') i+is, py%n(:,i), py%c(:,i), 0
enddo
! write the ny nodes
is = is + py%nn ! node offset
do i = 1, ny%nn
 write(6,'(i4,x,3f12.3,x,4i4)') i+is, ny%n(:,i), ny%c(:,i), 0
enddo
! write the pz nodes
is = is + ny%nn ! node offset
do i = 1, pz%nn
 write(6,'(i4,x,3f12.3,x,4i4)') i+is, pz%n(:,i), pz%c(:,i), 0
enddo
! write the nz nodes
is = is + pz%nn ! node offset
do i = 1, nz%nn
 write(6,'(i4,x,3f12.3,x,4i4)') i+is, nz%n(:,i), nz%c(:,i), 0
enddo 
! write the rlz nodes
is = is + nz%nn ! node offset
do i = 1, rlz%nn
 write(6,'(i4,x,3f12.3,x,4i4)') i+is, rlz%n(:,i), rlz%c(:,i), 0
enddo !}}}
is = px%nt + py%nt + pz%nt
is = is*2.0 + rlz%nt
write(6,'(A,x,i4)') "triangles",is !{{{
! write px triangles
is = 0 !index offset
do i = 1, px%nt
 write(6,'(3i5)') px%t(:,i)+is
enddo
! write nx triangles
is = px%nn !index offset
do i = 1, nx%nt
 write(6,'(3i5)') nx%t(:,i)+is
enddo
! write py triangles
is = is + nx%nn !index offset
do i = 1, py%nt
 write(6,'(3i5)') py%t(:,i)+is
enddo
! write ny triangles
is = is + py%nn !index offset
do i = 1, ny%nt
 write(6,'(3i5)') ny%t(:,i)+is
enddo
! write pz triangles
is = is + ny%nn !index offset
do i = 1, pz%nt
 write(6,'(3i5)') pz%t(:,i)+is
enddo
! write nz triangles
is = is + pz%nn !index offset
do i = 1, nz%nt
 write(6,'(3i5)') nz%t(:,i)+is
enddo 
! write rlz triangles
is = is + nz%nn !index offset
do i = 1, rlz%nt
 write(6,'(3i5)') rlz%t(:,i)+is
enddo !}}}

end program drawBG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
