! compile: gfortran mkPerlinPlane.f90 -o mkPerlinPlane
! reoutines converted from: https://thebookofshaders.com/13/
module pland

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function random( v )
implicit none
real(4), dimension(2), intent(in) :: v
real(4) :: random, b
real(4), dimension(2) :: a

!!!!!!!!!! function 1 !!!!!!!!!!!!!
a = (/ 12.9898, 78.233 /)
b = sin(dot_product(v,a))*43758.5453123

!!!!!!!!!! function 2 !!!!!!!!!!!!!
! hash() from https://www.shadertoy.com/view/4dS3Wd
!a = (/ 17.0, 13.0 /)
!b = 1.e4*sin(a(1)*v(1) +v(2)*0.1) *(0.1-abs( sin(v(2)*a(2) +v(1)) ))

!!!!!!!!!! function 3 !!!!!!!!!!!!!
!a = (/ 1.0, 0.5 /)
!b = 1234.5*dot_product(v,a)

random = b -floor(b)

end function random
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function noise( v, q )
implicit none
integer :: seed
real(4), dimension(2), intent(in) :: v
real(4) :: q, noise, a, b, c, d, pi, unt, x1, x2
real(4), dimension(2) :: i, f, u

pi = 3.141596535

unt = 1.0/q  ! unit 
i = floor(v/unt) ! integer of v
f = mod(v,unt)/unt     ! remainder [0:1)
u = 0.5*(1.0-cos(pi*f))  ! Perlin noise?
!u = 1.0-cos(2.0*pi*f)

!a = 1.0; b = 1.0; c = 1.0; d = 1.0
! four corners in 2D of a tile
a = random(i)
b = random(i +(/ 1.0, 0.0 /) )
c = random(i +(/ 0.0, 1.0 /) )
d = random(i +(/ 1.0, 1.0 /) )
!seed = int(dot_product(i,i))
!call random_seed(seed)
!call random_number(a)
!call random_number(b)
!call random_number(c)
!call random_number(d)

!u = f*f*(3.0 -2.0*f)
!u = f
x1 = a*(1.0-u(1)) +b*u(1) ! mix(a,b,u(1))
x2 = c*(1.0-u(1)) +d*u(1) ! mix(c,d,u(1))
noise = x1*(1.0-u(2)) +x2*u(2) ! mix(x1,x2,u(2))

!noise = a*(1.0-u(1)) +b*u(1) +(c-a)*u(2)*(1.0-u(1)) &
!    &   +(d-b)*u(1)*u(2) 

!a = 1.0-(a-0.5)*0.2 !/(seed+3.0)
!b = u(1)-(b-0.5)*0.2 !/(i(1)+12.0)
!c = u(2)-(c-0.5)*0.2  !/(i(2)+12.0)
!noise = a*sin(b*pi*2)*sin(c*pi*2)

end function noise
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function fbm( v, lac, gain, oct )
implicit none
real(4), dimension(2), intent(in) :: v
real(4), dimension(2) :: st
integer, intent(in) :: oct !number of octaves
real(4) :: fbm, val, amp, freq, lac, gain, normk
integer :: i

st = v
val = 0.0
amp = 1.0
freq = 4.0
normk=0.
! loop over octaves
do i = 1, oct
  val = val +amp*noise( st, freq )
  st = st*lac !lacunarity
  freq = freq*lac !lacunarity
  amp = amp*gain !gain
  normk = normk+amp
enddo
fbm = val/normk
!fbm = fbm*fbm*fbm*fbm
!fbm = fbm*fbm*fbm
fbm = fbm*fbm
!fbm = fbm
!fbm = sqrt(fbm)

end function fbm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module pland
!!!!!!!!!!!!!!!!!!!!!!!!!
program mkPlane
use pland
implicit none
integer :: i, j, oct, seed
real(4), dimension(2) :: x
real :: lac, gain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: Nx, Ny
integer :: n, rgb(3), t
real(4) :: dx
character(50) :: carg

if (iargc().ne.9) then
 write(0,*) "requires, x and y node numbers, grid spacing, and color"
 write(0,*) " Lacunarity, Gain, and Octives"
 write(0,*) "Usage:"
 write(0,*) "  mkPerlinPlane  Nx Ny dx R G B  L G O"
 STOP
endif

call getarg(1, carg); read(carg,*) Nx
call getarg(2, carg); read(carg,*) Ny
call getarg(3, carg); read(carg,*) dx
call getarg(4, carg); read(carg,*) rgb(1)
call getarg(5, carg); read(carg,*) rgb(2)
call getarg(6, carg); read(carg,*) rgb(3)

call getarg(7,carg); read(carg,*) lac   !2.0
call getarg(8,carg); read(carg,*) gain  !0.5
call getarg(9,carg); read(carg,*) oct   !8 ?

write(6,'(A)') "name PerlinPlane"
write(6,'(A)') "mode solid"
write(6,'(A)') "offset 0.0 0.0 0.0"
write(6,'(A)') "velocity 0.0 0.0 0.0"
write(6,'(A)') "orient 0.0 1.0 1.0 1.0"
write(6,'(A)') "spin 0.0 0.0 0.0 1.0"

write(6,'(A,x,i6)') "points ",Nx*Ny
! write points
n = 0
do j = 1, Ny
 do i = 1, Nx
  n = n + 1
  ! fractional coordinates
  x = (/ real(i)/real(Nx), real(j)/real(Ny) /)
  write(6,*)  n, dx*(real(i-1)-real(Nx)*0.5), dx*(real(j-1)-real(Ny)*0.5), &
     &  fbm( x, lac, gain, oct ),  rgb(:),  0
 enddo
enddo

write(6,'(A,x,i6)') "triangles ",2*(Nx-1)*(Ny-1)
n = 0; t = 0
do j = 1, Ny
 do i = 1, Nx
  n = n + 1 ! count the node ID, though

  if (j.eq.Ny) cycle
  if (i.eq.Nx) cycle

  ! upper triangle
  t = t + 1 ! triangle count
  write(6,*) n, n+1, n+Nx+1
  ! lower triangle
  t = t + 1 ! triangle count
  write(6,*) n, n+Nx, n+Nx+1

 enddo
enddo

end program mkPlane
