! compile: gfortran mkRoadCirc.f90 -o mkRoadCirc
program mkPlane
implicit none
integer :: Nx, Ny
integer :: i, j, n, rgb(3), t, Ntotal
real(4) :: dx, dist, rad, rw, x, y
character(50) :: carg

if (iargc().ne.3) then
 write(0,*) "requires, x and y node numbers, grid spacing, and color"
 write(0,*) "Usage:"
 write(0,*) "  mkPlane  Nx Ny dx" ! R G B"
 STOP
endif

call getarg(1, carg); read(carg,*) Nx
call getarg(2, carg); read(carg,*) Ny
call getarg(3, carg); read(carg,*) dx
!call getarg(4, carg); read(carg,*) rgb(1)
!call getarg(5, carg); read(carg,*) rgb(2)
!call getarg(6, carg); read(carg,*) rgb(3)

rad = 0.5*real(min(Nx,Ny))*dx ! max radius
rw = 0.2*real(min(Nx,Ny))*dx! road width ! quarter of min dist
!Nx = 16
!Ny = 10
!dx = 1.0 !distance between nodes
!rgb = (/ 200, 200, 200 /)

write(6,'(A)') "name Circle Road"
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
  x = dx*(real(i-1)-real(Nx)*0.5)
  y = dx*(real(j-1)-real(Ny)*0.5)
  dist = sqrt( x*x + y*y )
  if     (dist .gt. rad)    then; rgb = (/ 0, 200, 0 /) ! green grass
  elseif (dist .gt. rad-0.5*rw+0.5*dx) then; rgb = (/ 50, 50, 50 /) ! dark grey road color
  elseif (dist .gt. rad-0.5*rw-0.5*dx) then; rgb = (/ 255, 255, 0 /) ! yellow line color
  elseif (dist .gt. rad-rw) then; rgb = (/ 50, 50, 50 /) ! dark grey road color
  else;                           rgb = (/ 0, 200, 0 /) ! green grass
  endif
  write(6,*)  n, x, y, 0.0,  rgb(:),  1 ! need last flag to be 1 to get normals
 enddo
enddo
Ntotal = n

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

write(6,'(A,x,i6)') "VertexNormals ",Ntotal
do i = 1, Ntotal
 write(6,*) i, 0.0, 0.0, 1.0
enddo

end program mkPlane
