! compile: gfortran gts2obj.f90 -o gts2obj
program gts2obj
implicit none
character(len=80) :: inF, outF, Csc
integer :: Np, Ne, Nt, i
real(4), dimension(3) :: p
integer, allocatable, dimension(:,:) :: e
integer, dimension(3) :: t
integer, dimension(6) :: ii
real, dimension(6) :: tt
real :: scal

if (iargc().ne.3) then
  write(0,*) "need input file gts and outfile obj and scale factor"
  write(0,*) "defaults are hard coded into this program, all color is white"
endif

call getarg(1,inF)
call getarg(2,outF)
call getarg(3,Csc)

read(Csc,*) scal
open(10,file=trim(inF))
open(20,file=trim(outF))

! read header of gts file
read(10,*) Np, Ne, Nt
! write header to obj file
write(20,'(A)') "name gts2obj"
write(20,'(A)') "offset 10.0 0.0 0.0 #point offset to global coordinates"
write(20,'(A)') "velocity 0.0 0.0 0.0 #initial velocity of object"
write(20,'(A)') "orient 0.0 1.0 1.0 1.0 #angle axis notiation, converted to rotor quaternion"
write(20,'(A)') "spin 0.15708 0.0 0.0 1.0 #angle axis notation, angle per centisecond"
write(20,'(A)') "# point list with local coordinates and RGB for each vertex, and smoothing flag"
! loop over points
write(20,*) "points ",Np
do i = 1, Np
  read(10,*) p
  write(20,*) i, p*scal, 255, 255, 255, 0 
enddo

allocate( e(2,Ne) )
!write(20,*) "points ",Np
do i = 1, Ne
  read(10,*) e(:,i)
!  write(20,*) i, p, 255, 255, 255, 0 
enddo

write(20,*) "triangles ",Nt
do i = 1, Nt
  read(10,*) t
  ii(1:2) = e(:,t(1))
  ii(3:4) = e(:,t(2))
  ii(5:6) = e(:,t(3))
  tt = real(ii)
  call quicksort( tt )
  ii = nint(tt)
!  write(20,*) minval(e(:,t(1))), minval(e(:,t(2))), minval(e(:,t(3)))
  write(20,*) ii(1), ii(3), ii(5)
enddo

close(20)
close(10)
contains

! quicksort.f -*-f90-*-
! Author: t-nissie, some tweaks by 1AdAstra1
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
recursive subroutine quicksort(a)
  implicit none
  real :: a(:)
  real x, t
  integer :: first = 1, last
  integer i, j

  last = size(a, 1)
  x = a( (first+last) / 2 )
  i = first
  j = last
  
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  
  if (first < i - 1) call quicksort(a(first : i - 1))
  if (j + 1 < last)  call quicksort(a(j + 1 : last))
end subroutine quicksort

end program gts2obj
