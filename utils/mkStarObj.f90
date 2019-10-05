! DATE: July 28, 2019
! compile: gfortran ../quaternion.f90 starMod.f90 mkStarObj.f90 -o mkStarObj
! run: ./mkStarObj
program mkStarObj
use starMod
implicit none

starfilename = "stars.dat"
Nstars = 300

call read_star_file

call write_star_obj


end program mkStarObj
