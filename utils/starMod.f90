module starMod

character(80) :: starfilename
integer :: Nstars

type startype
 integer :: rank
 character(26) :: BayerName
 character(18) :: ProperName
 real :: RA, Dec !right ascension (convert from hours minutes to decimaal hours)
 real :: l, b  !galactic coordinates (degrees)
 character(13) :: spec !spectral type
 real :: visMag, absMag
 real :: para, perr !parallax and error (milli-arcseconds)
 integer :: dist ! distance in lightyear
 real :: cp(3), sp(3)
end type startype

type(startype), dimension(300) :: star


contains
subroutine read_star_file
use quaternion
implicit none
integer :: i, a, b, err
character(128) :: frmt


!Column 1: Approximate rank of the star in order of apparent magnitude.
!Column 2: Bayer name for the star.
!Column 3: Proper name of the star.  Usually a weird arabic name. :)
!Column 4: Right Ascension in hours and minutes for epoch 2000.
!Column 5: Declination in degrees for epoch 2000.
!Column 6: Galactic longitude of the star.
!Column 7: Galactic latitude of the star.
!Column 8: Spectral classification of the main stars in the system.
!Column 9: Apparent visual magnitude of the star.  A letter 'v' means the magnitude varies by
!          more than 0.1.  A letter 'e' means that the star is an eclipsing binary.
!Column 10: Absolute magnitude of the star.  A letter 'v' means the magnitude varies by more
!           than 0.1.
!Column 11: The Hipparcos parallax of the star (milli-arcseconds).
!Column 12: The error in the parallax (milli-arcseconds).
!Column 13: The distance in light years (=3.2616/parallax).

frmt = '(i3,2x,A26,A18,i2,x,i2,x,f5.1,x,f6.1,x,f5.1,2x,A13,f5.2,2x,'
frmt = trim(frmt)//'f5.2,2x,f6.2,x,f4.2,2x,i4)'

open(10, file=trim(starfilename))

 do i = 1, 29
  read(10,*) !scan through the header info
 enddo
 do i = 1, Nstars
  read(10,trim(frmt),iostat=err) star(i)%rank, star(i)%BayerName, &
   & star(i)%ProperName, a, b, star(i)%Dec, star(i)%l, star(i)%b, &
   & star(i)%spec, star(i)%visMag, star(i)%absMag, star(i)%para, star(i)%perr, &
   & star(i)%dist
   if (err.ne.0) then
     write(0,*) "read error ",err," in "//trim(starfilename)//" on line ", i+29
     write(0,frmt) star(i)%rank, star(i)%BayerName, &
   & star(i)%ProperName, a, b, star(i)%Dec, star(i)%l, star(i)%b, &
   & star(i)%spec, star(i)%visMag, star(i)%absMag, star(i)%para, star(i)%perr, &
   & star(i)%dist
     STOP
   endif
   star(i)%RA = real(a) + real(b)/60.0
   star(i)%sp(1) = real(star(i)%dist) !this distance doesn't matter
   star(i)%sp(2) = (star(i)%l-180.0)*3.1415926535/180.0
   star(i)%sp(3) = star(i)%b*3.1415926535/180.0
   call sph2cart( star(i)%sp, star(i)%cp )
 enddo
close(10)

end subroutine read_star_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_star_obj
implicit none
integer :: ic(3), i
real(4), dimension(3) :: sp

write(6,'(A)') "name stars"
write(6,'(A)') "mode point"
write(6,'(A)') "offset 0.0 0.0 0.0"
write(6,'(A)') "velocity 0.0 0.0 0.0"
write(6,'(A)') "orient 0.0 1.0 0.0 0.0"
write(6,'(A)') "spin 0.0 1.0 0.0 0.0"
write(6,'(A)') "points 300 "

do i = 1, Nstars
   ! linear relation,  cut off at mag=0 (looks better as pixels)
   ic(:) = int( 255.0/max(star(i)%visMag+1.0, 1.0) )
   ! True relation,  cut off at mag=0
   !ic(:) = int( 255.0*2.512**(1.0-max(star(i)%visMag+1.0, 1.0)) )
   write(6,*) i, star(i)%cp, ic(:), 0
enddo

write(6,'(A)') "triangles 0 "

end subroutine write_star_obj

end module starMod
