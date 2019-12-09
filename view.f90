! vim:fdm=marker
! Author: Ross J. Stewart
! DATE: December 2018 - October 2019
! compile: gcc -c ../../fcurses/kbhit.c
! compile: gfortran -O -fbounds-check -march=native -c ../../quaternion.f90 ../../fcurses/fcurses.f90 ../../fbMod/fbMod2.f90 ../../stringParseMods/lineParse.f90 shareMods.f90 renderMod.f90 ../../ArgsMod/ArgsMod.f90 view3.f90 
! compile: gfortran -lc -o view3 kbhit.o quaternion.o fcurses.o fbMod2.o lineParse.o shareMods.o renderMod.o ArgsMod.o view3.o -O -march=native
! clean: rm view3 kbhit.o quaternion.o starMod.o fcurses.o fbMod2.o lineParse.o shareMods.o renderMod.o ArgsMod.o view3.o 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
module screenstuff !{{{
use renderMod
real(4), dimension(3) :: axisX, axisY, axisZ !
real(4), dimension(3) :: caX, caY, caZ, caXs !local camera axes
end module screenstuff !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
program view !{{{
! external modules
use iso_c_binding !for usleep
use quaternion
use fbMod
use fcurses
use renderMod
! internal modules
use screenstuff
use controlMod
use options
!use starMod
implicit none
interface
 subroutine usleep(useconds) bind(C)
  use iso_c_binding !, only : usleep
  implicit none
  integer(c_int32_t), value :: useconds
 end subroutine
end interface
real(4), dimension(4) :: qrl, qrr, qru, qrd, qrw, qrc, qr, qtmp2
real(4), dimension(3) :: vt, normal
integer ::  t, err, i
integer, dimension(1) :: mj
!!! for fcurses
character, dimension(5) :: ch
character(5) :: Fc
character(120) :: line
real(4) :: start_time, end_time

call CPU_TIME(start_time)
! Default values !!!!!!{{{
call getOptions
! load init script commands
call loadScreenData( 10, FILscreen )

axisX = (/ 1.0, 0.0, 0.0 /) 
axisY = (/ 0.0, 1.0, 0.0 /) 
axisZ = (/ 0.0, 0.0, 1.0 /) 
!local camera axes
caX = (/ 1.0, 0.0, 0.0 /) 
caY = (/ 0.0, 1.0, 0.0 /) 
caZ = (/ 0.0, 0.0, 1.0 /) 
caXs = caX

if (trim(FILobj).ne."NULL") then
 Nobjects = Nobjects + 1
 call loadObject( 12, trim(FILobj), o(Nobjects) )
endif

! predefine quaternion rotors
! for left right up and down.
call quatrotor( (/ 0.0, 1.0, 0.0, 0.0 /), -scr%dtheta, qrl )
call quatrotor( (/ 0.0, 1.0, 0.0, 0.0 /),  scr%dtheta, qrr )
call quatrotor( (/ 0.0, 0.0, 1.0, 0.0 /), -scr%dtheta, qru )
call quatrotor( (/ 0.0, 0.0, 1.0, 0.0 /),  scr%dtheta, qrd )
call quatrotor( (/ 0.0, 0.0, 0.0, 1.0 /),  scr%dtheta, qrw )
call quatrotor( (/ 0.0, 0.0, 0.0, 1.0 /), -scr%dtheta, qrc )
!!!!!!!!!!!!!!!!!!!!!!!!}}}

! init the screen
call init_screen( scr%tmpdir )
call cls
! open and set framebuffer. write functins to fb%pxbuff before writting to device.
call fb%fbinit(10,scr%fbpath, scr%width, scr%height, scr%line, .true.)

! call one step dynamics to get global vectors uptodate
call objectDynamics

call fb%clear !clear the fb%pxbuff and fb%zbuff
if (.not.allocated(fb%zbuff)) write(0,*) "zbuffer not allocated"
call drawscreen

! can stop program here to get one snapshot of initial conditions before interactions
if (.not.scr%interactive) then
  if (scr%dumpName.ne."_") then
   call fb%save(trim(scr%dumpName),2) !dump the screen to PPM6 file
  endif
  write(6,*)  !for somereason I need this to exit
  call kill_screen( scr%tmpdir )
  call fb%close
  STOP
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

t=0
qr = (/ 1.0, 0.0, 0.0, 0.0 /) !basically a null rotation rotor
vt = 0.0 
do !main do loop for capturing keys
 call usleep(nint(scr%dt*1e6)) !save CPU cycle time by polling keyboard 100 times per second
 t = t + 1
 rotate = .false.
 translate = .false.
 !!! Key capture !!!!!! {{{
 if (kbhit().eq.1) then
  call getkey(ch)
  !!! Interpretation of ctrl chars and such
  select case (ch(1))
  case (achar(27)) !escape sequences
   if (ch(2).eq.'[') then     !usualy begin with '['
     Fc=ch(3)//ch(4)//ch(5)    !curser keys work
     ! rotation motions
     if (Fc.eq.'D  '    ) then; qr=qrl; rotate=.true. !LEFT
     elseif (Fc.eq.'C  ') then; qr=qrr; rotate=.true. !RIGHT
     elseif (Fc.eq.'A  ') then; qr=qru; rotate=.true. !UP
     elseif (Fc.eq.'B  ') then; qr=qrd; rotate=.true. !DOWN
     elseif (Fc.eq.'21~') then   !F10 Quit Exit
      write(6,*)  !for somereason I need this to exit
      call kill_screen( scr%tmpdir )
      call fb%close
      call CPU_TIME(end_time)
      write(6,*) "execution time: ",end_time-start_time
      STOP
     endif
   endif
  case ('q') ; qr=qrw; rotate=.true.   !left rudder, widdershins
  case ('e') ; qr=qrc; rotate=.true.   !right rudder, clockwise
  ! translation motion
  case ('W') ; vt=caZ; translate=.true.   !move up
  case ('S') ; vt=-caZ; translate=.true.   !move down
  case ('w') ; vt=caX; translate=.true.   !move forward
  case ('s') ; vt=-caX; translate=.true.   !move backward
  case ('a') ; vt=-caY; translate=.true.   !straf left
  case ('d') ; vt=caY; translate=.true.   !straf right
  case ('r') ; if (record) then; record=.false.; else; record = .true.; endif
  case ('p') ; scrot=.true.
  ! select node or triangle
  case (' ') ; ! fire something from ship
   if (scr%centerObj.gt.0) then !check if object under aim
    ! let's do something to the object, like change the node color
    call breakObject(scr%centerObj)
   endif
  case (':')  ! command mode
   call tput( ":", 1, lines ) 
   call fancygetrawline( line, achar(13), .true. ) !echo the keys
   call interpretLine( line, err ) !interpreter in "renderMod.f90"
   if (err.eq.1) then !command "exit" and "quit" terminate program
    write(6,*)  !for somereason I need this to exit
    call kill_screen( scr%tmpdir )
    call fb%close
    call CPU_TIME(end_time)
    write(6,*) "execution time: ",end_time-start_time
    STOP
   endif
  case ('Q')    !Quit Exit
   write(6,*)  !for somereason I need this to exit
   call kill_screen( scr%tmpdir )
   call fb%close
   call CPU_TIME(end_time)
   write(6,*) "execution time: ",end_time-start_time
   STOP
  case default !if key is unknown then just print it to screen, cuz why not?
   write(6,'(a1)',advance='no') ch(1)
  end select

  if (rotate) then
   if (impulseControl) then
    call quatmult( cam%spin, qr, qtmp2 )
    cam%spin = qtmp2
   else; cam%spin = qr
   endif
  endif 

  if (translate) then
   if (impulseControl) then
    cam%velocity = cam%velocity +vt
   else; cam%velocity = vt*scr%ds
   endif
   if (Ofollow.gt.0) then ! project velocity onto plane defined by normal
    ! closest vertex in object
    mj = minloc(o(Ofollow)%sp(1,:))
    ! get unit normal
    normal = o(Ofollow)%vertnorm(:,mj(1))
    ! vector rejection from normal unit vector (tanjential velocity)
    cam%velocity = cam%velocity - normal*dot_product(cam%velocity,normal)
   endif
  endif

 endif ! end keystroke !}}}

 !stuff that happens independent of key strokes, like environmental dynamics
 if (rotate.or.impulseControl) then
  ! rotate the global rotor by qr
  call quatmult( cam%spin, scr%globalrotor, qtmp2 )
  scr%globalrotor = qtmp2
  ! inverse of the globalrotor is a rotor that rotates back to global coordinates
  call quatinv( scr%globalrotor, cam%orient )
  !keep track of what direction the global axes point
  !call quatrotate( (/0.0,1.0,0.0,0.0/), scr%globalrotor, qtmp2)
  !axisX = qtmp2(2:4)
  call quat3rotate( (/1.0,0.0,0.0/), scr%globalrotor, axisX)
  !call quatrotate( (/0.0,0.0,1.0,0.0/), scr%globalrotor, qtmp2)
  !axisY = qtmp2(2:4)
  call quat3rotate( (/0.0,1.0,0.0/), scr%globalrotor, axisY)
  !call quatrotate( (/0.0,0.0,0.0,1.0/), scr%globalrotor, qtmp2)
  !axisZ = qtmp2(2:4)
  call quat3rotate( (/0.0,0.0,1.0/), scr%globalrotor, axisZ)
  !keep track of what direction the global axes point
  !call quatrotate( (/0.0,1.0,0.0,0.0/), cam%orient, qtmp2)
  !caX = qtmp2(2:4)
  call quat3rotate( (/1.0,0.0,0.0/), cam%orient, caX)
  call cart2sph( caX, caXs ) ! to get theta and phi values for the heading
  !call quatrotate( (/0.0,0.0,1.0,0.0/), cam%orient, qtmp2)
  !caY = qtmp2(2:4)
  call quat3rotate( (/0.0,1.0,0.0/), cam%orient, caY)
  !call quatrotate( (/0.0,0.0,0.0,1.0/), cam%orient, qtmp2)
  !caZ = qtmp2(2:4)
  call quat3rotate( (/0.0,0.0,1.0/), cam%orient, caZ)
 endif

 if (translate.or.impulseControl) then
   cam%po = cam%po +cam%velocity*scr%dt
   if (Ofollow.gt.0) then !make distance correction to follow surface
    ! closest vertex in object
    mj = minloc(o(Ofollow)%sp(1,:))
    ! get unit normal
    normal = o(Ofollow)%vertnorm(:,mj(1))
    ! project relative displacement onto the normal  
    ! vector rejection from normal unit vector (tanjential velocity)
    cam%velocity = cam%velocity - normal*dot_product(cam%velocity,normal)
    ! move down so that the cam is a radius from surface
    cam%po = cam%po +normal*(dot_product(o(Ofollow)%rp(:,mj(1)),normal) -cam%radius)
     
   endif
 endif

 if (run)     call objectDynamics
 if (contact) call objectContact( contactType, cameraContact )

 if (run.or.(rotate.or.translate)) then
  call fb%clear !clear the fb%pxbuff and fb%zbuff
  call drawscreen
 endif

 ! record or take screen shot
 if (record.or.scrot) then
  write(fc,'(i5.5)') t
  call fb%save(trim(scr%dumpName)//fc//".ppm",2) !dump the screen to PPM6 file
  if (scrot) scrot=.false. ! toggle off
 endif

 !!!! EXIT CONDITIONS !!!!!!!!!!!!!
 ! exit if collided with object
 if (cam%collidedWithObjID.ne.0) then !take damage?
   write(6,*) "YOU DIED!  Hit by Asteroid."
   exit !exiting the loop will quit program
 endif

 ! exit if timeout set
 if (scr%timeout.gt.0) then
  if (t.ge.scr%timeout) then
    write(6,*) "TIMED OUT"
    exit
  endif
 endif

 ! exit if there are no objects
 if (Nobjects.eq.0) then
   write(6,*) "YOU WON!  All Asteroids Destroyed!"
   exit
 endif

enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


call fb%close
call kill_screen( scr%tmpdir )
call CPU_TIME(end_time)
 write(6,*) "execution time: ",end_time-start_time

end program view !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine drawscreen !{{{
use fbMod
use screenstuff
use quaternion !, only : cart2sph, rotoraxis, quatrotate
implicit none
integer :: i, x(3), y(3), sr(2,2), myx, wi
character(4) :: px ! B, G, R, A?
real(4) :: ts(3), tt(3), pi, q(4), q2(4), pr(2,2), XY(2,Nobjects), r(3), v(3)
real(4) ::  ri, dx, r2, v2, rv
character(80) :: ctmp
type(PixBuffType) :: pb !proximity plot

!pi = 4.0*ATAN2(1.0,1.0) !3.1415926535
pi = ACOS(-1.0)
! draw box around the screen
px = char(255)//char(255)//char(255)//char(0)
call fb%rec(scr%sr(1,1),scr%sr(2,1), scr%sr(1,2),scr%sr(2,2), px)

! draw camera axis instrument
call fb%circle(scr%sr(1,2)+100, scr%sr(2,2)-100, 50, px)
px = char(0)//char(0)//char(255)//char(0)
x(1) = scr%sr(1,2) +100 +nint(axisX(2)*50.0)
y(1) = scr%sr(2,2) -100 -nint(axisX(1)*50.0)
call fb%line(scr%sr(1,2)+100, scr%sr(2,2)-100, x(1), y(1), px )
px = char(0)//char(255)//char(0)//char(0)
x(1) = scr%sr(1,2) +100 +nint(axisY(2)*50.0)
y(1) = scr%sr(2,2) -100 -nint(axisY(1)*50.0)
call fb%line(scr%sr(1,2)+100, scr%sr(2,2)-100, x(1), y(1), px )
px = char(255)//char(0)//char(0)//char(0)
x(1) = scr%sr(1,2) +100 +nint(axisZ(2)*50.0)
y(1) = scr%sr(2,2) -100 -nint(axisZ(1)*50.0)
call fb%line(scr%sr(1,2)+100, scr%sr(2,2)-100, x(1), y(1), px )

! Number of objects 
px = char(255)//char(255)//char(255)//char(0)
call fb%putString("Objects Remaining",scr%sr(1,2)+10, scr%sr(2,1)+10, 1, px)
call fb%putNumber(mod(Nobjects,1000)/100, scr%sr(1,2)+125, scr%sr(2,1)+7, 3, px )
call fb%putNumber(mod(Nobjects,100)/10, scr%sr(1,2)+131, scr%sr(2,1)+7, 3, px )
call fb%putNumber(mod(Nobjects,10), scr%sr(1,2)+137, scr%sr(2,1)+7, 3, px )

! Draw Heading, Position, Velocity and Rotation
!{{{
! HEADING
px = char(0)//char(255)//char(255)//char(0)  !yellow
write(ctmp,'(3(f5.2,x))') caX
call fb%putString("heading", scr%sr(1,2)+10, scr%sr(2,1)+53, 1, px)
call fb%putString("(xyz): "//trim(ctmp), &
       &  scr%sr(1,2)+15, scr%sr(2,1)+60, 1, px)
call cart2sph( caX, ts )
write(ctmp,'(2(f5.2,x))') ts(2:3)/pi
call fb%putString("( "//char(16)//char(17)//"):       "//trim(ctmp)//" "//char(18), &
       &  scr%sr(1,2)+15, scr%sr(2,1)+67, 1, px)
! POSITION
px = char(0)//char(255)//char(100)//char(0)  !greenish
write(ctmp,'(3(f8.2,x))') cam%po
call fb%putString("position", scr%sr(1,2)+10, scr%sr(2,1)+80, 1, px)
call fb%putString("(xyz): "//trim(ctmp), &
       &  scr%sr(1,2)+15, scr%sr(2,1)+87, 1, px)
call cart2sph( cam%po, ts )
write(ctmp,'(3(f8.2,x))') ts(1),ts(2:3)/pi
call fb%putString("(r"//char(16)//char(17)//"): "//trim(ctmp)//" "//char(18), &
       &  scr%sr(1,2)+15, scr%sr(2,1)+94, 1, px)
! VELOCITY
px = char(0)//char(200)//char(0)//char(0)   !green
write(ctmp,'(3(f8.2,x))') cam%velocity
call fb%putString("velocity", scr%sr(1,2)+10, scr%sr(2,1)+107, 1, px)
call fb%putString("(xyz): "//trim(ctmp), &
       &  scr%sr(1,2)+15, scr%sr(2,1)+114, 1, px)
call cart2sph( cam%velocity, ts )
write(ctmp,'(3(f8.2,x))') ts(1),ts(2:3)/pi
call fb%putString("(r"//char(16)//char(17)//"): "//trim(ctmp)//" "//char(18), &
       &  scr%sr(1,2)+15, scr%sr(2,1)+121, 1, px)
! ROTATION
px = char(100)//char(200)//char(0)//char(0)   !cyanish
call rotoraxis( cam%spin, q )
write(ctmp,'(4(f8.2,x))') q
call fb%putString("rotation", scr%sr(1,2)+10, scr%sr(2,1)+134, 1, px)
call fb%putString("(axyz):"//trim(ctmp), &
       &  scr%sr(1,2)+15, scr%sr(2,1)+141, 1, px)
call cart2sph( q(2:4), ts )
write(ctmp,'((f8.2,x),9x,2(f8.2,x))') q(1)/pi,ts(2:3)/pi
call fb%putString("(a "//char(16)//char(17)//"):"//trim(ctmp)//" "//char(18), &
       &  scr%sr(1,2)+15, scr%sr(2,1)+148, 1, px)
!}}}

! Draw the object Proximity map
!{{{
wi = 0; ri = 100.0
do i = 1, Nobjects
 XY(:,i) = 0.0
 r(1) = o(i)%po(1) - cam%po(1)
 r(2) = o(i)%po(2) - cam%po(2)
 r(3) = o(i)%po(3) - cam%po(3)
 r2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3) 

 v(1) = o(i)%velocity(1) - cam%velocity(1)
 v(2) = o(i)%velocity(2) - cam%velocity(2)
 v(3) = o(i)%velocity(3) - cam%velocity(3)
 v2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3) 
 rv = r(1)*v(1) + r(2)*v(2) + r(3)*v(3) 

 ! time to minimum distance
 if (v2.lt.1.e-6) then; XY(1:2,i) = 0.0
 else;                  XY(2,i) = -rv/v2
 ! (4ac-b*b)/(4a) !min value of quadratic
   XY(1,i) = (v2*r2 -rv*rv)/v2
 endif
 ! minimum distance; if this argument is negative it is a collision
 if (XY(1,i).le.0.0) then; XY(1,i) = max(maxval(scr%UniParms), 50.0)
 else;                     XY(1,i) = sqrt(XY(1,i)) ; endif
 if (XY(1,i).lt.ri) then
  ri = XY(1,i); wi = i !save worst condition
 endif
enddo
pr(1,:) = (/ 0.0, 0.5*max(maxval(scr%UniParms), 50.0) /) !distance range to plot
pr(2,:) = (/ 0.0, 2.0 /) !time range to plot
sr(1,:) = (/ 7, 99 /) !local pixel ranges
sr(2,:) = (/ 1, 93 /)

px = char(155)//char(155)//char(155)//char(0) !grey
call pb%init(100,100,100, .false.) !100x100 pixbuff
call pb%plot( XY, pr, sr, px )
! draw axis lines
call pb%line(sr(1,1),sr(2,1), sr(1,1),sr(2,2), px)
call pb%line(sr(1,1),sr(2,2), sr(1,2),sr(2,2), px)
! local pixel at x=R1+R2 is int(r/dx+pr(1,1)) + sr(1,1) +1
dx = (pr(1,2)-pr(1,1))/float(sr(1,2)-sr(1,1)-1)
myx = int((o(wi)%radius+cam%radius)/dx+pr(1,1)) +sr(1,1) +1
px = char(0)//char(0)//char(100)//char(0) !dark red
call pb%line(myx,sr(2,1), myx,sr(2,2), px) ! line of collision
px = char(155)//char(155)//char(155)//char(0) !grey
call pb%putString("t",1,40, 1, px)   !time label
call pb%putString("i",1,46, 1, px)   !time label
call pb%putString("m",1,52, 1, px)   !time label
call pb%putString("e",1,58, 1, px)   !time label
call pb%putString("min dist",20,95, 1, px)

px = char(255)//char(255)//char(255)//char(0) !white
call fb%putPixBuff( scr%sr(1,2)+50, scr%sr(2,1)+180, pb)
!}}}

call drawObjects

! draw a HUD or something on top of what is rendered
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  {{{
! center of view screen
x(1) = (scr%sr(1,2)-scr%sr(1,1))/2 +scr%sr(1,1)
y(1) = (scr%sr(2,2)-scr%sr(2,1))/2 +scr%sr(2,1)
px = char(0)//char(200)//char(0)//char(0)   !green
! draw a small cross-hairs
call fb%line( x(1), y(1)-5, x(1), y(1)-2, px )
call fb%line( x(1)-5, y(1), x(1)-2, y(1), px )
call fb%line( x(1), y(1)+2, x(1), y(1)+5, px )
call fb%line( x(1)+2, y(1), x(1)+5, y(1), px )

q(1) = 0.0; q(2:4) = cam%velocity
call quatrotate( q, scr%globalrotor, q2)
! get spherical velocity vector (r, theta, phi)
call cart2sph( q2(2:4), ts )
call cart2sph( -q2(2:4), tt )
if (abs(ts(1)).gt.1.e-12) then !nonzero velocity
 if (ts(3).gt.-0.5*scr%FOVx.and.ts(3).lt.0.5*scr%FOVx) then
  if (ts(2).gt.0.5*(pi-scr%FOVy).and.ts(2).lt.0.5*(scr%FOVy+pi)) then
   ! where on the screen is this vector located?
   x(1) = int((ts(3)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
   y(1) = int(scr%sr_y*((ts(2)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
   ! draw a small diamond
   call fb%line( x(1)-6, y(1), x(1), y(1)+6, px ) !upper V
   call fb%line( x(1)+6, y(1), x(1), y(1)+6, px )
   px = char(0)//char(100)//char(0)//char(0)   !dark green
   call fb%line( x(1)-6, y(1), x(1), y(1)-6, px ) !lower V
   call fb%line( x(1)+6, y(1), x(1), y(1)-6, px )
  endif
 elseif (tt(3).gt.-0.5*scr%FOVx.and.tt(3).lt.0.5*scr%FOVx) then
  if (tt(2).gt.0.5*(pi-scr%FOVy).and.tt(2).lt.0.5*(scr%FOVy+pi)) then
   ! where on the screen is this vector located?
   x(1) = int((tt(3)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
   y(1) = int(scr%sr_y*((tt(2)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
   ! draw a small diamond
   px = char(0)//char(0)//char(220)//char(0)   !dark red
   call fb%line( x(1)-6, y(1), x(1), y(1)+6, px ) !upper V
   call fb%line( x(1)+6, y(1), x(1), y(1)+6, px )
   px = char(0)//char(0)//char(120)//char(0)   !dark red
   call fb%line( x(1)-6, y(1), x(1), y(1)-6, px ) !lower V
   call fb%line( x(1)+6, y(1), x(1), y(1)-6, px )
  endif
 endif
endif

! get rotation axis
call rotoraxis( cam%spin, q )
if (abs(q(1)).gt.1.e-12) then ! there is finite rotation
 call cart2sph( q(2:4), ts ) ! axis vector in spherical
 call cart2sph( -q(2:4), tt ) ! axis vector in spherical
 if (ts(3).gt.-0.5*scr%FOVx.and.ts(3).lt.0.5*scr%FOVx) then
  if (ts(2).gt.0.5*(pi-scr%FOVy).and.ts(2).lt.0.5*(scr%FOVy+pi)) then
   ! where on the screen is this vector located?
   x(1) = nint((ts(3)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
   y(1) = nint(scr%sr_y*((ts(2)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
   px = char(0)//char(150)//char(0)//char(0)   !medium green
   call fb%circle( x(1), y(1), 7, px )
   call fb%line( x(1)+3, y(1)+3,  x(1)+5, y(1)+5,  px )
   call fb%line( x(1)-3, y(1)-3,  x(1)-5, y(1)-5,  px )
   call fb%line( x(1)+3, y(1)-3,  x(1)+5, y(1)-5,  px )
   call fb%line( x(1)-3, y(1)+3,  x(1)-5, y(1)+5,  px )
  endif
 elseif (tt(3).gt.-0.5*scr%FOVx.and.tt(3).lt.0.5*scr%FOVx) then
  if (tt(2).gt.0.5*(pi-scr%FOVy).and.tt(2).lt.0.5*(scr%FOVy+pi)) then
   ! where on the screen is this vector located?
   x(1) = nint((tt(3)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
   y(1) = nint(scr%sr_y*((tt(2)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
   px = char(0)//char(0)//char(170)//char(0)   !medium red (little brighter to see)
   call fb%circle( x(1), y(1), 7, px )
   call fb%line( x(1)+5, y(1)+5,  x(1)+7, y(1)+7,  px )
   call fb%line( x(1)-5, y(1)-5,  x(1)-7, y(1)-7,  px )
   call fb%line( x(1)+5, y(1)-5,  x(1)+7, y(1)-7,  px )
   call fb%line( x(1)-5, y(1)+5,  x(1)-7, y(1)+7,  px )
  endif
 endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!}}}
call fb%display !write fb%pxbuff to the framebuffer device

end subroutine drawscreen !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine breakObject( ob ) !{{{
use quaternion
use renderMod
implicit none
integer, intent(in) :: ob
integer :: t, t2, t1
real(4) :: cs
real(4), dimension(3) :: v1, v2, v3
logical :: found

found = .false.
! let's check this object size and how many times we've colored it red
do t=1, o(ob)%np;                  if (o(ob)%color(1,t).eq.255) then;  
  if (o(ob)%color(2,t).eq.0) then;  if (o(ob)%color(3,t).eq.0) then
    v2 = o(ob)%po - o(ob)%ip(:,t); found = .true.
  endif; endif
endif; enddo 
t1 = scr%centerNode
! make point under AIM red
o(scr%centerObj)%color(1:3,t1) = (/255, 0, 0/) !red

if (found) then
  ! then break this object in half!!!
   v1 = o(ob)%po - o(ob)%ip(:,t1) 
  ! plane made by colored points and center
   call cross_product( v1, v2, v3 )
   cs = sqrt(dot_product(v3,v3))
   if (cs.lt.1.e-12) then !very unlikely
     o(scr%centerObj)%color(1,t1) = 254 !red
     Return !this plane is undefined, the two red spots are aligned
   endif
   v3 = v3/cs !unit vector
  
   call o(ob)%scale(0.793698) ! make first of 2 fragments with half the mass
   do t = 1, o(ob)%np
     t2 = int(255.0*o(ob)%radius/2.5)
     o(ob)%color(1:3,t) = (/t2,t2,255/)
   enddo
  ! move the fragment over by a bit
   o(ob)%po = o(ob)%po -v3*o(ob)%radius*0.5
   o(ob)%velocity = o(ob)%velocity -v3*20.0
  ! copy this for the second fragment
   t2 = Nobjects + 1
   o(t2) = o(ob)  !copy
   o(t2)%po = o(t2)%po +v3*o(ob)%radius !offset
   o(t2)%velocity = o(t2)%velocity +v3*20.0
   Nobjects = t2  !update number of objects
else ! too small, so just destroy it
   if (o(ob)%radius.lt.1.5) o(ob)%np = 0 ! this will  not render it.
endif

end subroutine breakObject !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine getOptions !{{{
use ArgsMod
use options
implicit none
!real(4), dimension(1) :: realv
!integer, dimension(1) :: intv
character(len=80) :: cerr
character(len=80), dimension(1) :: charv


!set default commandline options
FILobj = "NULL"
FILscreen = "screen.rc"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! get any options from the command line 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! get the first argument to check for 'help' switch
call getarg( 1, cerr )
if (trim(cerr).eq."-h".or. trim(cerr).eq."--help") then
   call help
endif

! check for the options to the program
call getOpt( "-i", 1, charv, cerr )
if (cerr(1:1).eq.'0') FILobj = charv(1)
write(0,*) cerr
call getOpt( "-s", 1, charv, cerr )
if (cerr(1:1).eq.'0') FILscreen = charv(1)
write(0,*) cerr

end subroutine getOptions !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine help !{{{
use renderMod
write(0,*) "view [OPTIONS]"
write(0,*) "OPTIONS"
write(0,*) "  -i FILE   object file to view"
write(0,*) "  -s FILE   screen file with screen configuration data"
write(0,*) "  -h or --help  this help message"
write(0,*) 
write(0,*) "KEYS"
write(0,*) "  F10 or Q  quits the program"
write(0,*) "  UP        pitch down"
write(0,*) "  DOWN      pitch up"
write(0,*) "  LEFT      roll left"
write(0,*) "  RIGHT     roll right"
write(0,*) "  q         turn left"
write(0,*) "  e         turn right"
write(0,*) "  a         translate left"
write(0,*) "  d         translate right"
write(0,*) "  w         translate forward"
write(0,*) "  s         translate backward"
write(0,*) "  W         translate up"
write(0,*) "  S         translate down"
write(0,*) "  p         save entire framebuffer to ppm file defined in screen.rc"
write(0,*) "  r         toggle record mode, series of saves."
call commandHelp
STOP
end subroutine help !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
