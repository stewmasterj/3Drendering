! vim:fdm=marker
! requires lineParse.f90
! requires quaternion.f90
! DATE: the same as for the 'view' programs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
module renderMod

type object
 character(80) :: name
 integer :: np !number of points
 integer :: nt !number of triangles
 integer :: nvn !number of vertex normals
 character(5) :: mode !(point|wire|solid)
 real(4) :: mass, radius, stiffness, poisson, friction
 logical :: transparency ! flag for certain point-wise plotting
 real(4), allocatable, dimension(:,:) :: point,ip,rp,sp
 real(4), allocatable, dimension(:,:) :: vertnorm, ivn, rvn
 integer, allocatable, dimension(:,:) :: color, connect
 integer, allocatable, dimension(:)   :: smooth !average normals or not
 integer, allocatable, dimension(:,:) :: tri
 integer, allocatable, dimension(:)   :: sortID !sorted list of IDs
 ! changable values
 real(4), dimension(3) :: po !position offset wrt global coordinates
 real(4), dimension(3) :: velocity !changes the offset position
 real(4), dimension(4) :: orient !orientation quaternion
 real(4), dimension(4) :: spin !local rotation rate (changes orientation)
 integer :: collidedWithObjID !what object did this object collide with?
contains
 procedure :: ene => objectEnergy
 procedure :: momentum => objectMomentum
 procedure :: scale => objectScale
 procedure :: reallocate => objectReallocate
 procedure :: addPoint => objectAddPoint
 procedure :: findVertexNormal 
 procedure :: findTriangleNormal 
 procedure :: closestTriangle
 procedure :: mkconnections
 procedure :: sortNodes !only populates 'sortID'
end type object

type screen
 character(80) :: fbpath, tmpdir, dumpName
 integer :: width, height, line, timeout, sr_x, sr_y
 real(4) :: FOVx, FOVy
 integer, dimension(2,2) :: sr
 real(4) :: dt, dtheta, ds
! real(4), dimension(3) :: camera_p, camera_vel
 real(4), dimension(4) :: globalrotor !, camerarotor, camera_spin
 logical :: interactive
 integer :: centerObj, centerTri, centerNode, seed(1)
 ! virtual universe data
 character(10) :: UniGeom !sphere or box or null
 logical :: periodic   !fold positions or rebound
 real(4), dimension(3) :: UniParms ! parameters (r for sphere, L(x,y,z) for  box)
 character(4) :: bgpx
end type


! this module also has these variables
type(object), dimension(:), allocatable :: o
type(object) :: cam !camera object, not part of all other objects
type(screen) :: scr
integer :: ObjectBufferSize, Nobjects

!!!!!!!!!!!!!!!!!
! must first call loadScreenData 
!    this allocates necessary arrays and initialized variables
!!!!!!!!!!!!!!!!!!!

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine loadObject( FD, fl, obj ) !{{{
! read an object from file descriptor FD and file name fl
use lineParse
use quaternion
implicit none
character(*) :: fl
character(80) :: line, word, tmp
integer :: FD, err, i, ln, id, sm, c(3)
real(4) :: p(3), q1(4), q2(4)
type(object), intent(inout) :: obj

obj%name = "Object_Name"
obj%mode = "point"
obj%po(:) = 0.0
obj%velocity(:) = 0.0
obj%orient(:) = 0.0; obj%orient(2:4) = 1.0
obj%spin(:) = 0.0; obj%spin(2:4) = 1.0
! for contact
obj%mass = 1.0
obj%radius = 1.0
! for hertzian contact
obj%stiffness = 1.0
obj%poisson = 0.3
obj%friction = 1.0

open(FD,file=trim(fl),iostat=err)
if (err.ne.0) then
  STOP "some problem opening object file:"//trim(fl)
endif
do 
 line = getLine( FD, ln, err )
 if (err.ne.0) exit
 word = s_get_word(1, line)
 select case (trim(word))
  case ("name"); obj%name = s_get_word(2, line)
  case ("mode"); obj%mode = trim(s_get_word(2, line))
  case ("mass"); call s_get_val(2, line, obj%mass )
  case ("radius"); call s_get_val(2, line, obj%radius )
  case ("stiffness"); call s_get_val(2, line, obj%stiffness )
  case ("poisson"); call s_get_val(2, line, obj%poisson )
  case ("friction"); call s_get_val(2, line, obj%friction )
  case ("offset"); 
   call s_get_val(2, line, obj%po(1) )
   call s_get_val(3, line, obj%po(2) )
   call s_get_val(4, line, obj%po(3) )
  case ("velocity"); 
   call s_get_val(2, line, obj%velocity(1) )
   call s_get_val(3, line, obj%velocity(2) )
   call s_get_val(4, line, obj%velocity(3) )
  case ("orient"); 
   call s_get_val(2, line, obj%orient(1) )
   call s_get_val(3, line, obj%orient(2) )
   call s_get_val(4, line, obj%orient(3) )
   call s_get_val(5, line, obj%orient(4) )
   q1(1) = 0.0; q1(2:4) = obj%orient(2:4) !axis portion
   call quatrotor( q1, obj%orient(1), q2 )
   obj%orient = q2
  case ("spin"); 
   call s_get_val(2, line, obj%spin(1) )
   call s_get_val(3, line, obj%spin(2) )
   call s_get_val(4, line, obj%spin(3) )
   call s_get_val(5, line, obj%spin(4) )
   q1(1) = 0.0; q1(2:4) = obj%spin(2:4) !axis portion
   call quatrotor( q1, obj%spin(1), q2 )
   obj%spin = q2
  case ("points"); 
   call s_get_val(2, line, obj%np )
   allocate( obj%point(3,obj%np), obj%color(3,obj%np), obj%smooth(obj%np) )
   allocate( obj%ip(3,obj%np), obj%rp(3,obj%np), obj%sp(3,obj%np) )
   do i = 1, obj%np
    read(FD,*,iostat=err) id, p(1:3), c(1:3), sm
    if (IS_IOSTAT_EOR(err)) then
     STOP "ERROR: premature End of Record reading point data"
    elseif (IS_IOSTAT_END(err)) then
     STOP "ERROR: premature End of File reading point data"
    endif
    obj%point(1:3,id) = p(1:3)
    obj%color(1:3,id) = c(1:3)
    obj%smooth(id)    = sm !zero is sharp (no averaging), one will average normals
!    obj%ip(:,i) = obj%point(:,i) +obj%po(:) !global coordinates
   enddo
  case ("triangles"); tmp = s_get_word(2, line); read(tmp,*) obj%nt
   allocate( obj%tri(3,obj%nt) )
   do i = 1, obj%nt
    read(FD,*,iostat=err) c(1:3)
    if (IS_IOSTAT_EOR(err)) then
     STOP "ERROR: premature End of Record reading triangle data"
    elseif (IS_IOSTAT_END(err)) then
     STOP "ERROR: premature End of File reading triangle data"
    endif
    obj%tri(1:3,i) = c(1:3)
   enddo
  case ("VertexNormals"); 
   call s_get_val(2, line, obj%nvn )
   allocate( obj%vertnorm(3,obj%np), obj%ivn(3,obj%np), obj%rvn(3,obj%np) ) !allocate for all points, usually all or none are smooth.
   do i = 1, obj%nvn
    read(FD,*,iostat=err) id, p(1:3)
    if (IS_IOSTAT_EOR(err)) then
     STOP "ERROR: premature End of Record reading VertexNormals data"
    elseif (IS_IOSTAT_END(err)) then
     STOP "ERROR: premature End of File reading VertexNormals data"
    endif
    if (obj%smooth(id).eq.1) then !zero is sharp (no averaging), one will average normals
     obj%vertnorm(1:3,id) = p(1:3)
    else
     write(0,*) "WARNING: Vertex Normal line:",i," ID:",id," not flagged for smoothing"
    endif
   enddo
  case default !ignores lines that don't match anything
 end select
enddo
close(FD)

! if object has a local rotation WRT global, apply it here.
do i = 1, obj%np
  obj%ip(:,i) = obj%point(:,i) +obj%po(:) !global coordinates
enddo


end subroutine loadObject !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine xyz2Object( FD, fl, obj ) !{{{
! read an object from file descriptor FD and file name fl
! col : column array from file corresponding to
! col(1:3) x,y,z columns
! col(4)   point radius
! col(5)   pseudo transparency (current color)*(t)+(point color)*(1-t)
! col(6:8) either RGB or HSV columns
! rng(2,8) min and max limits for each column value in col(:)
!          if min=max & col(i)=0 then do not read a column, but use the rng value
!          if min=max=0.0 & col(i)>0 then auto calculate range to fit
! set mode either sphere or circle or fcircle (filled circle)
!use lineParse
use quaternion
use inputShareMod ! for col, rng, h
use fbMod !only HSV2RGB
implicit none
character(*) :: fl
!integer, intent(in) :: col(8)
!real(4), intent(in) :: rng(2,8)
type(object), intent(inout) :: obj
!logical, intent(in) :: h ! is the input color columns for HSV?
! internal
character(80) :: line
integer :: FD, err, i, j, ln, c(4), r(3), skip
real(4) :: p(4), mr(2,8), q1(4), q2(4)
real(4), allocatable, dimension(:) :: dum

! allocate dummy for max needed columns
allocate( dum( maxval(col) ) )

obj%name = trim(fl)
obj%mode = trim(md)
obj%po(:) = 0.0
obj%velocity(:) = 0.0
!obj%orient(:) = 0.0; obj%orient(2:4) = 1.0
   q1(1) = 0.0; q1(2:4) = 1.0 !obj%orient(2:4) !axis portion
   call quatrotor( q1, 0.0, q2 )
   obj%orient = q2
!obj%spin(:) = 0.0; obj%spin(2:4) = 1.0
   q1(1) = 0.0; q1(2:4) = 1.0 !obj%orient(2:4) !axis portion
   call quatrotor( q1, 0.0, q2 )
   obj%spin = q2
! for contact
obj%mass = 1.0
obj%radius = 1.0
! for hertzian contact
obj%stiffness = 1.0
obj%poisson = 0.3
obj%friction = 1.0
!
obj%nt = 0 ! no triangle data

open(FD,file=trim(fl),iostat=err)
if (err.ne.0) then
  STOP "some problem opening xyz file:"//trim(fl)
endif

! first word on first line is number of points
!read(FD,*,iostat=err) obj%np 
!if (err.ne.0) exit
!read(FD,*,iostat=err) ! empty line or comments

! scan the file to get number of points
i = 0; skip = 0
do
 read(FD,'(A)',iostat=err) line
 if (IS_IOSTAT_END(err)) exit
 if (line(1:1).eq."#") then; skip = skip + 1; cycle; endif
 i = i + 1
enddo
obj%np = i
write(6,*) "Read xyz file. comments:", skip, " lines:", i
write(6,*) "rereading..."

rewind(FD)
do i = 1, skip; read(FD,*); enddo ! skip the commented header

allocate( obj%point(3,obj%np), obj%color(3,obj%np), obj%smooth(obj%np), stat=err )
if (err.ne.0) STOP "ERROR: allocating points color smooth from data file"
allocate( obj%ip(3,obj%np), obj%rp(3,obj%np), obj%sp(3,obj%np), stat=err )
if (err.ne.0) STOP "ERROR: allocating ip rp sp from data file"

obj%nvn = obj%np ! just make nvn=np to have a float number
allocate( obj%vertnorm(3,obj%np), stat=err )
if (err.ne.0) STOP "ERROR: allocating vertnorm from data file"

mr = rng ! copy this here since 'mr' gets updated

ln = skip
do i = 1, obj%np
 ln = ln+1
 read(FD,*,iostat=err) dum(:)
 if (err.ne.0) then
    if (IS_IOSTAT_EOR(err)) then
     write(0,*) "ERROR: premature End of Record reading point data in xyz file, line:", ln
    elseif (IS_IOSTAT_END(err)) then
     write(0,*) "ERROR: premature End of File reading point data in xyz file, line:", ln
    else
     write(0,*) "ERROR: reading point data in xyz file, line: ", ln
    endif
    STOP
 endif
 obj%point(1,i) = dum(col(1))
 obj%point(2,i) = dum(col(2))
 obj%point(3,i) = dum(col(3))
 !obj%smooth(i)  = 1 !must be flagged otherwise vertex normals won't save in obj file
 obj%smooth(i)  = 0 !must be zero so that dynamics don't rotate the normals
 if (col(4).eq.0) then
  obj%vertnorm(1,i) = rng(1,4) 
 else
  obj%vertnorm(1,i) = abs(dum(col(4))) ! point radius
 endif
 if (col(5).eq.0) then
  obj%vertnorm(2,i) = rng(1,5) ! point transparency
 else
  obj%vertnorm(2,i) = abs(dum(col(5))) ! point transparency
 endif
 obj%vertnorm(3,i) = 0.0

 ! calculate mins and maxes from the data
 do j = 1, 8
  if (col(j).le.0) cycle 
  if (i.eq.1) mr(:,j) = dum(col(j)) ! initialize
  if (dum(col(j)).lt.mr(1,j)) mr(1,j) = dum(col(j))
  if (dum(col(j)).gt.mr(2,j)) mr(2,j) = dum(col(j))
 enddo
enddo


! set global ranges if not specified in 'rng'
do j = 1, 8
  if (col(j).le.0) cycle ! it is set to a value
  if (abs(rng(1,j)-rng(2,j)).lt.1.e-12) then ! it needs to be set
    rng(:,j) = mr(:,j)
  endif
enddo

! if reading transparency values from file column or if it's set to nonzero.
if (col(5).ne.0 .or. abs(rng(1,5)-rng(2,5)).gt.1.e-12) then
  obj%transparency = .true.
endif

write(6,'(a)') "dataRanges read from file if not specified:"
write(6,*) "x: ",rng(:,1);     write(6,*) "y: ",rng(:,2)
write(6,*) "z: ",rng(:,3);     write(6,*) "r: ",rng(:,4)
write(6,*) "t: ",rng(:,5);     write(6,*) "R,h: ",rng(:,6)
write(6,*) "G,s: ",rng(:,7);   write(6,*) "B,v: ",rng(:,8)


rewind(FD)
do i = 1, skip; read(FD,*); enddo ! skip the commented header
! if object has a local rotation WRT global, apply it here.
do i = 1, obj%np
  obj%ip(:,i) = obj%point(:,i) +obj%po(:) !global coordinates

 read(FD,*,iostat=err) dum(:)
 ! color stuff
 ! normalize column value to color range
 do j = 1, 4
  if (col(j+4).eq.0.or.abs(rng(1,j+4)-rng(2,j+4)).lt.1.e-12) then
   p(j) = rng(1,j+4)
   c(j) = nint(p(j))
  else
   if (rng(1,4+j).lt.rng(2,4+j)) then
     p(j) = max( min(dum(col(4+j)),rng(2,4+j)) ,rng(1,4+j)) ! fit the value within the range
   else
     p(j) = max( min(dum(col(4+j)),rng(1,4+j)) ,rng(2,4+j)) ! fit the value within the range
   endif
   p(j) = (p(j)-rng(1,4+j))/(rng(2,4+j)-rng(1,4+j))  ! fraction within the range
   c(j) = nint(p(j)*255.0)
  endif
 enddo ! j
 if (h) then ! provided as HSV values, need them in RGB
  call HSV2RGB( c(2),c(3),c(4), r(1),r(2),r(3) )
  obj%color(1:3,i) = r(1:3)
 else ! they are already in RGB
  obj%color(1:3,i) = c(2:4)
 endif
 obj%vertnorm(2,i) = real(c(1))

enddo

close(FD)

end subroutine xyz2Object !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine saveObject( FD, fl, i ) !{{{
use quaternion
implicit none
integer, intent(in) :: FD, i
character(*) :: fl
integer :: err, j
real(4), dimension(4) :: q

open(FD,file=trim(fl),iostat=err)
if (err.ne.0) then
  STOP "some problem opening output file:"//trim(fl)
endif

write(FD,'(A)') "name "//o(i)%name
write(FD,'(A)') "mode "//o(i)%mode
write(FD,'(A,E12.5)') "mass ",o(i)%mass
write(FD,'(A,E12.5)') "radius ",o(i)%radius
write(FD,'(A,E12.5)') "stiffness ",o(i)%stiffness
write(FD,'(A,E12.5)') "poisson ",o(i)%poisson
write(FD,'(A,E12.5)') "friction ",o(i)%friction
write(FD,'(A,3E12.5)') "offset ",o(i)%po
write(FD,'(A,3E12.5)') "velocity ",o(i)%velocity
! convert back to angle axis notiation from rotor
call rotoraxis( o(i)%orient, q )
write(FD,'(A,4E12.5)') "orient ",q
call rotoraxis( o(i)%spin, q )
write(FD,'(A,4E12.5)') "spin ",q
! write points
write(FD,'(A)') "points"
do j = 1, o(i)%np
   write(FD,*,iostat=err) j, o(i)%point(1:3,j), o(i)%color(1:3,j), o(i)%smooth(j)
   if (err.ne.0) then
      write(0,*) "ERROR writing points in object file"
      return
   endif
enddo
! write triangles
write(FD,'(A)') "triangles"
do j = 1, o(i)%nt
   write(FD,*,iostat=err) o(i)%tri(1:3,j)
   if (err.ne.0) then
      write(0,*) "ERROR writing triangles in object file"
      return
   endif
enddo
! write Vertex Normals
write(FD,'(A)') "VertexNormals"
do j = 1, o(i)%np
   if (o(i)%smooth(j).eq.1) then !write it
      write(FD,*,iostat=err) j, o(i)%vertnorm(1:3,j)
      if (err.ne.0) then
         write(0,*) "ERROR writing vertex normals in object file"
         return
      endif
   endif
enddo

close(FD)

end subroutine saveObject !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine loadScreenData( FD, fl ) !{{{
! read screen data from file descriptor FD and file name fl
use lineParse
use quaternion
!use starMod, only : starfilename, Nstars
use controlMod
use inputShareMod !for col, rng, h, md
implicit none
character(80) :: fl, line
integer :: FD, err, ln
real(4) :: q1(4), q2(4)

scr%fbpath = "/dev/fb0"
scr%tmpdir = "/tmp"
scr%width = 1440 !these three are just copied into the framebuffer object
scr%height = 900
scr%line = 1472
scr%FOVx = 1.74533 ! field of view (radians)
scr%FOVy = 1.5708
scr%sr(1,:) = (/200, 700/)  ! local screen position for the rendering
scr%sr(2,:) = (/100, 550/)
scr%dt = 0.01   !time step
scr%dtheta = 0.261799  ! angular increment for movement
scr%ds = 1.0    ! spatial increment for movement

scr%interactive = .false.
scr%dumpName = "_"
scr%timeout = 0  !stop program after this many time steps
! universe data
scr%UniGeom = "NULL" !no virtual position boundaries
scr%periodic = .false. !not periodic (would rebound if geom set
scr%UniParms(:) = 0.0

scr%seed(1) = 123456789
scr%globalrotor = (/ 0.0, 1.0, 1.0, 1.0 /) !initial angle-axis notation
scr%bgpx = char(0)//char(0)//char(0)//char(0)

! camera object data
cam%name = "Camera"
cam%po  = (/ -10.0, 0.0, 0.0 /) !default camera position
cam%velocity  = (/ 0.0, 0.0, 0.0 /) ! yoru personal velocity
cam%orient = (/ 0.0, 1.0, 1.0, 1.0 /) !initial angle-axis notation
cam%spin = (/ 0.0, 1.0, 1.0, 1.0 /) !initial angle-axis notation
q1(1) = 0.0; q1(2:4) = cam%orient(2:4) !axis portion
call quatrotor( q1, cam%orient(1), q2 )
cam%orient = q2
q1(1) = 0.0; q1(2:4) = cam%spin(2:4) !axis portion
call quatrotor( q1, cam%spin(1), q2 )
cam%spin = q2
! for camera contact
cam%mass = 1.0
cam%radius = 1.0
! for hertzian camera contact
cam%stiffness = 1.0
cam%poisson = 0.3
cam%friction = 1.0


ObjectBufferSize = 1 ! load only one object by default
Nobjects = 0 ! everytime LoadObject is called, this will increment

impulseControl = .false. !control is spatial rather than rate modulating
Ofollow = 0 ! object to follow, zero for no follow
contact = .false. !contact detection off by default
cameraContact = .false. !objects don't interact with camera by default
!starfilename = "stars.dat"
!Nstars = 300

record = .false. ! dump PPM files every frame?
scrot = .false.  ! single screen shot flag

! stuff that persists in input file only used for plotting
col(:) = 0
rng(:,:) = 0.0
rng(:,4) = 1.0 ! radius
h = .true.
md = "point"

open(FD,file=trim(fl),iostat=err)
if (err.ne.0) then
  STOP "ERROR opening scene file:"//trim(fl)
endif
ln = 0
do 
  line = getLine( FD, ln, err )
!write(6,'(i3,x,A)') ln,trim(line) !echo the read line
  if (err.ne.0) then ; err=0; exit ; endif
  if (len(trim(line)).le.0) cycle
  call interpretLine( line, err )
  if (err.gt.0) then
    write(0,*) "ERROR: interpreting line: ",ln," of file:"//trim(fl)
    exit
  endif
  if (err.lt.0) exit !user defined exit
enddo
close(FD)
call flush(0)
call flush(6)
if (err.ne.0) STOP "done reading input config."

if (.not.allocated(o)) allocate( o(ObjectBufferSize) )

! convert angle axis notation into actual quaternion rotors
!q1(1) = 0.0; q1(2:4) = scr%camerarotor(2:4) !axis portion
!call quatrotor( q1, scr%camerarotor(1), q2 )
!scr%camerarotor = q2

!q1(1) = 0.0; q1(2:4) = scr%camera_spin(2:4) !axis portion
!call quatrotor( q1, scr%camera_spin(1), q2 )
!scr%camera_spin = q2


!scr%sr_x = scr%sr(1,2)-scr%sr(1,1) ! x width for 3d screen
!scr%sr_y = scr%sr(2,2)-scr%sr(2,1) ! y width


end subroutine loadScreenData !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine interpretLine( s, error ) !{{{
! this routine could be used to set variables or perform complex things that a
! single input character can't
use lineParse
use quaternion
use fbMod
!use starMod, only : starfilename, Nstars
use controlMod
use inputShareMod !for col, rng, h, md
implicit none
integer :: wc
integer, intent(out) :: error
character(*) :: s
character(80) :: word, tmp, word2
integer :: i, j, k, c(4), mo, mj
real(4) :: p(3), scl, pi
real(4) :: q1(4), q2(4), val(8)
character :: CR, LF

CR = achar(13)
LF = achar(10)
pi = 4.0*ATAN2(1.0,1.0) !=log(-1.0)/i
error = 0 !positive errors are fatal and kill the program, negative are warnings

wc   = s_word_count( s )
word = s_get_word( 1, s )


select case(trim(word))
  case ("help"); call commandHelp; run=.false. 
  case ("q","exit","quit"); error=-1 ! this will flag as error and kill program
  case ("echo"); write(6,'(A)') s(6:); call flush(6)
  case ("exec"); call execute_command_line( s(6:), WAIT=.true. ) 
  case ("fbpath");      if (wc>1) then; scr%fbpath = s_get_word(2, s)
   else; write(6,*) trim(scr%fbpath); run=.false.; endif
  case ("dumpFile");    if (wc>1) then; scr%dumpName = s_get_word(2, s)
   else; write(6,*) trim(scr%dumpName); run=.false.; endif
!  case ("starFile","starfile");    if (wc>1) then; starfilename = s_get_word(2, s)
!   else; write(6,*) trim(starfilename); run=.false.; endif
  case ("tmpdir");      if (wc>1) then; scr%tmpdir = s_get_word(2, s)
   else; write(6,*) trim(scr%tmpdir); run=.false.; endif
  case ("interactive"); scr%interactive = .true.
  case ("FOV");         if (wc>1) then
    call s_get_val(2, s, scr%FOVx )
    call s_get_val(3, s, scr%FOVy )
   else; write(6,*) scr%FOVx, scr%FOVy; run=.false.; endif
  case ("subscreen");   if (wc>1) then
   call s_get_val(2, s, scr%sr(2,1) )
   call s_get_val(3, s, scr%sr(1,1) )
   call s_get_val(4, s, scr%sr(2,2) )
   call s_get_val(5, s, scr%sr(1,2) )
   scr%sr_x = scr%sr(1,2)-scr%sr(1,1) ! x width for 3d screen
   scr%sr_y = scr%sr(2,2)-scr%sr(2,1) ! y width
   else; write(6,*) scr%sr(2,1), scr%sr(1,1), scr%sr(2,2), scr%sr(1,2); run=.false.; endif
  case ("timestep");    if (wc>1) then; call s_get_val(2, s, scr%dt )
   else; write(6,*) scr%dt; run=.false.; endif
!  case ("Nstars");      if (wc>1) then; call s_get_val(2, s, Nstars )
!   else; write(6,*) Nstars; run=.false.; endif
  case ("dtheta");      if (wc>1) then; call s_get_val(2, s, scr%dtheta )
   else; write(6,*) scr%dtheta; run=.false.; endif
  case ("dspace");      if (wc>1) then; call s_get_val(2, s, scr%ds )
   else; write(6,*) scr%ds; run=.false.; endif
  case ("width");       if (wc>1) then; call s_get_val(2, s, scr%width )
   else; write(6,*) scr%width; run=.false.; endif
  case ("height");      if (wc>1) then; call s_get_val(2, s, scr%height )
   else; write(6,*) scr%height; run=.false.; endif
  case ("timeout");     if (wc>1) then; call s_get_val(2, s, scr%timeout )
   else; write(6,*) scr%timeout; run=.false.; endif
  case ("linelength");  if (wc>1) then; call s_get_val(2, s, scr%line )
   else; write(6,*) scr%line; run=.false.; endif
  case ("NobjectBuff"); if (wc>1) then; call s_get_val(2, s, ObjectBufferSize )
    if (.not.allocated(o)) then; allocate( o(0:ObjectBufferSize) )
     o(0)%np = 0 !no background object
    else; deallocate(o); allocate( o(0:ObjectBufferSize) ); endif
   else; write(6,*) ObjectBufferSize; run=.false.; endif
  case ("contact");     if (wc>1) then; call s_get_val(2, s, contactType )
      if (contactType.gt.0) then; contact=.true.
      else;                       contact=.false.; endif
      if (wc>2) call s_get_val(3, s, cameraContact )
   else; write(6,*) "Contact Type:",contactType; run=.false.; endif
!!!! Load data !!!!
  case ("dataColumns","dataCols")
   if (wc==1) then; write(6,'(a,8i3)') "dataColumns (xyzrthsv):",col
    run=.false.
    return; endif
   do i = 2, wc
    call s_get_vali( i, s, col(i-1) ); enddo
  case ("dataRange","dataRanges")
   if (wc==1) then; write(6,'(a)') "dataRanges:"
    write(6,*) "x: ",rng(:,1);     write(6,*) "y: ",rng(:,2)
    write(6,*) "z: ",rng(:,3);     write(6,*) "r: ",rng(:,4)
    write(6,*) "t: ",rng(:,5);     write(6,*) "R,h: ",rng(:,6)
    write(6,*) "G,s: ",rng(:,7);   write(6,*) "B,v: ",rng(:,8)
    run=.false.
    return; endif
   tmp = s_get_word(2, s)
   if     (trim(tmp).eq."1".or.trim(tmp).eq."x") then
    call s_get_valr4( 3, s, rng(1,1) ); call s_get_valr4( 4, s, rng(2,1) )
   elseif (trim(tmp).eq."2".or.trim(tmp).eq."y") then
    call s_get_valr4( 3, s, rng(1,2) ); call s_get_valr4( 4, s, rng(2,2) )
   elseif (trim(tmp).eq."3".or.trim(tmp).eq."z") then
    call s_get_valr4( 3, s, rng(1,3) ); call s_get_valr4( 4, s, rng(2,3) )
   elseif (trim(tmp).eq."4".or.trim(tmp).eq."r") then
    call s_get_valr4( 3, s, rng(1,4) ); call s_get_valr4( 4, s, rng(2,4) )
   elseif (trim(tmp).eq."5".or.trim(tmp).eq."t") then
    call s_get_valr4( 3, s, rng(1,5) ); call s_get_valr4( 4, s, rng(2,5) )
   elseif (trim(tmp).eq."6".or.trim(tmp).eq."h".or.trim(tmp).eq."R") then
    call s_get_valr4( 3, s, rng(1,6) ); call s_get_valr4( 4, s, rng(2,6) )
   elseif (trim(tmp).eq."7".or.trim(tmp).eq."s".or.trim(tmp).eq."G") then
    call s_get_valr4( 3, s, rng(1,7) ); call s_get_valr4( 4, s, rng(2,7) )
   elseif (trim(tmp).eq."8".or.trim(tmp).eq."v".or.trim(tmp).eq."B") then
    call s_get_valr4( 3, s, rng(1,8) ); call s_get_valr4( 4, s, rng(2,8) )
   endif
  case ("dataColorHSV"); h = .true.
  case ("dataColorRGB"); h = .false.
  case ("dataPoint"); md = trim(s_get_word(2, s))
  case ("LoadObject","LoadData"); 
   if (wc==1) then; write(6,*) "Nobjects:",Nobjects; run=.false.; return; endif
   Nobjects = Nobjects + 1
   if (Nobjects.gt.ObjectBufferSize) then
    write(0,*) "ERROR: not loading this object as there is no more space in object buffer"
    error=1; return; endif
   if (.not.allocated(o)) then
    write(0,*) "ERROR: must allocate object buffer before loading object, see directive: NobjectBuff"
    error=1; return; endif
   tmp = s_get_word(2, s)
   if (trim(word).eq."LoadObject") then
    call loadObject( 12, trim(tmp), o(Nobjects) )
   elseif (trim(word).eq."LoadData") then
    ! this requires col, rng, h, md to be set
    if (maxval(col).le.1) then
     write(0,*) "ERROR: Please set 'dataColumns' to read from the data file"
     error=1; return; endif
    call xyz2Object( 12, trim(tmp), o(Nobjects) )
   endif
  case ("LoadBackground");
   if (.not.allocated(o)) then
    write(0,*) "ERROR: must allocate object buffer before loading object, see directive: NobjectBuff"
    error=1; return
   endif
   tmp = s_get_word(2, s)
   call loadObject( 12, trim(tmp), o(0) )
   ! local points rotation to global coordinates
   ! if it's not spinning, this wouldbe fine to do once
   do i = 1, o(0)%np
    q1(1) = 0.0; q1(2:4) = o(0)%point(1:3,i)
    call quatrotate( q1, o(0)%orient, q2)
    o(0)%ip(:,i) = q2(2:4) +o(0)%po(:)
    ! background objects don't use shading so no need to rotate vertex normals
   enddo
!!!!!!!!!!!!!!!!!!!!!!!!!
  case ("camera_pos"); 
   !if (wc==1) then; write(6,*) scr%camera_p; run=.false.; return; endif
   !call s_get_val(2, s, scr%camera_p(1) )
   !call s_get_val(3, s, scr%camera_p(2) )
   !call s_get_val(4, s, scr%camera_p(3) )
   if (wc==1) then; write(6,*) cam%po; run=.false.; return; endif
   call s_get_val(2, s, cam%po(1) )
   call s_get_val(3, s, cam%po(2) )
   call s_get_val(4, s, cam%po(3) )
  case ("camera_orient"); 
   if (wc==1) then; write(6,*) scr%globalrotor; run=.false.; return; endif
   call s_get_val(2, s, scr%globalrotor(1) )
   call s_get_val(3, s, scr%globalrotor(2) )
   call s_get_val(4, s, scr%globalrotor(3) )
   call s_get_val(5, s, scr%globalrotor(4) )
   q1(1) = 0.0; q1(2:4) = scr%globalrotor(2:4) !axis portion
   call quatrotor( q1, scr%globalrotor(1), q2 )
   scr%globalrotor = q2
  case("dynamics") ; call objectDynamics
  case("draw") ; call drawObjects
  case("bgColor") ; c=0; if (wc>1) call s_get_val(2, s, c(3))
    if (wc>2)  call s_get_val(3, s, c(2))
    if (wc>3)  call s_get_val(4, s, c(1))
    if (wc.eq.1) then; write(6,*) c(3),c(2),c(1), " ", c(4); endif
    if (wc>1) scr%bgpx = char(c(1))//char(c(2))//char(c(3))//char(c(4))
  case("fbinit") ; call fb%fbinit(10,scr%fbpath, scr%width, scr%height, scr%line, .true. )
  case("universe") ; scr%UniGeom = trim(s_get_word(2, s))
   if (trim(scr%UniGeom).eq."sphere") then
    call s_get_val(3, s, scr%UniParms(1))
   elseif (trim(scr%UniGeom).eq."box") then
    call s_get_val(3, s, scr%UniParms(1))
    call s_get_val(4, s, scr%UniParms(2))
    call s_get_val(5, s, scr%UniParms(3))
   endif
  case ("seed");  if (wc>1) then; call s_get_val(2, s, scr%seed(1) )
   else; write(6,*) scr%seed(1); run=.false.; endif
  case ("follow");  if (wc>1) then; call s_get_val(2, s, Ofollow )
    if (Ofollow.gt.0) call o(Ofollow)%mkconnections
   else; write(6,*) Ofollow; run=.false.; endif
  case("periodic"); scr%periodic=.true. !set universe to be periodic
  case("fbclear") ; call fb%clear !clear the fb%pxbuff and fb%zbuff
  case("pause"); run=.false. !stop animation
  case("run"); run=.true. !stop animation
  case("record"); if (record) then; record=.false.; else; record = .true.; endif
  case("picture","scrot"); scrot=.true. !record every frame
    if (.not.scr%interactive.and.wc.ge.2) then
      tmp = s_get_word(2, s)
      call fb%save(trim(tmp),2)
      scrot=.false.
    endif
  case("impulseControl"); impulseControl=.true. 
  case("spatialControl"); impulseControl=.false. 
  case("redraw"); call fb%display !redraw the screen based on current and last rendering
  case("cpo"); !copy object to empty objects
   call s_get_val(2, s, mo)
   call s_get_val(3, s, j)
   k = j
   if (wc.eq.4) call s_get_val(4, s, k) !specify range to copy to
   if (k>Nobjects) Nobjects=k
   if (Nobjects>ObjectBufferSize) then
    write(0,*) CR//"Warning: copied objects exceed object buffer size"
    k = ObjectBufferSize
    Nobjects = ObjectBufferSize
   endif
   do i = j, k
    o(i) = o(mo)  ! simple copy object
   enddo
  case("orand"); !randomize object range properties
   tmp = s_get_word( 2, s )
   call s_get_val( 3, s, i)
   call s_get_val( 4, s, j)
   do k = 5, wc
    call s_get_val( k, s, val(k-4))
   enddo
   call ObjectRandom( tmp, i, j, val )
  case("nearest","closest"); !report which object is closest
   scl = o(1)%sp(1,1) !min point location
   mo = 1 !object with min point location
   mj = 1
   do i = 1, Nobjects; do j = 1, o(i)%np
     if (o(i)%sp(1,j).lt.scl) then; mo = i; mj = j; scl = o(i)%sp(1,j); endif
   enddo; enddo
   write(6,*) CR//trim(word)//" object,vertex ID:",mo,",",mj
   write(6,*) CR//" relPos (dx, dy, dz):",o(mo)%rp(:,mj)
   write(6,*) CR//" sphPos (dist, phi/pi, theta/pi):",o(mo)%sp(1,mj),o(mo)%sp(2:3,mj)/pi
   run=.false.; return
  case("underAim"); !report object ID and node ID being aimed at.
   write(6,*) CR//" Aiming upon object ID:",scr%centerObj, &
       &   " triangle ID:",scr%centerTri," vertex ID:",scr%centerNode
   if (scr%centerObj.gt.0) then
    write(6,*) CR//" sphPos (dist, phi/pi, theta/pi):", &
     &  o(scr%centerObj)%sp(1,scr%centerNode), &
     &  o(scr%centerObj)%sp(2:3,scr%centerNode)/pi
    write(6,*) CR//" color:",o(scr%centerObj)%color(:,scr%centerNode)
   endif
   run=.false.; return
  case ("saveObject","saveObj"); ! save this object as a separate object file
     if (wc==1) then; write(6,*) "Must specify an object ID number to save"
      run=.false.; return; endif
     call s_get_val( 2, s, i)
     if (wc==2) then; call saveObject( 20, "objectDump.obj", i )
     elseif (wc==3) then; call saveObject( 20, s_get_word(3,s), i ); endif
  case("cam","camera"); ! regarding the camera object
   if (wc==1) then
    write(6,*) CR//"name:"//trim( cam%name )
    write(6,*) CR//"mass:",       cam%mass 
    write(6,*) CR//"radius:",     cam%radius 
    write(6,*) CR//"stiffness:",  cam%stiffness
    write(6,*) CR//"poisson:",    cam%poisson
    write(6,*) CR//"friction:",   cam%friction
    !write(6,*) CR//"Npoints:",    cam%np       !number of points
    !write(6,*) CR//"Ntriangles:", cam%nt       !number of triangles
    !write(6,*) CR//"mode:"//trim( cam%mode )   !(point|wire|solid)
    write(6,*) CR//"offset:",     cam%po       !position offset wrt global coordinates
    write(6,*) CR//"velocity:",   cam%velocity !changes the offset position wtf global
    write(6,*) CR//"orientation:",scr%globalrotor   !orientation quaternion
    write(6,*) CR//"spin:",       cam%spin     !local rotation rate (changes orientation)
    run=.false.; return
   endif
   word2 = s_get_word(2, s)
   select case(trim(word2)) !{{{
    case ("name"); cam%name = s_get_word(3, s)
    case ("mass");   call s_get_val(3, s, cam%mass )
    case ("radius"); call s_get_val(3, s, cam%radius )
    case ("stiffness"); call s_get_val(3, s, cam%stiffness )
    case ("poisson"); call s_get_val(3, s, cam%poisson )
    case ("friction"); call s_get_val(3, s, cam%friction )
    case ("offset","pos","position"); 
     if (wc==1) then; write(6,*) cam%po; run=.false.; return; endif
     call s_get_val(3, s, cam%po(1) )
     call s_get_val(4, s, cam%po(2) )
     call s_get_val(5, s, cam%po(3) )
    case ("velocity","vel"); 
     if (wc==1) then; write(6,*) cam%velocity; run=.false.; return; endif
     call s_get_val(3, s, cam%velocity(1) )
     call s_get_val(4, s, cam%velocity(2) )
     call s_get_val(5, s, cam%velocity(3) )
    case ("orient"); !cam%orient should be constant, only globablrotor changes
     if (wc==1) then; write(6,*) scr%globalrotor; run=.false.; return; endif
     call s_get_val(3, s, scr%globalrotor(1) )
     call s_get_val(4, s, scr%globalrotor(2) )
     call s_get_val(5, s, scr%globalrotor(3) )
     call s_get_val(6, s, scr%globalrotor(4) )
     q1(1) = 0.0; q1(2:4) = scr%globalrotor(2:4) !axis portion
     call quatrotor( q1, scr%globalrotor(1), q2 )
     scr%globalrotor = q2
!     call s_get_val(4, s, cam%orient(1) )
!     call s_get_val(5, s, cam%orient(2) )
!     call s_get_val(6, s, cam%orient(3) )
!     call s_get_val(7, s, cam%orient(4) )
!     q1(1) = 0.0; q1(2:4) = cam%orient(2:4) !axis portion
!     call quatunit( q1, q2 )
!     call quatrotor( q2, cam%orient(1), q1 )
!     cam%orient = q1
    case ("spin"); 
     if (wc==1) then; write(6,*) cam%spin; run=.false.; return; endif
     call s_get_val(3, s, cam%spin(1) )
     call s_get_val(4, s, cam%spin(2) )
     call s_get_val(5, s, cam%spin(3) )
     call s_get_val(6, s, cam%spin(4) )
     q1(1) = 0.0; q1(2:4) = cam%spin(2:4) !axis portion
     call quatunit( q1, q2 )
     call quatrotor( q2, cam%spin(1), q1 )
     cam%spin = q1
   end select !}}}
  case("o"); ! regarding an object
   if (wc==1) then; write(6,*) "Nobjects:",Nobjects; run=.false.; return; endif
   call s_get_val(2, s, i) !object ID
   if (wc==2) then
    write(6,*) CR//"name:"//trim( o(i)%name )
    write(6,*) CR//"mass:",       o(i)%mass 
    write(6,*) CR//"radius:",     o(i)%radius 
    write(6,*) CR//"stiffness:",  o(i)%stiffness
    write(6,*) CR//"poisson:",    o(i)%poisson
    write(6,*) CR//"friction:",   o(i)%friction
    write(6,*) CR//"Npoints:",    o(i)%np       !number of points
    write(6,*) CR//"Ntriangles:", o(i)%nt       !number of triangles
    write(6,*) CR//"mode:"//trim( o(i)%mode )   !(point|wire|solid)
    write(6,*) CR//"offset:",     o(i)%po       !position offset wrt global coordinates
    !write(6,*) CR//"localPos:",   o(i)%point    !initial position wrt local coordinates
    !write(6,*) CR//"initPos:",    o(i)%ip       !initial position wrt global coordinates
    !write(6,*) CR//"relPos:",     o(i)%rp       !relative position wrt camera
    !write(6,*) CR//"sphPos:",     o(i)%sp/pi    !relative spherical position wrt camera
    write(6,*) CR//"velocity:",   o(i)%velocity !changes the offset position wtf global
    write(6,*) CR//"orientation:",o(i)%orient   !orientation quaternion
    write(6,*) CR//"spin:",       o(i)%spin     !local rotation rate (changes orientation)
    run=.false.; return
   endif
   word2 = s_get_word(3, s)
   select case(trim(word2)) !{{{
    case ("save"); ! save this object as a separate object file
     if (wc==3) then; call saveObject( 20, "objectDump.obj", i )
     elseif (wc==4) then; call saveObject( 20, s_get_word(4,s), i ); endif
    case ("name"); o(i)%name = s_get_word(4, s)
    case ("mode"); o(i)%mode = trim(s_get_word(4, s))
    case ("mass");   call s_get_val(4, s, o(i)%mass )
    case ("radius"); call s_get_val(4, s, o(i)%radius )
    case ("stiffness"); call s_get_val(4, s, o(i)%stiffness )
    case ("poisson"); call s_get_val(4, s, o(i)%poisson )
    case ("friction"); call s_get_val(4, s, o(i)%friction )
    case ("offset","pos","position"); 
     call s_get_val(4, s, o(i)%po(1) )
     call s_get_val(5, s, o(i)%po(2) )
     call s_get_val(6, s, o(i)%po(3) )
    case ("velocity","vel"); 
     call s_get_val(4, s, o(i)%velocity(1) )
     call s_get_val(5, s, o(i)%velocity(2) )
     call s_get_val(6, s, o(i)%velocity(3) )
    case ("orient"); 
     call s_get_val(4, s, o(i)%orient(1) )
     call s_get_val(5, s, o(i)%orient(2) )
     call s_get_val(6, s, o(i)%orient(3) )
     call s_get_val(7, s, o(i)%orient(4) )
     q1(1) = 0.0; q1(2:4) = o(i)%orient(2:4) !axis portion
     call quatunit( q1, q2 )
     call quatrotor( q2, o(i)%orient(1), q1 )
     o(i)%orient = q1
    case ("spin"); 
     call s_get_val(4, s, o(i)%spin(1) )
     call s_get_val(5, s, o(i)%spin(2) )
     call s_get_val(6, s, o(i)%spin(3) )
     call s_get_val(7, s, o(i)%spin(4) )
     q1(1) = 0.0; q1(2:4) = o(i)%spin(2:4) !axis portion
     call quatunit( q1, q2 )
     call quatrotor( q2, o(i)%spin(1), q1 )
     o(i)%spin = q1
    case ("scale"); !scale all object positions
     call s_get_val(4, s, scl )
     call o(i)%scale( scl )
     !do j = 1, o(i)%np
     ! o(i)%point(:,j) = o(i)%point(:,j)*scl
     !enddo
     !! also scale radius and mass (assuming constant density)
     !o(i)%radius = o(i)%radius*scl
     !o(i)%mass = o(i)%mass*scl*scl*scl
    case ("add"); !add a new point to triangle under AIM or ID specified
     if (wc==4) then !specified a triangle ID
      call s_get_val(4, s, j ) !point ID
     else ! try for the triangle under AIM
      j = scr%centerTri
      if (j.le.0) write(6,'(A)') CR//"no triangle under aim. Please specify triangle ID."
     endif
     if (j.le.o(i)%nt) then; call o(i)%addPoint( j ) 
     else; write(6,'(A)') CR//"triangle ID too large for this object"; endif
     run=.false.; return
    case ("point"); !change point location
     if (wc==3) then; write(6,'(A)') CR//"missing point ID for this object"
      run=.false.; return
     elseif (wc==4) then;
      call s_get_val(4, s, j ) !point ID
      write(6,*) CR//"local:",o(i)%point(:,j), o(i)%color(:,j), o(i)%smooth(j)
      write(6,*) CR//"initPos:",o(i)%ip(:,j)
      write(6,*) CR//"relPos:",o(i)%rp(:,j)
      write(6,*) CR//"sphPos:",o(i)%sp(1,j), o(i)%sp(2:3,j)/pi
      run=.false.; return
     elseif (wc>4) then;
      call s_get_val(4, s, j ) !point ID
      call s_get_val(5, s, p(1) )
      call s_get_val(6, s, p(2) )
      call s_get_val(7, s, p(3) )
      o(i)%point(:,j) = p
     endif
    case ("color"); !change point color
     if (wc==3) then; write(6,'(A)') CR//"missing point ID for this object (0=all)"
      run=.false.; return
     elseif (wc==4) then;
      call s_get_val(4, s, j ) !point ID
      write(6,*) o(i)%point(:,j), o(i)%color(:,j), o(i)%smooth(j)
      run=.false.; return
     elseif (wc>4) then;
      call s_get_val(4, s, j ) !point ID
      call s_get_val(5, s, c(1) )
      call s_get_val(6, s, c(2) )
      call s_get_val(7, s, c(3) )
      if (j.gt.0) then; o(i)%color(1:3,j) = c(1:3)
      else !change al point's color
       do j = 1, o(i)%np
        o(i)%color(1:3,j) = c(1:3)
       enddo
      endif
     endif
    case ("smooth"); !change point to have vertex normal
     if (wc==3) then; write(6,'(A)') CR//"missing point ID for this object (0=all)"
      run=.false.; return
     elseif (wc==4) then;
      call s_get_val(4, s, j ) !point ID
      write(6,*) o(i)%point(:,j), o(i)%color(:,j), o(i)%smooth(j)
      run=.false.; return
     elseif (wc>4) then;
      call s_get_val(4, s, j ) !point ID
      call s_get_val(5, s, k ) ! set 0 for sharp, set 1 for smooth
      if (j.gt.0) then; 
       if (k.eq.1.and.o(i)%smooth(j).ne.1) then
        o(i)%nvn = o(i)%nvn + 1
        call o(i)%findVertexNormal(j)
       elseif (k.ne.1.and.o(i)%smooth(j).eq.1) then
        o(i)%nvn = o(i)%nvn - 1
       endif
       o(i)%smooth(j) = k
      else ! for every point
       do j = 1, o(i)%np
        if (k.eq.1.and.o(i)%smooth(j).ne.1) then
         o(i)%nvn = o(i)%nvn + 1
         call o(i)%findVertexNormal(j)
        elseif (k.ne.1.and.o(i)%smooth(j).eq.1) then
         o(i)%nvn = o(i)%nvn - 1
        endif
        o(i)%smooth(j) = k
       enddo
      endif
     endif
    case ("disp"); !change point location via added displacement vector
     if (wc==3) then; write(6,'(A)') CR//"missing point ID for this object"
      run=.false.; return
     elseif (wc==4) then;
      call s_get_val(4, s, j ) !point ID
      write(6,*) o(i)%point(:,j), o(i)%color(:,j), o(i)%smooth(j)
      run=.false.; return
     elseif (wc>4) then;
      call s_get_val(4, s, j ) !point ID
      call s_get_val(5, s, p(1) )
      call s_get_val(6, s, p(2) )
      call s_get_val(7, s, p(3) )
      o(i)%point(:,j) = o(i)%point(:,j) + p
     endif
    case default
     write(6,'(A)') CR//"object option not recognized:"//trim(word2)
   end select !}}}
 
  case default; write(6,'(A)') "don't have interpretation for: '"//trim(word)//"'"//achar(13)
end select

end subroutine interpretLine !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
function objectEnergy( o ) result(EE) !{{{
! returns the translational and angular kinetic energy (1:2)
use quaternion
class(object), intent(in) :: o
real(4),  dimension(2) :: EE !translational and angular kinetic energy
real(4), dimension(4) :: q

call rotoraxis( o%spin, q ) !get angle rate [and axis]
EE(1) = 0.5*o%mass*abs(dot_product(o%velocity,o%velocity))
EE(2) = 0.5*0.4*o%mass*o%radius*o%radius*q(1)*q(1) 

end function objectEnergy !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
function objectMomentum( o ) result(EE) !{{{
! returns the translational and angular kinetic energy, and total (1:8)
use quaternion
class(object), intent(in) :: o
real(4),  dimension(8) :: EE !translational and angular kinetic energy
real(4), dimension(4) :: q
real(4), dimension(3) :: v

call rotoraxis( o%spin, q ) !get angle rate [and axis]
v = o%velocity
EE(1:3) = o%mass*v !m*v
EE(7) = sqrt(dot_product(v,v)) !magntude
EE(8) = 0.4*o%mass*o%radius*o%radius*q(1) !I*omega
EE(4:6) = EE(8)*q(2:4) ! angular vector

end function objectMomentum !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine objectScale( o, scl )  !{{{
! scales object lengths by scl and mass by scl**3
!use quaternion
class(object) :: o
real(4), intent(in) :: scl
integer :: j
do j = 1, o%np
 o%point(:,j) = o%point(:,j)*scl
enddo
! also scale radius and mass (assuming constant density)
o%radius = o%radius*scl
o%mass = o%mass*scl*scl*scl
end subroutine objectScale !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine sortNodes( o ) !{{{
! sorts the node distances far to near, sorted IDs in sortID
class(object) :: o
integer :: i !, swap1, t
!integer, dimension(1) :: swap
!real(4) :: temp
real(4), allocatable, dimension(:) :: p

if (.not.allocated(o%sortID)) allocate( o%sortID(o%np) )

allocate( p(o%np) )
p = o%sp(1,:) ! copy distances to this temp array, p
do i = 1, o%np
  o%sortID(i) = i ! set the index ID
enddo

call hpsort_eps_epw( o%np, p, o%sortID, 1.e-10)

! Sort algorithm by John Mahaffy March 10, 1995
!do i = 1, o%np ! sort minimum distance to greatest
!  !swap = MINLOC( o%sp(1,i:o%np) ) ! sort by distance
!  swap = MAXLOC( p(i:o%np) ) ! sort by distance
!  swap1 = swap(1) +i -1
!  if (swap1.ne.i) then
!    temp = p(i)
!    p(i) = p(swap1)
!    p(swap1) = temp
!    t = o%sortID(i)
!    o%sortID(i) = o%sortID(swap1)
!    o%sortID(swap1) = t
!  endif
!enddo

!do i = 1, o%np
!  write(0,*) o%sp(1,i), p(i), o%sortID(i)
!enddo

deallocate( p )
end subroutine sortNodes !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine findTriangleNormal( o, t, norm ) !{{{
! calculates the triangle normal for triangle t in local initial object coordinates
implicit none
class(object) :: o
integer, intent(in) :: t
integer :: i, j, k
real(4) :: dist
real(4), dimension(3), intent(out) :: norm
real(4), dimension(3) :: aa, bb

 i = o%tri(1,t) !vertex ID
 j = o%tri(2,t) !vertex ID
 k = o%tri(3,t) !vertex ID
  ! calculate normal vector for this triangular face (cross product)
  aa(:) = o%point(:,i) - o%point(:,j) !vertex i to j
  bb(:) = o%point(:,i) - o%point(:,k) !vertex i to k
  norm(1) = aa(2)*bb(3) - aa(3)*bb(2)
  norm(2) = aa(3)*bb(1) - aa(1)*bb(3)
  norm(3) = aa(1)*bb(2) - aa(2)*bb(1)
  dist = sqrt(dot_product(norm,norm)) !make unit vector
  if (dist.lt.1.e-12) then
   norm = 0.0; Return
  endif ! collapsed triangle (i.e. line) undefined normal
  norm(:) = norm(:)/dist !make unit vector

end subroutine findTriangleNormal !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine findVertexNormal( o, id ) !{{{
! calculates the vertex normal for point ID
! during shading this acts to smooth the vertex
! NOTE: this routine is not designed to be run during animation, too costly!
implicit none
class(object) :: o
integer, intent(in) :: id
integer :: t, nt
real(4), dimension(3) :: norm, nnorm, normsum
real(4) :: a1, a2
logical :: doit

normsum = 0.0
nt = 0 ! number of adjacent triangles
! must loop over triangles to find those that share this vertex
do t = 1, o%nt
 doit = .false.
 ! this construction saves a little time searching all three points of t
 if (.not.( id .eq. o%tri(1,t))) then !vertex ID
  if (.not.( id .eq. o%tri(2,t))) then !vertex ID
   if (.not.( id .eq. o%tri(3,t))) then !vertex ID
    cycle
   else; doit = .true.; endif
  else; doit = .true.; endif
 else; doit = .true.; endif

 if (doit) then
  nt = nt + 1
  ! calculate normal vector for this triangular face (cross product)
  call o%findTriangleNormal( t, nnorm )
  if (nt.gt.1) then ! compare angle
   a1 = dot_product(norm,nnorm)
   a2 = dot_product(norm,-nnorm)
   if (a2.gt.a1) nnorm = -nnorm
  endif
  norm = nnorm ! this is th enew old normal
  normsum = normsum + norm
 endif
enddo

! normal average
o%vertnorm(:,id) = normsum/real(nt)

! this below works for a sphere
!norm=o%point(:,id)
!o%vertnorm(:,id) = norm/sqrt(dot_product(norm,norm))

end subroutine findVertexNormal !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine objectReallocate( o, n ) !{{{
class(object) :: o
integer,intent(in) :: n
integer :: nn, dn
real(4), allocatable, dimension(:,:) :: p
integer, allocatable, dimension(:,:) :: c
integer, allocatable, dimension(:) :: s

! new points
dn = n -o%np

! point positions
allocate( p(3,n) )
 p(1:3,1:o%np) = o%point(1:3,1:o%np)
 p(1:3,o%np+1:n) = 0.0
call move_alloc(p,o%point)

if (o%nvn.gt.0) then
 allocate( p(3,n) )
  p(1:3,1:o%np) = o%vertnorm(1:3,1:o%np)
  p(1:3,o%np+1:n) = 0.0
 call move_alloc(p,o%vertnorm)
 allocate( p(3,n) )
  p(1:3,1:o%np) = o%ivn(1:3,1:o%np)
  p(1:3,o%np+1:n) = 0.0
 call move_alloc(p,o%ivn)
 allocate( p(3,n) )
  p(1:3,1:o%np) = o%rvn(1:3,1:o%np)
  p(1:3,o%np+1:n) = 0.0
 call move_alloc(p,o%rvn)
endif

allocate( p(3,n) )
 p(1:3,1:o%np) = o%ip(1:3,1:o%np)
 p(1:3,o%np+1:n) = 0.0
call move_alloc(p,o%ip)

allocate( p(3,n) )
 p(1:3,1:o%np) = o%rp(1:3,1:o%np)
 p(1:3,o%np+1:n) = 0.0
call move_alloc(p,o%rp)

allocate( p(3,n) )
 p(1:3,1:o%np) = o%sp(1:3,1:o%np)
 p(1:3,o%np+1:n) = 0.0
call move_alloc(p,o%sp)

! point color
allocate( c(3,n) )
 c(1:3,1:o%np) = o%color(1:3,1:o%np)
 c(1:3,o%np+1:n) = 0
call move_alloc(c,o%color)

! point smooth parameter (or Aux)
allocate( s(n) )
 s(1:o%np) = o%smooth(1:o%np)
 s(o%np+1:n) = 0
call move_alloc(s,o%smooth)

! object triangles
nn = o%nt +2*dn
allocate( c(3,nn) ) !2 new triangles per new point
 c(1:3,1:o%nt) = o%tri(1:3,1:o%nt)
 c(1:3,o%nt+1:nn) = 0
call move_alloc(c,o%tri)

end subroutine objectReallocate !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine objectAddPoint( o, tri ) !{{{
! add a new node to the center of triangle tri.
!      A
!      d        d is new node
!   C     B
! triangle ABC becomes ABd, with new triangles: BCd and CAd
implicit none
class(object) :: o
integer, intent(in) :: tri
integer :: A, b, C, d, t

d = o%np+1
call o%reallocate( d ) !reallocate arrays to accomodate this new node
o%np = d !new number of total points


if (d.gt.size(o%point(1,:))) then
 write(0,*) "ERROR: point ID too large",d, o%np, size(o%point(1,:))
 STOP
endif
! node IDs of the triangle to subdivide
A = o%tri(1,tri)
B = o%tri(2,tri)
C = o%tri(3,tri)

! find position of 'd', average of triangle nodes
o%point(:,d) = 0.3333333*(o%point(:,A) + o%point(:,B) + o%point(:,C))
o%ip(:,d) = 0.3333333*(o%ip(:,A) + o%ip(:,B) + o%ip(:,C))
o%rp(:,d) = 0.3333333*(o%rp(:,A) + o%rp(:,B) + o%rp(:,C))
o%sp(:,d) = 0.3333333*(o%sp(:,A) + o%sp(:,B) + o%sp(:,C))
o%smooth(d) = nint(0.3333333*(o%smooth(A) + o%smooth(B) + o%smooth(C)))
o%color(:,d) = nint(0.3333333*(o%color(:,A) + o%color(:,B) + o%color(:,C)))

! construct triangles
o%tri(3,tri) = d !ABd
t = o%nt+1
o%tri(:,t) = (/ B, C, d /)
t = o%nt+2
o%tri(:,t) = (/ C, A, d /)

o%nt = t !new number of total triangles

end subroutine objectAddPoint !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine PhongTess( p1,p2,p3, n1,n2,n3, u, b ) !{{{
! Point Normal triangle surface coordinate interpolation 
! p(1:3) = positions vecotrs
! n(1:3) = normal vectors
! u(1:3) = barycentric coordinates of point to interpolate
! b(1:3) = output the interpolated surface position
implicit none
real(4), dimension(3), intent(in) :: p1, p2, p3, n1, n2, n3, u
real(4), dimension(3), intent(out) :: b
real(4), dimension(3) :: p12, p21, p23, p32, p31, p13 

p12 = p2-dot_product(p2-p1,n1)*n1
p21 = p1-dot_product(p1-p2,n2)*n2
p23 = p3-dot_product(p3-p2,n2)*n2
p32 = p2-dot_product(p2-p3,n3)*n3
p31 = p1-dot_product(p1-p3,n3)*n3
p13 = p3-dot_product(p3-p1,n1)*n1

b = u(1)*u(1)*p1 +u(2)*u(2)*p2 +u(3)*u(3)*p3 &
  & +u(1)*u(2)*(p12+p21) +u(2)*u(3)*(p23+p32) &
  & +u(3)*u(1)*(p31+p13)

end subroutine PhongTess !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine PhongNorm( n1,n2,n3, u, n ) !{{{
! Point Normal triangle surface coordinate interpolation 
! n(1:3) = normal vectors
! u(1:3) = barycentric coordinates of point to interpolate
! b(1:3) = output the interpolated surface position
implicit none
real(4), dimension(3), intent(in) :: n1, n2, n3, u
real(4), dimension(3), intent(out) :: n
real(4), dimension(3) :: p12, p21, p23, p32, p31, p13 

p12 = n2-dot_product(n2-n1,n1)*n1
p21 = n1-dot_product(n1-n2,n2)*n2
p23 = n3-dot_product(n3-n2,n2)*n2
p32 = n2-dot_product(n2-n3,n3)*n3
p31 = n1-dot_product(n1-n3,n3)*n3
p13 = n3-dot_product(n3-n1,n1)*n1

n = u(1)*u(1)*n1 +u(2)*u(2)*n2 +u(3)*u(3)*n3 &
  & +u(1)*u(2)*(p12+p21) +u(2)*u(3)*(p23+p32) &
  & +u(3)*u(1)*(p31+p13)

! some reason this isn't normalized
n = n/sqrt(dot_product(n,n))

end subroutine PhongNorm !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine PNtriPos( p1,p2,p3, n1,n2,n3, u, b ) !{{{
! Point Normal triangle surface coordinate interpolation 
! p(1:3) = positions vecotrs
! n(1:3) = normal vectors
! u(1:3) = barycentric coordinates of point to interpolate
! b(1:3) = output the interpolated surface position
implicit none
real(4), dimension(3), intent(in) :: p1, p2, p3, n1, n2, n3, u
real(4), dimension(3), intent(out) :: b
real(4), dimension(3) :: b300, b030, b003
real(4), dimension(3) :: b210, b120, b021
real(4), dimension(3) :: b012, b102, b201, b111
real(4), dimension(3) :: E, V
real(4), dimension(3) :: w12, w21, w23, w32, w31, w13

! wij = (pj-pi)*Ni
w12 = (p2-p1)*n1; w21 = (p1-p2)*n2
w23 = (p3-p2)*n2; w32 = (p2-p3)*n3
w31 = (p1-p3)*n3; w13 = (p3-p1)*n1

b300 = p1; b030 = p2; b003 = p3

b210 = (2.0*p1 +p2 -w12*n1)/3.0
b120 = (2.0*p2 +p1 -w21*n2)/3.0
b021 = (2.0*p2 +p3 -w23*n2)/3.0
b012 = (2.0*p3 +p2 -w32*n3)/3.0
b102 = (2.0*p3 +p1 -w31*n3)/3.0
b201 = (2.0*p1 +p3 -w13*n1)/3.0

E = (b210 +b120 +b021 +b012 +b102 +b201)/6.0
V = (p1 +p2 +p3)/3.0
b111 = 0.5*(E-V)

! ideally to raster an interpolated surface these coefficients should be constant
! and you'd put a loop here for the local barycentric coordinates
b = b300*u(3)*u(3)*u(3) + b030*u(1)*u(1)*u(1) + b003*u(2)*u(2)*u(2) + &
  & b210*3.0*u(3)*u(3)*u(1) + b120*3.0*u(3)*u(1)*u(1) + b201*3.0*u(3)*u(3)*u(2) + &
  & b021*3.0*u(1)*u(1)*u(2) + b102*3.0*u(3)*u(2)*u(2) + b111*6.0*u(3)*u(2)*u(1)


end subroutine PNtriPos !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine PNtriNorm( p1,p2,p3, n1,n2,n3, u, n ) !{{{
! Point Normal triangle surface coordinate interpolation 
! p1:3 = positions vecotrs
! n1:3 = normal vectors
! u(1:3) = barycentric coordinates of point to interpolate
! n(1:3) = output the interpolated surface normal
implicit none
real(4), dimension(3), intent(in) :: p1, p2, p3, n1, n2, n3, u
real(4), dimension(3), intent(out) :: n
real(4), dimension(3) :: n200, n020, n002
real(4), dimension(3) :: n110, n011, n101
real(4) :: v12, v23, v31

! vij = 2.0*(pj-pi).cdot.(ni+nj)/(pj-pi).cdot.(pj-pi)
v12 = 2.0*dot_product((p2-p1),(n1+n2))/dot_product((p2-p1),(p2-p1))
v23 = 2.0*dot_product((p3-p2),(n2+n3))/dot_product((p3-p2),(p3-p2))
v31 = 2.0*dot_product((p1-p3),(n1+n2))/dot_product((p2-p1),(p2-p1))

n200 = n1; n020 = n2; n002 = n3
n110 = (n1 +n2 -v12*(p2-p1))
n011 = (n2 +n3 -v23*(p3-p2))
n101 = (n3 +n1 -v31*(p1-p3))
n110 = n110/sqrt(dot_product(n110,n110))
n011 = n011/sqrt(dot_product(n011,n011))
n101 = n101/sqrt(dot_product(n101,n101))

! ideally to raster an interpolated surface these coefficients should be constant
! and you'd put a loop here for the local barycentric coordinates
n =  n200*u(3)*u(3) +n020*u(1)*u(1) +n002*u(2)*u(2) & 
  & +n110*u(3)*u(1) +n011*u(1)*u(2) +n101*u(3)*u(2)

! some reason this isn't normalized
n = n/sqrt(dot_product(n,n))

end subroutine PNtriNorm !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine closestTriangle( o, tri, w ) !{{{
implicit none
class(object) :: o
integer, intent(out) :: tri
integer :: ca(1), cv, i, t
integer :: iA
real(4) :: dA, ra(1), x1,y1,x2,y2,x3,y3
real(4), dimension(3), intent(out) :: w

! closest vertex in this object
ca = minloc(o%sp(1,:)) 
cv = ca(1)

! check for vertec to triangle connexions
if (.not.allocated(o%connect)) then
 call o%mkconnections
endif 

! iA is the first vertex ID of the first triangle that has cv
!iA = o%tri(1,o%connect(cv,1)); iB=0
!dA = o%sp(1,cv); dB=0.0

iA = o%connect(cv,1) !; iB = 0.0
dA = 0.0
! loop over the connecting triangles
do i = 1, o%connect(cv,0)
 t = o%connect(cv,i) ! triangle ID
 x1 = o%ip(1,o%tri(1,t))
 y1 = o%ip(2,o%tri(1,t))
 x2 = o%ip(1,o%tri(2,t))
 y2 = o%ip(2,o%tri(2,t))
 x3 = o%ip(1,o%tri(3,t))
 y3 = o%ip(2,o%tri(3,t))
 call getBary( x1,y1, x2,y2, x3,y3, cam%po(1), cam%po(2), w )
 ra = minval(w) 
!write(0,*) "  t,w:",t,w, char(13)
 ! find triangle that has the largest barycentric coordinates
 ! negative coordinates are outside the triangle
 if (ra(1).gt.dA) then
   dA= ra(1); iA = t
 endif

! do j = 1, 3
!  v = o%tri(j,t) ! vector ID
!  if (v.eq.cv) cycle 
!  ! collect the next closest tris
!  if (o%sp(1,v).lt.dA) then
!   dB = dA; iB = iA ! save current tri
!   dA = o%sp(1,v); iA = t ! update
!  endif
! enddo
enddo
tri = iA

if (dA.lt.0.0) then !still an issue
  write(0,*) "ERROR: problem with closest triangle.",char(13)
  write(0,*) "point not projected inside any nearest triangle"
endif

end subroutine closestTriangle !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine getBary3d( p1,p2,p3, t, w ) !{{{
! provide a 3d vector position, t
! provide triangle vertec positions, p1,p2,p3
! returns Barycentric coordinate of, t, as w
implicit none
real(4), dimension(3), intent(in) :: p1, p2, p3, t
real(4), dimension(3), intent(out) :: w
real(4), dimension(3) :: v0, v1, v2
real(4) :: d00, d01, d11, denom, d20, d21
v0 = p1-p2 
v1 = p3-p2
v2 = t-p2

d00 = dot_product( v0, v0 )
d01 = dot_product( v0, v1 )
d11 = dot_product( v1, v1 )
denom = d00*d11 -d01*d01
d20 = dot_product( v2, v0 )
d21 = dot_product( v2, v1 )

w(1) = (d11*d20 -d01*d21)/denom
w(2) = (d00*d21 -d01*d20)/denom
w(3) = 1 -w(1) -w(2)
end subroutine getBary3d !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine getBary( x1,y1, x2,y2, x3,y3,  x,y, w ) !{{{
implicit none
real(4), intent(in) :: x1,x2,x3,y1,y2,y3,x,y
real(4), dimension(3), intent(out) :: w
real(4) :: y23, x32, x13, y13, y31, denom, x03, y03

y23 = y2-y3
x32 = x3-x2
x13 = x1-x3
y13 = y1-y3
y31 = y3-y1
denom = y23*x13 + x32*y13
x03 = x-x3
y03 = y-y3

w(1) = y23*x03 + x32*y03
w(1) = w(1)/denom
w(2) = y31*x03 + x13*y03
w(2) = w(2)/denom
w(3) = 1.0 -w(1) -w(2)

end subroutine getBary !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine mkconnections( o ) !{{{
implicit none
class(object) :: o
integer :: i, t, id, j

if (.not.allocated(o%connect)) allocate( o%connect(o%np,0:6) )
o%connect = 0
do t = 1, o%nt ! loop triangles
 do i = 1, 3 !triangle points
  id = o%tri(i,t) !vertex ID
  j = o%connect(id,0) + 1
  if (j.gt.6) then
    write(0,*) "ERROR: overflow, too many triangles connected to vertex: ",id
  endif
  o%connect(id,0) = j ! record current connected triangles
  o%connect(id,j) = t ! save triangle ID
 enddo
enddo

end subroutine mkconnections !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine objectDynamics !{{{
use quaternion
implicit none
integer :: i, n
real(4) :: dist, rp(3), nv(3), dv(3)
real(4), dimension(4) :: qtmp2

!stuff that happens independent of key strokes, like environmental dynamics
! rotate the cumuative local rotor by the rotor rate qlrr

! deal with background object first and separatly since it doesn't translate
if (o(0)%np.gt.0) then !there is a background object
  n = 0 !background object ID
! if you want backgrounds to spin on their own...
!  call quatmult( o(n)%spin, o(n)%orient, qtmp2)
!  o(n)%orient = qtmp2
  do i = 1, o(n)%np
   ! local points rotation to global coordinates
!   qtmp1(1) = 0.0; qtmp1(2:4) = o(n)%point(1:3,i)
!   call quatrotate( qtmp1, o(n)%orient, qtmp2)
!   o(n)%ip(:,i) = qtmp2(2:4) +o(n)%po(:)
   
   ! points global coordinates to camera coordinates
   !o(n)%rp(:,i) = o(n)%ip(:,i) - scr%camera_p(:)
   !qtmp1(1) = 0.0; qtmp1(2:4) = o(n)%ip(1:3,i)
   !call quatrotate( qtmp1, scr%globalrotor, qtmp2)
   call quat3rotate( o(n)%ip(:,i), scr%globalrotor, rp)
   !o(n)%rp(:,i) = qtmp2(2:4)
   !call cart2sph( qtmp2(2:4), o(n)%sp(:,i))
   call cart2sph( rp, o(n)%sp(:,i))
  enddo !end point loop
endif

! deal with all other objects as near enough for perceptable traslation
do n = 1, Nobjects

  call quatmult( o(n)%spin, o(n)%orient, qtmp2)
  o(n)%orient = qtmp2

  ! reorient velocity vector
  dv = o(n)%velocity
  ! integrate position of object in space.
  o(n)%po(:)   = o(n)%po(:) + dv*scr%dt
  do i = 1, o(n)%np
   ! if the object its self is rotating it can be implemented here
   !  based on its original coordinates.
   !   ip=qlrr*lp*qlrr^-1 + po
   !   rp=qgr*(ip-cam)*qgr^-1
   ! local points rotation to global coordinates
   call quat3rotate( o(n)%point(:,i), o(n)%orient, rp)
   o(n)%ip(:,i) = rp +o(n)%po(:)
   
   ! points global coordinates to camera coordinates
   call quat3rotate( o(n)%ip(:,i)-cam%po, scr%globalrotor, o(n)%rp(:,i))
   
   call cart2sph(o(n)%rp(:,i), o(n)%sp(:,i))
   ! transform the vertex normals as well.
   if (o(n)%smooth(i).eq.1) then
    ! local points rotation to global coordinates
    call quat3rotate( o(n)%vertnorm(:,i), o(n)%orient, o(n)%ivn(:,i))
    
    ! points global coordinates to camera coordinates
    call quat3rotate( o(n)%ivn(:,i), scr%globalrotor, o(n)%rvn(:,i))
   endif
  enddo !end point loop
enddo !end object loop

! test for universe geometry
select case (trim(scr%UniGeom))
 case ("sphere") 
  do n = 1, Nobjects
   !rp = o(n)%po -scr%camera_p
   rp = o(n)%po -cam%po
   !dist = dot_product(o(n)%po,o(n)%po) ! distance from origin
   dist = dot_product(rp,rp) !relative position distance
   if (scr%periodic) then
    if (dist.gt.scr%UniParms(1)*scr%UniParms(1)) then !fold back
     !o(n)%po = -o(n)%po*scr%UniParms(1)**2/dist !negate and scale down
     !o(n)%po = -rp*scr%UniParms(1)**2/dist !negate and scale down
     !o(n)%po = -(rp)*scr%UniParms(1)**2/dist +scr%camera_p !negate and scale down
     o(n)%po = -(rp)*scr%UniParms(1)**2/dist +cam%po !negate and scale down
    endif
   else !reflect
    if (dist.gt.scr%UniParms(1)*scr%UniParms(1)) then !fold back
     ! r = d-2(d.n).n
     nv(:) = -rp/sqrt(dist) !normal vector to reflexion plane
     dv(:) = o(n)%velocity(:)
     o(n)%velocity(:) = dv -2.0*dot_product(dv,nv)*nv
     ! scale position back to the boundary
     !o(n)%po = rp*(scr%UniParms(1)*scr%UniParms(1))/(dist) +scr%camera_p
     o(n)%po = rp*(scr%UniParms(1)*scr%UniParms(1))/(dist) +cam%po
    endif
   endif
  enddo
 case ("box") 
  do n = 1, Nobjects
   !rp = o(n)%po -scr%camera_p
   rp = o(n)%po -cam%po
   if (scr%periodic) then
   ! this needs to be a displacement from the camera position
    if (abs(rp(1)).gt.scr%UniParms(1)) then !fold back
     !o(n)%po(1) = o(n)%po(1) -2.0*scr%UniParms(1)*o(n)%po(1)/abs(o(n)%po(1))
     o(n)%po(1) = o(n)%po(1) -2.0*scr%UniParms(1)*rp(1)/abs(rp(1))
    elseif (abs(rp(2)).gt.scr%UniParms(2)) then !fold back
     !o(n)%po(2) = o(n)%po(2) -2.0*scr%UniParms(2)*o(n)%po(2)/abs(o(n)%po(2))
     o(n)%po(2) = o(n)%po(2) -2.0*scr%UniParms(2)*rp(2)/abs(rp(2))
    elseif (abs(rp(3)).gt.scr%UniParms(3)) then !fold back
     !o(n)%po(3) = o(n)%po(3) -2.0*scr%UniParms(3)*o(n)%po(3)/abs(o(n)%po(3))
     o(n)%po(3) = o(n)%po(3) -2.0*scr%UniParms(3)*rp(3)/abs(rp(3))
    endif
   else !reflect
    if (abs(rp(1)).gt.scr%UniParms(1)) then 
     o(n)%velocity(1) = -o(n)%velocity(1) 
    elseif (abs(rp(2)).gt.scr%UniParms(2)) then 
     o(n)%velocity(2) = -o(n)%velocity(2) 
    elseif (abs(rp(3)).gt.scr%UniParms(3)) then 
     o(n)%velocity(3) = -o(n)%velocity(3) 
    endif
   endif
  enddo
end select

end subroutine objectDynamics !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine objectContact( ct, lc ) !{{{
! checks for contact among every object in object array: o(:)
! contact type (ct)
!  1  billiard ball contact. no angular momentum transfer, elastic contact
!  2  hertzian contact. angular momentum transfer, elastic contact
! DATE: May 19-24, 2019
! camera contact (cc)
!  0  true: include camera as contact object
!  1  false: do not include camaera object in contact
! DATE: Sept 23, 2019
use quaternion
implicit none
integer :: i, j, ct, cc
real :: d2, ri, rj, mi, mj, si, sj, dt, avgF, depth, Emod, rad, dwi, dwj, sfric
real :: Efric, dtimp, dti, dtj, COF
real, dimension(3) :: r2, n, ui, uj, vi, vj, vi2, vj2, fi, fj, fric, faxis
real, dimension(4) :: qi, qj, q1, q2, EE
real, dimension(8,2) :: MM
logical :: lc

! contact requires that objects have mass and an average object radius for a
! spherical approximation of the object. works fine for asteroids and moons, etc.

if (lc) then; cc = 0
else;         cc = 1; endif

do i = 1, Nobjects
   ri = o(i)%radius
   o(i)%collidedWithObjID = 0
   do j = 1, i-cc !if camera contact, use j for camera
      if (i.eq.j) then
       rj = cam%radius
       cam%collidedWithObjID = 0
       !displacement vector
       n(1) = o(i)%po(1)-cam%po(1)
       n(2) = o(i)%po(2)-cam%po(2)
       n(3) = o(i)%po(3)-cam%po(3)
      else
       rj = o(j)%radius
       o(j)%collidedWithObjID = 0
       !displacement vector
       n(1) = o(i)%po(1)-o(j)%po(1)
       n(2) = o(i)%po(2)-o(j)%po(2)
       n(3) = o(i)%po(3)-o(j)%po(3)
      endif
      ! square components
      r2(1) = n(1)*n(1)
      r2(2) = n(2)*n(2)
      r2(3) = n(3)*n(3)
      ! distance squared
      d2    = r2(1) +r2(2) +r2(3)
      ! if they're in contact
      if (d2.lt.(ri+rj)*(ri+rj)) then
        mi = o(i)%mass
        if (i.ne.j) then
          mj = o(j)%mass
          o(i)%collidedWithObjID = j
          o(j)%collidedWithObjID = i
        else;         
          mj =  cam%mass
          o(i)%collidedWithObjID = -1
          cam%collidedWithObjID = i
write(0,*) "collision w/cam:", i, sqrt(d2), ri, rj
        endif
        n  = n/sqrt(d2) !unit normal vector
        !perpendicular speed (used for collision dynamics)
        si = dot_product(o(i)%velocity, n)
        if (i.ne.j) then; sj = dot_product(o(j)%velocity, n)
        else;             sj = dot_product( cam%velocity, n); endif
        !perpendicular velocity vector
        vi = si*n
        vj = sj*n
        !tangential velocity vector (unchanged during frictionless elastic collision)
        ui = o(i)%velocity -vi
        if (i.ne.j) then; uj = o(j)%velocity -vj
        else;             uj =  cam%velocity -vj; endif
        select case(ct)
         !###################################################################
         case(1) !billiard ball
          ! perform a simple rebound: conservation of momentun (in every component)
          !  mi*vi() +mj*vj() = mi*vi() +mj*vj()
          ! conservation of kinetic energy
          !  mi*vi*vi + mj*vj*vj = mi*vi*vi + mj*vj*vj
          ! solve for final velocity vectors
          ! the tanjential contact velocity remains constant for both bodies
          !  ui() = ui()  &&  uj() = uj()
          ! only need to solve for the "head-on" component
           vi2 = (mi-mj)/(mi+mj)*vi + 2.0*mj/(mi+mj)*vj
           vj2 =  2.0*mi/(mi+mj)*vi +(mj-mi)/(mi+mj)*vj
! momentum debugging
!write(0,*) "before", mi*vi, mj*vj
!write(0,*) "after",  mi*vi2, mj*vj2
!write(0,*) "diff",  mi*(vi2-vi), mj*(vj2-vj)
          ! add the perpendicular and tangential vectors back together
           o(i)%velocity = vi2 +ui
           if (i.ne.j) then; o(j)%velocity = vj2 +uj
           else;             cam%velocity = vj2 +uj; endif
         !###################################################################
         case(2) !frictional Hertzian contact, angular momentum transfer
          ! head-on collision components don't affect friction
           vi2 = (mi-mj)/(mi+mj)*vi + 2.0*mj/(mi+mj)*vj
           vj2 =  2.0*mi/(mi+mj)*vi +(mj-mi)/(mi+mj)*vj
          ! find friction vector common to both bodies
          ! cross the normal vector with the rotation axis to get -friction direction
           call rotoraxis( o(i)%spin, qi )
           if (i.ne.j) then; call rotoraxis( o(j)%spin, qj )
           else;             call rotoraxis(  cam%spin, qj ); endif
!write(0,*) "iaxis, jaxis",qi(2:4), qj(2:4)
          ! effective friction coefficient
           if (i.ne.j) then
                 COF = 2.0/( 1.0/o(i)%friction + 1.0/o(j)%friction )
           else; COF = 2.0/( 1.0/o(i)%friction + 1.0/cam%friction ); endif
           call cross_product(n, qi(1)*qi(2:4), fi);   fi = -fi
           call cross_product(n, qj(1)*qj(2:4), fj);   fj = -fj
           fric = fi + fj ! mutual friction vector
          ! calculate friction axis, fric-cross-n
           call cross_product(n, fric, Faxis)
           sfric = sqrt(dot_product(Faxis,Faxis))
if (abs(sfric).lt.1e-12) then
 write(0,*) "ERROR: frictional axis zero.:", Faxis, sfric, i, j 
 write(0,*) "normal,:", n
 write(0,*) "fric,:", fric
 write(0,*) "qi:", qi
 write(0,*) "qj:", qj
 write(0,*) "spini:", o(i)%spin
 if (i.ne.j) then; write(0,*) "spinj:", o(j)%spin
 else;             write(0,*) "spinj:", cam%spin; endif
 call flush(0)
endif
           Faxis = Faxis/sfric !unit vector
!write(0,*) "Faxis(3)",Faxis
           sfric = sqrt(dot_product(fric,fric))
           fric  = fric/sfric !unit vector
!write(0,*) "fric(3)",fric
          ! effective stiffness
           if (i.ne.j) then
           Emod = 1.0/( (1.0-o(i)%poisson*o(i)%poisson)/o(i)%stiffness &
                     & +(1.0-o(j)%poisson*o(j)%poisson)/o(j)%stiffness  )
           else
           Emod = 1.0/( (1.0-o(i)%poisson*o(i)%poisson)/o(i)%stiffness &
                     & +(1.0-cam%poisson*cam%poisson)/cam%stiffness  )
           endif
          ! effective radius
           rad = 1.0/( 1.0/ri + 1.0/rj )
          ! max depth of hertzian contact, found from equating kinetic energy to deformation.
          !  = 15^(2/5)/(4*2^(3/5))*(KE/EsqrtR)^(2/5)
           depth = 0.487257*( 0.5*(mi*si*si+mj*sj*sj)/(Emod*sqrt(rad)) )**(0.4)
           if (depth.gt.min(ri,rj)) then !too much energy
              depth = min(ri,rj) !or use rad
           endif
          ! average normal force of impulse. 16/15*(..), from force-distance integration
           avgF = 0.533333*Emod*sqrt(rad*depth*depth*depth)
          ! give friction the correct amplitude
           sfric = avgF*COF
           fric = fric*sfric
          ! duration of impulse, dt
          ! change in linear momentum / average normal force of impulse
           dtimp = mi*sqrt(dot_product(vi2-vi,vi2-vi)) / avgF
           !dtimp = mj*sqrt(dot_product(vj2-vj,vj2-vj)) / avgF !j and i should be equal.
EE(1:2)=o(i)%ene()
!EE(3:4)=o(j)%ene()
!write(0,*) achar(13)//"Elin Eang (i,j) tot:",EE(1),EE(3),EE(2),EE(4),sum(EE)
MM(:,1)=o(i)%momentum()
!MM(:,2)=o(j)%momentum()
!write(0,*) achar(13)//"Mlin (i,j) mag, tot:",MM(1:3,1),MM(1:3,2),MM(7,1),MM(7,2),MM(1:3,1)+MM(1:3,2)
!write(0,*) achar(13)//"Mang (i,j) mag, tot:",MM(4:6,1),MM(4:6,2),MM(8,1),MM(8,2),MM(4:6,1)+MM(4:6,2)
           if (sfric.gt.1e-12) then
            ! calculate energy transfer due to friction, Ffric*dist
             call rotoraxis( o(i)%spin, q1 )
             if (i.ne.j) then; call rotoraxis( o(j)%spin, q2 )
             else;             call rotoraxis( cam%spin, q2 ); endif
             dti = 0.4*mi*ri*abs(q1(1))/sfric
             dtj = 0.4*mj*rj*abs(q2(1))/sfric
            ! duration of frictive forces (could be less than impact duration)
             dt = min(dtimp, max(dti,dtj) )
            ! frictional energy is equivalent to the change in angular energy.
             Efric = 0.5*sfric*dt*abs(ri*q1(1) +rj*q2(1))
            ! subtract the tangential linear KE
             Efric = Efric -0.5*sfric*sfric*dt*dt*(1.0/mi+1.0/mj)
!write(0,*) achar(13)//"Frictional energy: ",Efric, "dt_imp:",dtimp, "dt(i:j):",dti,dtj
             if (dtimp.gt.scr%dt) then
              write(0,*) achar(13)//"WARNING: collision of objects:",i,j, &
                &  " duration longer than timestep:",dtimp,scr%dt
              write(0,*) achar(13)//" collision depth:",depth, "stiffness likely too low."
             endif
            ! add the frictional force to the tangential linear velocities
             ui = ui + fric/mi*dt 
             uj = uj - fric/mj*dt 
            ! change in angular velocity for i and j
            !  dw = tourque*dt/inertialmoment
             dwi = -dt/(0.4*mi*ri)*sfric ! angular rotation
             dwj = -dt/(0.4*mj*rj)*sfric ! angular rotation
             q1(1) = 0.0; q1(2:4) = Faxis
             call quatrotor(q1, dwi, qi) ! frictional rotor for i
             call quatrotor(q1, dwj, qj) ! frictional rotor for j
            ! multiply the frictional rotors to the objects' spin rotors
             call quatmult( o(i)%spin, qi, q1)  
            ! update the spin rotors
             o(i)%spin = q1
             if (i.ne.j) then
                call quatmult( o(j)%spin, qj, q2)  
                o(j)%spin = q2
             else;              
                call quatmult( cam%spin, qj, q2)
                cam%spin = q2
             endif
           endif
          ! add the perpendicular and tangential vectors back together
           o(i)%velocity = vi2 +ui
           if (i.ne.j) then; o(j)%velocity = vj2 +uj
           else;             cam%velocity = vj2 +uj; endif
!EE(1:2)=o(i)%ene()
!EE(3:4)=o(j)%ene()
!write(0,*) achar(13)//"Elin Eang (i,j) tot:",EE(1),EE(3),EE(2),EE(4),sum(EE)
!MM(:,1)=o(i)%momentum()
!MM(:,2)=o(j)%momentum()
!write(0,*) achar(13)//"Mlin (i,j) mag, tot:",MM(1:3,1),MM(1:3,2),MM(7,1),MM(7,2),MM(1:3,1)+MM(1:3,2)
!write(0,*) achar(13)//"Mang (i,j) mag, tot:",MM(4:6,1),MM(4:6,2),MM(8,1),MM(8,2),MM(4:6,1)+MM(4:6,2)
         !###################################################################
         case default
          write(0,*) "ERROR: unsupported contact type:",ct
          Return
        end select
      endif
   enddo
enddo

end subroutine objectContact !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine ObjectRandom( typ, i, j, v ) !{{{
use quaternion
implicit none
character(*) :: typ
integer, intent(in) :: i, j
real(4), dimension(8) :: v
integer :: n, nn, ms !,m
!integer, allocatable,i dimension(:) :: ss
real(4), allocatable, dimension(:) :: rands
real(4), dimension(4) :: q1, q2
!real(4) :: dist

!mass r          set object's total mass (for collision)"
!radius r        set object's average radius (for contact)"
!offset x y z    set position offset of object 'I'"
!velocity u v w  set velocity of object 'I'"
!orient a u v w  set angle and rotation axis of object 'I'"
!spin a u v w    set angle rate and rotation axis of object 'I'"
!call random_seed(SIZE=k)
!allocate( ss(k) )
!ss = 0
!ss(1) = scr%seed(1)
ms = scr%seed(1)
call random_seed(ms)

nn = j-i
select case(trim(typ))
 case("mass");
  allocate( rands(nn+1) )
  call random_number(rands)
  do n = i, j
   o(n)%mass = rands(n-i+1)*(v(2)-v(1)) + v(1)
  enddo
 case("radius");
  allocate( rands(nn+1) )
  call random_number(rands)
  do n = i, j
   o(n)%radius = rands(n-i+1)*(v(2)-v(1)) + v(1)
  enddo
 case("offset");
  allocate( rands(nn*3+1) )
  call random_number(rands)
  do n = i, j
   o(n)%po(1) = rands(n-i+1)*(v(2)-v(1)) + v(1)
   o(n)%po(2) = rands(n-i+1+nn)*(v(4)-v(3)) + v(3)
   o(n)%po(3) = rands(n-i+1+nn*2)*(v(6)-v(5)) + v(5)
   ! check the positions already created for overlap
   !do m = i, n-1
   ! dist = dot_product(o(n)%po-o(m)%po,o(n)%po-o(m)%po)
   ! if (dist .lt. o(n)%radius*o(m)%radius) then !overlap
   ! endif
   !enddo
  enddo
 case("velocity");
  allocate( rands(nn*3+1) )
  call random_number(rands)
  do n = i, j
   o(n)%velocity(1) = rands(n-i+1)*(v(2)-v(1)) + v(1)
   o(n)%velocity(2) = rands(n-i+1+nn)*(v(4)-v(3)) + v(3)
   o(n)%velocity(3) = rands(n-i+1+nn*2)*(v(6)-v(5)) + v(5)
  enddo
 case("orient");
  allocate( rands(nn*4+1) )
  call random_number(rands)
  do n = i, j
   o(n)%orient(1) = rands(n-i+1)*(v(2)-v(1)) + v(1)
   o(n)%orient(2) = rands(n-i+1+nn)*(v(4)-v(3)) + v(3)
   o(n)%orient(3) = rands(n-i+1+nn*2)*(v(6)-v(5)) + v(5)
   o(n)%orient(4) = rands(n-i+1+nn*3)*(v(8)-v(7)) + v(7)
   q1(1) = 0.0; q1(2:4) = o(i)%orient(2:4) !axis portion
   call quatunit( q1, q2 )
   call quatrotor( q2, o(i)%orient(1), q1 )
   o(n)%orient = q1
  enddo
 case("spin");
  allocate( rands(nn*4+1) )
  call random_number(rands)
  do n = i, j
   o(n)%spin(1) = rands(n-i+1)*(v(2)-v(1)) + v(1)
   o(n)%spin(2) = rands(n-i+1+nn)*(v(4)-v(3)) + v(3)
   o(n)%spin(3) = rands(n-i+1+nn*2)*(v(6)-v(5)) + v(5)
   o(n)%spin(4) = rands(n-i+1+nn*3)*(v(8)-v(7)) + v(7)
   q1(1) = 0.0; q1(2:4) = o(i)%spin(2:4) !axis portion
   call quatunit( q1, q2 )
   call quatrotor( q2, o(i)%spin(1), q1 )
   o(n)%spin = q1
  enddo
 case default
  write(0,*) "WARNING: unknown object randomizer type: "//trim(typ)
end select

deallocate(rands)

end subroutine ObjectRandom !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine drawObjects !{{{
use fbMod
implicit none
integer :: k, n, ii, i, t, v, j, ic(3), x(3), y(3), cx, cy
character(4) :: px, px1, px2, px3
real(4) :: tp(3,3), aa(3), bb(3), norm(3), trinorm(3), pi, ozb, dist
logical :: sorted

!pi = 4.0*ATAN2(1.0,1.0) !=log(-1.0)/i
pi = ACOS(-1.0)

scr%centerObj = 0
scr%centerTri = 0
scr%centerNode = 0
!screen center location
cx = scr%sr_x/2 + scr%sr(1,1)
cy = scr%sr_y/2 + scr%sr(2,1)
ozb = fb%zbuff(cx,cy) !old zbuff value at center of screen

! clear the subscreen with the background color
call fb%fillrec(scr%sr(1,1),scr%sr(2,1), scr%sr(1,2),scr%sr(2,2), scr%bgpx )

do n = 0, Nobjects

 if (o(n)%np.eq.0) cycle !object has no points and thus no triangles to render
 sorted = .false.
 !if (.not.scr%interactive) then ! will only order points for nonInteractive mode
   if (o(n)%transparency) then
     call o(n)%sortNodes
     sorted = .true.
   endif
 !endif

 select case(trim(o(n)%mode))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
 case("point") !{{{
  ! print points to screen if within Field of View (FOV)
  do ii = 1, o(n)%np
    if (sorted) then
      i = o(n)%sortID(o(n)%np-ii+1) ! sorted far to near
    else; i = ii; endif
    if (o(n)%sp(3,i).lt.-0.5*scr%FOVx.or.o(n)%sp(3,i).gt.0.5*scr%FOVx) cycle
    if (o(n)%sp(2,i).lt.0.5*(pi-scr%FOVy).or.o(n)%sp(2,i).gt.0.5*(scr%FOVy+pi)) cycle
    !ic(:) = int(o(n)%color(:,i)*1.0/max(o(n)%sp(1,i)-10.0,1.0)) !get contrast based on distance
    ic(:) = int(o(n)%color(:,i)) 
    if (o(n)%transparency) then
      px = char(ic(3))//char(ic(2))//char(ic(1))//char(nint(o(n)%vertnorm(2,i)))
    else
      px = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
    endif
    ! get screen position
    x(1) = int((o(n)%sp(3,i)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
    y(1) = int(scr%sr_y*((o(n)%sp(2,i)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
    !if (x(1).lt.0.or.x(1).gt.fb%width)  then; error=.true.; cycle; endif
    !if (y(1).lt.0.or.y(1).gt.fb%height) then; error=.true.; cycle; endif
    call fb%putPixel( x(1), y(1), px)
  enddo !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
 case("circ") !{{{
  ! print points to screen if within Field of View (FOV)
  do ii = 1, o(n)%np
    if (sorted) then
      i = o(n)%sortID(o(n)%np-ii+1) ! sorted far to near
    else; i = ii; endif
    if (o(n)%sp(3,i).lt.-0.5*scr%FOVx.or.o(n)%sp(3,i).gt.0.5*scr%FOVx) cycle
    if (o(n)%sp(2,i).lt.0.5*(pi-scr%FOVy).or.o(n)%sp(2,i).gt.0.5*(scr%FOVy+pi)) cycle
    !ic(:) = int(o(n)%color(:,i)*1.0/max(o(n)%sp(1,i)-10.0,1.0)) !get contrast based on distance
    ic(:) = int(o(n)%color(:,i)) 
    if (o(n)%transparency) then
      px = char(ic(3))//char(ic(2))//char(ic(1))//char(nint(o(n)%vertnorm(2,i)))
    else
      px = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
    endif
    ! get screen position
    x(1) = int((o(n)%sp(3,i)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
    y(1) = int(scr%sr_y*((o(n)%sp(2,i)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
    !if (x(1).lt.0.or.x(1).gt.fb%width)  then; error=.true.; cycle; endif
    !if (y(1).lt.0.or.y(1).gt.fb%height) then; error=.true.; cycle; endif
    !call fb%putPixel( x(1), y(1), px)
    ! calculate visual pixel distance corresponding to its radius
    !  ang*Dist=r so ang=r/dist
    ! point radius saved as vertnorm(1,i), divide by distance to camera
    x(2) = nint(o(n)%vertnorm(1,i)/(o(n)%sp(1,i)*scr%FOVy)*real(scr%sr_y))
!write(69,*) x(2)
    call fb%circle( x(1), y(1), x(2), px)
  enddo !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
 case("fcirc") !{{{
  ! print points to screen if within Field of View (FOV)
  do ii = 1, o(n)%np
    if (sorted) then
      i = o(n)%sortID(o(n)%np-ii+1) ! sorted far to near
    else; i = ii; endif
    if (o(n)%sp(3,i).lt.-0.5*scr%FOVx.or.o(n)%sp(3,i).gt.0.5*scr%FOVx) cycle
    if (o(n)%sp(2,i).lt.0.5*(pi-scr%FOVy).or.o(n)%sp(2,i).gt.0.5*(scr%FOVy+pi)) cycle
    !ic(:) = int(o(n)%color(:,i)*1.0/max(o(n)%sp(1,i)-10.0,1.0)) !get contrast based on distance
    ic(:) = int(o(n)%color(:,i)) 
    if (o(n)%transparency) then
      px = char(ic(3))//char(ic(2))//char(ic(1))//char(nint(o(n)%vertnorm(2,i)))
    else
      px = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
    endif
    ! get screen position
    x(1) = int((o(n)%sp(3,i)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
    y(1) = int(scr%sr_y*((o(n)%sp(2,i)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
    !if (x(1).lt.0.or.x(1).gt.fb%width)  then; error=.true.; cycle; endif
    !if (y(1).lt.0.or.y(1).gt.fb%height) then; error=.true.; cycle; endif
    !call fb%putPixel( x(1), y(1), px)
    ! calculate visual pixel distance corresponding to its radius
    !  ang*Dist=r so ang=r/dist
    ! point radius saved as vertnorm(1,i), divide by distance to camera
    x(2) = nint(o(n)%vertnorm(1,i)/(o(n)%sp(1,i)*scr%FOVy)*real(scr%sr_y))
    call fb%fillcircle( x(1), y(1), x(2), px, scr%sr, o(n)%sp(1,i) )
  enddo !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
 case("sph") !{{{
  ! print points to screen if within Field of View (FOV)
  do ii = 1, o(n)%np
    if (sorted) then
      i = o(n)%sortID(o(n)%np-ii+1) ! sorted far to near
    else; i = ii; endif
    if (o(n)%sp(3,i).lt.-0.5*scr%FOVx.or.o(n)%sp(3,i).gt.0.5*scr%FOVx) cycle
    if (o(n)%sp(2,i).lt.0.5*(pi-scr%FOVy).or.o(n)%sp(2,i).gt.0.5*(scr%FOVy+pi)) cycle
    !ic(:) = int(o(n)%color(:,i)*1.0/max(o(n)%sp(1,i)-10.0,1.0)) !get contrast based on distance
    ic(:) = int(o(n)%color(:,i)) 
    if (o(n)%transparency) then
      px = char(ic(3))//char(ic(2))//char(ic(1))//char(nint(o(n)%vertnorm(2,i)))
    else
      px = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
    endif
    ! get screen position
    x(1) = int((o(n)%sp(3,i)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
    y(1) = int(scr%sr_y*((o(n)%sp(2,i)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
    !if (x(1).lt.0.or.x(1).gt.fb%width)  then; error=.true.; cycle; endif
    !if (y(1).lt.0.or.y(1).gt.fb%height) then; error=.true.; cycle; endif
    !call fb%putPixel( x(1), y(1), px)
    ! calculate visual pixel distance corresponding to its radius
    !  ang*Dist=r so ang=r/dist
    ! point radius saved as vertnorm(1,i), divide by distance to camera
    x(2) = nint(o(n)%vertnorm(1,i)/(o(n)%sp(1,i)*scr%FOVy)*real(scr%sr_y))
    call fb%sphere( x(1), y(1), x(2), px, scr%sr, o(n)%sp(1,i) )
  enddo !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
 case("cloud") !{{{
  ! print points to screen if within Field of View (FOV)
  do ii = 1, o(n)%np
    if (sorted) then
      i = o(n)%sortID(o(n)%np-ii+1) ! sorted far to near
    else; i = ii; endif
    if (o(n)%sp(3,i).lt.-0.5*scr%FOVx.or.o(n)%sp(3,i).gt.0.5*scr%FOVx) cycle
    if (o(n)%sp(2,i).lt.0.5*(pi-scr%FOVy).or.o(n)%sp(2,i).gt.0.5*(scr%FOVy+pi)) cycle
    !ic(:) = int(o(n)%color(:,i)*1.0/max(o(n)%sp(1,i)-10.0,1.0)) !get contrast based on distance
    ic(:) = int(o(n)%color(:,i)) 
    if (o(n)%transparency) then
      px = char(ic(3))//char(ic(2))//char(ic(1))//char(nint(o(n)%vertnorm(2,i)))
    else
      px = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
    endif
    ! get screen position
    x(1) = int((o(n)%sp(3,i)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
    y(1) = int(scr%sr_y*((o(n)%sp(2,i)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
    !if (x(1).lt.0.or.x(1).gt.fb%width)  then; error=.true.; cycle; endif
    !if (y(1).lt.0.or.y(1).gt.fb%height) then; error=.true.; cycle; endif
    !call fb%putPixel( x(1), y(1), px)
    ! calculate visual pixel distance corresponding to its radius
    !  ang*Dist=r so ang=r/dist
    ! point radius saved as vertnorm(1,i), divide by distance to camera
    x(2) = nint(o(n)%vertnorm(1,i)/(o(n)%sp(1,i)*scr%FOVy)*real(scr%sr_y))
    call fb%cloud( x(1), y(1), x(2), px, scr%sr, o(n)%sp(1,i) )
  enddo !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
 case("wire") !{{{
  ! print edges(lines) to screen if within Field of View (FOV)
  do t = 1, o(n)%nt
   do v = 1, 3
    i = o(n)%tri(v,t) !vertex ID, i
    if (v.eq.3) then; j = o(n)%tri(1,t) !vertex ID, j
    else; j = o(n)%tri(v+1,t) !vertex ID, j
    endif
    !if (j.le.i) cycle !only draw lines once [this won't always work...]
    ! only cycle if both points are out of view, line routine will cycle pixels
    if (max(o(n)%sp(3,i),o(n)%sp(3,j)).lt.-0.5*scr%FOVx.or. &
      & min(o(n)%sp(3,i),o(n)%sp(3,j)).gt.0.5*scr%FOVx) cycle
    if (max(o(n)%sp(2,i),o(n)%sp(2,j)).lt.0.5*(pi-scr%FOVy).or. &
      & min(o(n)%sp(2,i),o(n)%sp(2,j)).gt.0.5*(scr%FOVy+pi)) cycle
    if (max(o(n)%sp(3,i),o(n)%sp(3,j))-min(o(n)%sp(3,i),o(n)%sp(3,j)).gt.3.14) cycle !check for face wrapping
    if (max(o(n)%sp(2,i),o(n)%sp(2,j))-min(o(n)%sp(2,i),o(n)%sp(2,j)).gt.3.14) cycle !check for face wrapping
    ic(:) = int(o(n)%color(:,i)) !*1.0/max(o(n)%sp(1,i)-10.0,1.0)) !get contrast based on distance
    px1 = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
    ic(:) = int(o(n)%color(:,j)) !*1.0/max(o(n)%sp(1,j)-10.0,1.0)) !get contrast based on distance
    px2 = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
    ! get screen position
    x(1) = int((o(n)%sp(3,i)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
    y(1) = int(scr%sr_y*((o(n)%sp(2,i)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
    x(2) = int((o(n)%sp(3,j)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
    y(2) = int(scr%sr_y*((o(n)%sp(2,j)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
    !if (minval(x).lt.0.or.maxval(x).gt.fb%width)  then; error=.true.; cycle; endif
    !if (minval(y).lt.0.or.maxval(y).gt.fb%height) then; error=.true.; cycle; endif
    call fb%line2c( x(1),y(1), x(2),y(2), px1,px2, scr%sr, o(n)%sp(1,i),o(n)%sp(1,j))
   enddo
  enddo !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
  case("solid") !{{{
  ! print triangles to screen if within Field of View (FOV)
  do t = 1, o(n)%nt
    i = o(n)%tri(1,t) !vertex ID
    j = o(n)%tri(2,t) !vertex ID
    k = o(n)%tri(3,t) !vertex ID
    tp(:,1) = o(n)%sp(:,i)
    tp(:,2) = o(n)%sp(:,j)
    tp(:,3) = o(n)%sp(:,k)
    !if (j.le.i) cycle !only draw lines once [this won't always work...]
    ! only cycle if both points are out of view, line routine will cycle pixels
    if (maxval(tp(3,:)).lt.-0.5*scr%FOVx.or. &
      & minval(tp(3,:)).gt.0.5*scr%FOVx) cycle
    if (maxval(tp(2,:)).lt.0.5*(pi-scr%FOVy).or. &
      & minval(tp(2,:)).gt.0.5*(scr%FOVy+pi)) cycle
    if (maxval(tp(3,:))-minval(tp(3,:)).gt.3.14) cycle !check for face wrapping
    if (maxval(tp(2,:))-minval(tp(2,:)).gt.3.14) cycle 
    ! calculate normal vector for this triangular face (cross product)
    aa(:) = o(n)%rp(:,i) - o(n)%rp(:,j) !vertex i to j
    bb(:) = o(n)%rp(:,i) - o(n)%rp(:,k) !vertex i to k
    norm(1) = aa(2)*bb(3) - aa(3)*bb(2)
    norm(2) = aa(3)*bb(1) - aa(1)*bb(3)
    norm(3) = aa(1)*bb(2) - aa(2)*bb(1)
    dist = sqrt(dot_product(norm,norm)) !make unit vector
    if (dist.lt.1.e-12) cycle ! collapsed triangle (i.e. line) undefined normal
    trinorm(:) = norm(:)/dist !make unit vector
    ! use normal vector and view vector dot product of each vertex to calculate brightness of surface
    if (o(n)%smooth(i).eq.0) then; norm = trinorm
    else; norm = o(n)%rvn(:,i); endif
    aa(:) = o(n)%rp(:,i)/sqrt(dot_product(o(n)%rp(:,i),o(n)%rp(:,i))) !make unit vector
    ic(:) = int(o(n)%color(:,i)*abs(dot_product(aa,norm))) 
    px1 = char(ic(3))//char(ic(2))//char(ic(1))//char(0)

    if (o(n)%smooth(j).eq.0) then; norm = trinorm
    else; norm = o(n)%rvn(:,j); endif
    aa(:) = o(n)%rp(:,j)/sqrt(dot_product(o(n)%rp(:,j),o(n)%rp(:,j))) !make unit vector
    ic(:) = int(o(n)%color(:,j)*abs(dot_product(aa,norm)))
    px2 = char(ic(3))//char(ic(2))//char(ic(1))//char(0)

    if (o(n)%smooth(k).eq.0) then; norm = trinorm
    else; norm = o(n)%rvn(:,k); endif
    aa(:) = o(n)%rp(:,k)/sqrt(dot_product(o(n)%rp(:,k),o(n)%rp(:,k))) !make unit vector
    ic(:) = int(o(n)%color(:,k)*abs(dot_product(aa,norm)))
    px3 = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
    ! FOG: scale vertex colours by the brightness of the surface
!    ic(:) = int(o%color(:,i)/max(tp(1,1)-10.0,1.0)) !get contrast based on distance
!    px1 = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
!    ic(:) = int(o%color(:,j)/max(tp(1,2)-10.0,1.0)) !get contrast based on distance
!    px2 = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
!    ic(:) = int(o%color(:,k)/max(tp(1,3)-10.0,1.0)) !get contrast based on distance
!    px3 = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
    ! get screen position
    x(1) = int((tp(3,1)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
    y(1) = int(scr%sr_y*((tp(2,1)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
    x(2) = int((tp(3,2)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
    y(2) = int(scr%sr_y*((tp(2,2)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
    x(3) = int((tp(3,3)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
    y(3) = int(scr%sr_y*((tp(2,3)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
    !if (minval(x).lt.0.or.maxval(x).gt.fb%width)  then; error=.true.; cycle; endif
    !if (minval(y).lt.0.or.maxval(y).gt.fb%height) then; error=.true.; cycle; endif
    call fb%fillTriangle3c( x(1),y(1), x(2),y(2), x(3),y(3), px1,px2,px3, &
        &  scr%sr, tp(1,1),tp(1,2),tp(1,3) )

    if ( (fb%zbuff(cx,cy).gt.0.0.and.ozb.lt.0.0) .or.   &
       & (fb%zbuff(cx,cy).lt.ozb) ) then 
     ozb = fb%zbuff(cx,cy)
     scr%centerObj = n
     scr%centerTri = t !this triangle ID
     ic(1) = (cx-x(1))*(cx-x(1)) +(cy-y(1))*(cy-y(1))
     ic(2) = (cx-x(2))*(cx-x(2)) +(cy-y(2))*(cy-y(2))
     ic(3) = (cx-x(3))*(cx-x(3)) +(cy-y(3))*(cy-y(3))
     !node of object 'n' with minimum angular distance to center of screen
     ic(1) = minloc(ic,DIM=1)
     scr%centerNode = o(n)%tri(ic(1),t) 
    endif

  enddo !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
 end select


enddo


end subroutine drawObjects !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine commandHelp !{{{
integer :: i, n
character(80), dimension(73) :: c
n = 0
c(:) =""
n=n+1;c(n)="Commands that take parameters show their values when executed sans parms"
n=n+1;c(n)="Command Parameters  Description"
n=n+1;c(n)=" help               Displays this help message"
n=n+1;c(n)=" q,quit,exit        terminates the program"
n=n+1;c(n)=" echo               displays anything else on line to the screen"
n=n+1;c(n)=" exec               executes anything else on line as shell command"
n=n+1;c(n)=" fbpath F           path to framebuffer device, must have write permissions"
!n=n+1;c(n)=" starfile F         path to star file for background stars"
n=n+1;c(n)=" dumpfile F         path for dump file for keyboard screenshots"
n=n+1;c(n)=" tmpdir F           path to a temporary directory for scratch files"
n=n+1;c(n)=" interactive        puts program into interactive rather than single render"
n=n+1;c(n)=" contact N          calculates contact conditions for dynamic objects; N is the"
n=n+1;c(n)="                     type of contact algorithm, 0=off, 1=elastic-Frictionless"
n=n+1;c(n)=" FOV x y            field of view angles (radians)"
n=n+1;c(n)=" subscreen x1 y1 x2 y2 pixel ranges, top-left and bottom right, (row column)"
n=n+1;c(n)=" timestep r         time between frames for animation"
!n=n+1;c(n)=" Nstars n           number of stars to put on background"
n=n+1;c(n)=" dtheta r           angle change for key press"
n=n+1;c(n)=" dspace r           distance increment for key press"
n=n+1;c(n)=" width n            pixel width of screen"
n=n+1;c(n)=" height n           pixel height of screen"
n=n+1;c(n)=" timeout n          animation steps to run"
n=n+1;c(n)=" linelength N       frame buffer line length Bytes"
n=n+1;c(n)=" NobjectBuff N      allocate the object buffer to this many objects. "
n=n+1;c(n)="                     WARNING will delete all existing objects"
n=n+1;c(n)=" LoadObject FILE    will load object from file FILE"
n=n+1;c(n)=" dataColumns i(8)   list the column numbers of a data file that correspond to"
n=n+1;c(n)="           x,y,z,radius,transparency,color(3). color depends on dataColor"
n=n+1;c(n)=" dataRange [xyzrthsvRGB] min max    corresponding to the data in the data file."
n=n+1;c(n)=" dataColorHSV       indicate that color(3) specified in dataColumns are HSV"
n=n+1;c(n)=" dataColorRGB       indicate that color(3) specified in dataColumns are RGB"
n=n+1;c(n)=" dataPoint [point,circ,fcirc,sph,cloud]    object mode type for data points"
n=n+1;c(n)=" LoadData FILE      FILE is a data file with columns, e.g. CSV, SSV."
n=n+1;c(n)="                     loads into the next available object."
n=n+1;c(n)=" saveObject I F     Save object ID=I to file=F. Same as: o I save F"
n=n+1;c(n)=" camera_pos x y z   location of camera" 
n=n+1;c(n)=" camera_orient a u v w  angle-axis orientation of camera" 
n=n+1;c(n)=" universe c x y z   set a universal shape (c=box|sphere), size parameters"
n=n+1;c(n)=" periodic           sets the universal size to be periodic else will rebound"
n=n+1;c(n)=" fbclear            clear the pixel buffer"
n=n+1;c(n)=" dynamics           run one step of object dynamics to update relative vectors"
n=n+1;c(n)=" fbinit             initializes the framebuffer, needed for 'draw'"
n=n+1;c(n)=" bgColor C(3)       draw a filled rectange in subscreen with RGB colour C(3)"
n=n+1;c(n)=" draw               draws current objects to the display buffer (not screen)"
n=n+1;c(n)=" redraw             writes the current pixel buffer to the frame buffer (screen)"
n=n+1;c(n)=" pause              stops rendering animation"
n=n+1;c(n)=" run                resumes rendering animation"
n=n+1;c(n)=" picture [file]     take a screen shot, if 'file' is not given, file='dumpfile'"
n=n+1;c(n)=" record             toggle recording mode (video)"
n=n+1;c(n)=" follow I           camera follows surface of object I"
n=n+1;c(n)=" impulseControl     set input key mode to impulse"
n=n+1;c(n)=" spatialControl     set input key mode to spatial"
n=n+1;c(n)=" cpo I J (K)        copy object ID I to J (through K)"
n=n+1;c(n)=" nearest,closest    reports the nearest vertex relative sph position and obj"
n=n+1;c(n)=" underAim           reports the obj, tri, and node at center of screen"
n=n+1;c(n)=" o I c              object ID 'I' definition 'c'"
n=n+1;c(n)="    object command 'c' can be:"
n=n+1;c(n)="     save s          save this object as file name 's'"
n=n+1;c(n)="     name s          set name of object 'I'"
n=n+1;c(n)="     mode s          set render mode of object 'I', (point|wire|solid)"
n=n+1;c(n)="     mass r          set object's total mass (for collision)"
n=n+1;c(n)="     radius r        set object's average radius (for contact)"
n=n+1;c(n)="     offset x y z    set position offset of object 'I'"
n=n+1;c(n)="     velocity u v w  set velocity of object 'I'"
n=n+1;c(n)="     orient a u v w  set angle and rotation axis of object 'I'"
n=n+1;c(n)="     spin a u v w    set angle rate and rotation axis of object 'I'"
n=n+1;c(n)="     scale r         scale all local positions of object 'I' by multiple 'r'"
n=n+1;c(n)="     add J           add a node to triangle J of object 'I', default J=underAIM"
n=n+1;c(n)="     point J x y z   set local position of point 'J' in object 'I'"
n=n+1;c(n)="     disp J u v w    displace local position of point 'J' in object 'I'"
n=n+1;c(n)="     color J r g b   set color of point 'J' in object 'I' (J=0 for all points)"
n=n+1;c(n)="     smooth J k      set smoothness of point J in object I (J=0 for all points)"
n=n+1;c(n)=" orand c i j  v(:)   randomizes objects i-j property 'c' with ranges v(:)"
n=n+1;c(n)="..Type 'run' to continue interactive animation"
  
do i = 1, n
  write(6,'(A)') trim(c(i))//char(13)
enddo
end subroutine commandHelp !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

end module renderMod
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
