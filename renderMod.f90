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
 real(4), allocatable, dimension(:,:) :: point,ip,rp,sp
 real(4), allocatable, dimension(:,:) :: vertnorm, ivn, rvn
 integer, allocatable, dimension(:,:) :: color
 integer, allocatable, dimension(:)   :: smooth !average normals or not
 integer, allocatable, dimension(:,:) :: tri
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

! camera object data
cam%name = "Camera"
cam%po  = (/ -10.0, 0.0, 0.0 /) !default camera position
cam%velocity  = (/ 0.0, 0.0, 0.0 /) ! yoru personal velocity
cam%orient = (/ 0.0, 1.0, 1.0, 1.0 /) !initial angle-axis notation
cam%spin = (/ 0.0, 1.0, 1.0, 1.0 /) !initial angle-axis notation
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
contact = .false. !contact detection off by default
cameraContact = .false. !objects don't interact with camera by default
!starfilename = "stars.dat"
!Nstars = 300

record = .false. ! dump PPM files every frame?
scrot = .false.  ! single screen shot flag

open(FD,file=trim(fl),iostat=err)
ln = 0
do 
  line = getLine( 10, ln, err )
!write(6,'(i3,x,A)') ln,trim(line) !echo the read line
  if (err.ne.0) exit 
  if (len(trim(line)).le.0) cycle
  call interpretLine( line, err )
  if (err.gt.0) then
    write(0,*) "ERROR: interpreting line: ",ln," of file:"//trim(fl)
    exit
  endif
enddo
close(FD)
if (err.gt.0) STOP

if (.not.allocated(o)) allocate( o(ObjectBufferSize) )

! convert angle axis notation into actual quaternion rotors
!q1(1) = 0.0; q1(2:4) = scr%camerarotor(2:4) !axis portion
!call quatrotor( q1, scr%camerarotor(1), q2 )
!scr%camerarotor = q2

!q1(1) = 0.0; q1(2:4) = scr%camera_spin(2:4) !axis portion
!call quatrotor( q1, scr%camera_spin(1), q2 )
!scr%camera_spin = q2

q1(1) = 0.0; q1(2:4) = cam%orient(2:4) !axis portion
call quatrotor( q1, cam%orient(1), q2 )
cam%orient = q2

q1(1) = 0.0; q1(2:4) = cam%spin(2:4) !axis portion
call quatrotor( q1, cam%spin(1), q2 )
cam%spin = q2

scr%sr_x = scr%sr(1,2)-scr%sr(1,1) ! x width for 3d screen
scr%sr_y = scr%sr(2,2)-scr%sr(2,1) ! y width


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
implicit none
integer :: wc, error
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
  case ("q","exit","quit"); error=1 ! this will flag as error and kill program
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
  case ("LoadObject"); 
   if (wc==1) then; write(6,*) "Nobjects:",Nobjects; run=.false.; return; endif
   Nobjects = Nobjects + 1
   if (Nobjects.gt.ObjectBufferSize) then
    write(0,*) "ERROR: not loading this object as there is no more space in object buffer"
    error=1; return
   endif
   if (.not.allocated(o)) then
    write(0,*) "ERROR: must allocate object buffer before loading object, see directive: NobjectBuff"
    error=1; return
   endif
   tmp = s_get_word(2, s)
   call loadObject( 12, trim(tmp), o(Nobjects) )
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
  case("periodic"); scr%periodic=.true. !set universe to be periodic
  case("fbclear") ; call fb%clear !clear the fb%pxbuff and fb%zbuff
  case("pause"); run=.false. !stop animation
  case("run"); run=.true. !stop animation
  case("record"); if (record) then; record=.false.; else; record = .true.; endif
  case("picture","scrot"); scrot=.true. !record every frame
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
  !qtmp1(1) = 0.0; qtmp1(2:4) = o(n)%velocity(:)
  !call quatrotate( qtmp1, scr%globalrotor, qtmpv)
!  call quat3rotate( o(n)%velocity, scr%globalrotor, dv)
  dv = o(n)%velocity
  ! integrate position of object in space.
  !o(n)%po(:)   = o(n)%po(:) + qtmpv(2:4)*scr%dt
  o(n)%po(:)   = o(n)%po(:) + dv*scr%dt
  do i = 1, o(n)%np
   ! if the object its self is rotating it can be implemented here
   !  based on its original coordinates.
   !   ip=qlrr*lp*qlrr^-1 + po
   !   rp=qgr*(ip-cam)*qgr^-1
   ! local points rotation to global coordinates
   !qtmp1(1) = 0.0; qtmp1(2:4) = o(n)%point(1:3,i)
   !call quatrotate( qtmp1, o(n)%orient, qtmp2)
   call quat3rotate( o(n)%point(:,i), o(n)%orient, rp)
   o(n)%ip(:,i) = rp +o(n)%po(:)
   
   ! points global coordinates to camera coordinates
   !o(n)%rp(:,i) = o(n)%ip(:,i) - scr%camera_p(:)
   !qtmp1(1) = 0.0; qtmp1(2:4) = o(n)%rp(1:3,i)
   !call quatrotate( qtmp1, scr%globalrotor, qtmp2)
   !o(n)%rp(:,i) = qtmp2(2:4)
   !call quat3rotate( o(n)%ip(:,i)-scr%camera_p, scr%globalrotor, o(n)%rp(:,i))
   call quat3rotate( o(n)%ip(:,i)-cam%po, scr%globalrotor, o(n)%rp(:,i))
   
   call cart2sph(o(n)%rp(:,i), o(n)%sp(:,i))
   ! transform the vertex normals as well.
   if (o(n)%smooth(i).eq.1) then
    ! local points rotation to global coordinates
    !qtmp1(1) = 0.0; qtmp1(2:4) = o(n)%vertnorm(1:3,i)
    !call quatrotate( qtmp1, o(n)%orient, qtmp2)
    !o(n)%ivn(:,i) = qtmp2(2:4) !+o(n)%po(:)
    call quat3rotate( o(n)%vertnorm(:,i), o(n)%orient, o(n)%ivn(:,i))
    
    ! points global coordinates to camera coordinates
    !o(n)%rvn(:,i) = o(n)%ivn(:,i) !- scr%camera_p(:)
    !qtmp1(1) = 0.0; qtmp1(2:4) = o(n)%ivn(1:3,i)
    !call quatrotate( qtmp1, scr%globalrotor, qtmp2)
    !o(n)%rvn(:,i) = qtmp2(2:4)
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
integer :: k, n, i, t, v, j, ic(3), x(3), y(3), cx, cy
character(4) :: px, px1, px2, px3
real(4) :: tp(3,3), aa(3), bb(3), norm(3), trinorm(3), pi, ozb, dist

!pi = 4.0*ATAN2(1.0,1.0) !=log(-1.0)/i
pi = ACOS(-1.0)

scr%centerObj = 0
scr%centerTri = 0
scr%centerNode = 0
!screen center location
cx = scr%sr_x/2 + scr%sr(1,1)
cy = scr%sr_y/2 + scr%sr(2,1)
ozb = fb%zbuff(cx,cy) !old zbuff value at center of screen

do n = 0, Nobjects
 if (o(n)%np.eq.0) cycle !object has no points and thus no triangles to render

 select case(trim(o(n)%mode))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
 case("point") !{{{
  ! print points to screen if within Field of View (FOV)
  do i = 1, o(n)%np
    if (o(n)%sp(3,i).lt.-0.5*scr%FOVx.or.o(n)%sp(3,i).gt.0.5*scr%FOVx) cycle
    if (o(n)%sp(2,i).lt.0.5*(pi-scr%FOVy).or.o(n)%sp(2,i).gt.0.5*(scr%FOVy+pi)) cycle
    !ic(:) = int(o(n)%color(:,i)*1.0/max(o(n)%sp(1,i)-10.0,1.0)) !get contrast based on distance
    ic(:) = int(o(n)%color(:,i)) 
    px = char(ic(3))//char(ic(2))//char(ic(1))//char(0)
    ! get screen position
    x(1) = int((o(n)%sp(3,i)+0.5*scr%FOVx)/scr%FOVx*scr%sr_x) + scr%sr(1,1)
    y(1) = int(scr%sr_y*((o(n)%sp(2,i)-0.5*(pi-scr%FOVy))/scr%FOVy)) + scr%sr(2,1)
    !if (x(1).lt.0.or.x(1).gt.fb%width)  then; error=.true.; cycle; endif
    !if (y(1).lt.0.or.y(1).gt.fb%height) then; error=.true.; cycle; endif
    call fb%putPixel( x(1), y(1), px)
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
character(80), dimension(57) :: c
n = 0
c(:) =""
n=n+1;c(n)="Commands that take parameters show their values when executed sans parms"
n=n+1;c(n)="Command Parameters  Description"
n=n+1;c(n)=" help               Displays this help message"
n=n+1;c(n)=" q,quit,exit        terminates the program"
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
n=n+1;c(n)=" camera_pos x y z   location of camera" 
n=n+1;c(n)=" camera_orient a u v w  angle-axis orientation of camera" 
n=n+1;c(n)=" universe c x y z   set a universal shape (c=box|sphere), size parameters"
n=n+1;c(n)=" periodic           sets the universal size to be periodic else will rebound"
n=n+1;c(n)=" fbclear            clear the pixel buffer"
n=n+1;c(n)=" redraw             writes the current pixel buffer to the frame buffer"
n=n+1;c(n)=" pause              stops rendering animation"
n=n+1;c(n)=" run                resumes rendering animation"
n=n+1;c(n)=" picture            take a screen shot"
n=n+1;c(n)=" record             toggle recording mode (video)"
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
