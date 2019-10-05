! vim:fdm=marker
!         type 'za' to refold the folds
module ArgsMod

!============================================================================80
!
! Subroutines for easy commandline parsing
!
! getOpt( switch, n, value, c80err )
!       Returns the value array of n elements that are provided to the
!        command line switch. Any error is returned in the character
!        array c80err.
!
! Example:
!       For a given commandline set of arguments:
!             program  -a 3.1  -p T F T  -s Name  -n 30
!       
!       call getOpt( "-s", 1, name, cerr )
!       call getOpt( "-n", 1, num, cerr )
!       call getOpt( "-a", 1, dist, cerr )
!       call getOpt( "-p", 3, PBC, cerr )
!
! Return values: c80err
!       c80err="-1: Init"
!       c80err="0: Successful"
!       c80err="1: No Value for given switch"
!       c80err="2: Value: "//trim(opt2)//" is not an integer number"
!       c80err="3: This switch requires "//cn//" values"
!       c80err="4: Switch "//switch//" not found among arguments"
!
! Be sure to define default values in case a switch is not specified on the
!  command line.
!
! The above subroutine is overloaded with the following four routines.
! This allows the programmer to supply either a Character, integer, real(4)
!  or logical variable as the returning value.
!
! getCOpt( switch, n, value, c80err )
! getIOpt( switch, n, value, c80err )
! getR4Opt( switch, n, value, c80err )
! getLOpt( switch, n, value, c80err )
!
!
! Author: Ross J. Stewart
! Date: Thursday, May 17, 2012           
! email: rossjohnsonstewart@gmail.com
!
!============================================================================80

! example of function overloading
interface getOpt
 module procedure getCOpt, getIOpt, getLOpt, getR4Opt
end interface

contains
subroutine getCOpt(switch,n,value,c80err) !{{{
implicit none
integer,               intent(in)  :: n
integer                            :: i, j
character(len=2),      intent(in)  :: switch
character(len=2)                   :: opt, cn
character(len=80),     intent(out) :: c80err
character(len=80)                  :: opt2
character(len=80), dimension(n), intent(out) :: value

value=" "
c80err="-1: Init"

do i=1,iargc()
 call getarg(i,opt)
 if ( opt .eq. switch ) then
   ! if this is the switch read next for a value
   do j=1,n
     if ((i+j).gt.iargc()) then
      write(cn,'(I2)') n
      c80err="3: This switch requires "//cn//" values"
      Return
     endif
     call getarg(i+j,opt2)
     if ( opt2(1:1) .eq. '-' ) then
     ! this is a switch too!
       opt2=" " ! this option doesn't have or take a value
       c80err="1: No Value for given switch"
       Return
     endif
     value(j)=trim(opt2)
     c80err="0: Successful"
   enddo
  if (n.eq.0) c80err="0: Successful"
 endif
 if (i.eq.iargc().and.c80err(1:1).ne."0") then
  c80err="4: Switch "//switch//" not found among arguments"
  Return
 elseif (i.eq.iargc()) then
  c80err="0: Successful"
  Return
 endif ! end, not a switch
enddo

end subroutine !}}}
subroutine getIOpt(switch,n,value ,c80err) !{{{
implicit none
integer,               intent(in)  :: n
integer                            :: i, j, k
character(len=2),      intent(in)  :: switch
character(len=2)                   :: opt, cn
character(len=80),     intent(out) :: c80err
character(len=80)                  :: opt2
character(len=14)                  :: integerNumbers
integer, dimension(n), intent(out) :: value
logical                            :: good

good=.true.
integerNumbers="-+0123456789 "//char(09)
value=0
c80err="-1: Init"

do i=1,iargc()
 call getarg(i,opt)
 if ( opt .eq. switch ) then
   ! if this is the switch read next for a value
   do j=1,n
     if ((i+j).gt.iargc()) then
      write(cn,'(I2)') n
      c80err="3: This switch requires "//cn//" values"
      Return
     endif
     call getarg(i+j,opt2)
     if ( opt2(1:1) .eq. '-' ) then
     ! this is a switch too!
       opt2=" " ! this option doesn't have or take a value
       c80err="1: No Value for given switch"
       Return
     endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
     do k=1,80
      !Check if this is a real value
      good = good .and. (index(integerNumbers,opt2(k:k)).ne.0)
      if (opt2(k:k).eq." ".or.opt2(k:k).eq.char(09)) exit !end of string
     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
     if (good) then
     ! if it's a number write it as a real
      read(opt2,*) value(j)
      c80err="0: Successful"
     else
      c80err="2: Value: "//trim(opt2)//" is not an integer number"
      Return
     endif
   enddo
  if (n.eq.0) c80err="0: Successful"
 endif
 if (i.eq.iargc().and.c80err(1:1).ne."0") then
  c80err="4: Switch "//switch//" not found among arguments"
  Return
 elseif (i.eq.iargc()) then
  c80err="0: Successful"
  Return
 endif ! end, not a switch
enddo

end subroutine !}}}
subroutine getR4Opt(switch,n,value,c80err) !{{{
implicit none
integer,               intent(in)  :: n
integer                            :: i, j, k
character(len=2),      intent(in)  :: switch
character(len=2)                   :: opt, cn
character(len=80),     intent(out) :: c80err
character(len=80)                  :: opt2
character(len=19)                  :: realNumbers
real(4), dimension(n), intent(out) :: value
logical                            :: good

good=.true.
realNumbers="-+0123456789.END "//char(09)
value=0.0
c80err="-1: init"

do i=1,iargc()
 call getarg(i,opt)
 if ( opt .eq. switch ) then
   ! if this is the switch read next for a value
   do j=1,n
     if ((i+j).gt.iargc()) then
      write(cn,'(I2)') n
      c80err="3: This switch requires "//cn//" values"
      Return
     endif
     call getarg(i+j,opt2)
   !  if ( opt2(1:1) .eq. '-' ) then
   !  ! this is a switch too!
   !    opt2=" " ! this option doesn't have or take a value
   !    c80err="1: No Value for given switch"
   !    Return
   !  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
     do k=1,80
      !Check if this is a real value
      good = good .and. (index(realNumbers,opt2(k:k)).ne.0)
      if (opt2(k:k).eq." ".or.opt2(k:k).eq.char(09)) exit !end of string
     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
     if (good) then
     ! if it's a number write it as a real
      read(opt2,*) value(j)
      c80err="0: Successful"
     else
      c80err="2: Value: "//trim(opt2)//" is not a real number"
      Return
     endif
   enddo
  if (n.eq.0) c80err="0: Successful"
 endif
 if (i.eq.iargc().and.c80err(1:1).ne."0") then
  c80err="4: Switch "//switch//" not found among arguments"
  Return
 elseif (i.eq.iargc()) then
  c80err="0: Successful"
  Return
 endif ! end, not a switch
enddo

end subroutine !}}}
subroutine getLOpt(switch,n,value,c80err) !{{{
implicit none
integer,               intent(in)  :: n
integer                            :: i, j
character(len=2),      intent(in)  :: switch
character(len=2)                   :: opt, cn
character(len=80),     intent(out) :: c80err
character(len=80)                  :: opt2
logical, dimension(n), intent(out) :: value

value=.false.
c80err="-1: Init"
do i=1,iargc()
 call getarg(i,opt)
 if ( opt .eq. switch ) then
   ! if this is the switch read next for a value
   do j=1,n
     if ((i+j).gt.iargc()) then
      write(cn,'(I2)') n
      c80err="3: This switch requires "//cn//" values"
      Return
     endif
     call getarg(i+j,opt2)
     if ( opt2(1:1) .eq. '-' ) then
     ! this is a switch too!
       opt2=" " ! this option doesn't have or take a value
       c80err="1: No Value for given switch"
       Return
     endif
     read(opt2,*) value(j)
     c80err="0: Successful"
   enddo
  if (n.eq.0) c80err="0: Successful"
 endif
 if (i.eq.iargc().and.c80err(1:1).ne."0") then
  c80err="4: Switch "//switch//" not found among arguments"
  Return
 elseif (i.eq.iargc()) then
  c80err="0: Successful"
 Return
 endif ! end, not a switch
enddo

end subroutine !}}}

end module ArgsMod
