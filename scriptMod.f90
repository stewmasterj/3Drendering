! vim: fdm=marker
! Module part of the LineParseMod folder
! DATE: Wed Jun  3 20:15:43 EDT 2020
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
module scriptMod
implicit none
integer :: Nscript, Nvars, Nlabels
character(80), dimension(:), allocatable :: script
!character(40), dimension(:), allocatable :: vnam ! variable name
!character(40), dimension(:), allocatable :: vval ! variable value
character(40), dimension(:), allocatable :: labels  ! goto labels
integer, dimension(:), allocatable :: snum  ! script lines for labels

type scriptvariable
 character(40) :: n ! variable name
 integer :: nv ! number of array elements
 character(80), dimension(:), allocatable :: v ! value
end type scriptvariable
type(scriptvariable), dimension(:), allocatable :: svar

! call loadScript( FD, FN )
! call runScript

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine loadScript( FD, FN ) !{{{
use lineParse
implicit none
integer, intent(in) :: FD
character(*), intent(in) :: FN
integer :: err, nn, nl, wc, i, j, ns
character(80) :: line, word

open(FD,file=trim(FN),iostat=err)
if (err.ne.0) then
  write(0,*) "ERROR: reading script file:"//trim(FN)
  write(0,*) "ERROR: error number: ", err
  STOP
endif

nn = 0; nl = 0
do
  read(FD,'(A)',end=200) line
  call left_of( "#", line )
  wc = s_word_count(line)
  if (wc.gt.0) then
  ! if not empty line
    nn = nn + 1
    nn = nn + s_pat_count(";",line) ! semicolon denotes new command line
    word = s_get_word(1,line)
    if (word.eq."LBL".or.word.eq."lbl") nl = nl + 1 ! must be first word
  endif
enddo
200 continue

Nscript = nn
Nlabels = nl
!write(0,*) "read noncommented lines from script file: ", Nscript

allocate( script(Nscript), labels(nl), snum(nl) )
rewind(FD)

nn = 0; nl =0
do
  read(FD,'(A)',end=300) line
  call left_of( "#", line )
  wc = s_word_count(line)
  if (wc.gt.0) then
  ! if not empty line
    ns = s_pat_count(";",line) ! semicolon denotes new command line
    if (ns.eq.0) then
     nn = nn + 1
     script(nn) = trim(adjustl(line)) ! save this line of the script
    else
     j = 0
     do i = 0, ns
      nn = nn + 1
      word = line(j+1:80) ! save this line of the script
      j = index(trim(word(j+1:80)),";") ! when does ; appear in word?
      call left_of(";",word)
      script(nn) = trim(word)
     enddo
    endif
    word = s_get_word(1,line)
    if (word.eq."LBL".or.word.eq."lbl") then
     nl = nl + 1
     if (s_word_count(script(nn-ns)).ne.2) then
      write(0,*) "ERROR: sline:",nn-ns,"bad command: "//trim(script(nn-ns))
      write(0,*) "ERROR: the LBL command requires an argument:"
      write(0,*) " Usage:  LBL LABEL"
      write(0,*) "   LABEL  can be an arithmetic expression or string."
      STOP
     endif 
     labels(nl) = trim(s_get_word(2,line))   ! save label name (counts duplicats so beware)
     snum(nl)   = nn-ns          ! save script line
!     write(0,*) "LBL "//trim(labels(nl)),nn
    endif
  endif
enddo
300 continue

close(FD)

!do i = 1, nn
! write(0,'(A)') trim(script(i))
!enddo
!call flush(0)

end subroutine loadScript !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
subroutine runScript
use lineParse
use renderMod
implicit none
integer :: i, j, ai, wc, n, lastcall, err
character(80) :: line, w(20), tmp, tmp2
real(4) :: rv ! real value from evaluate
logical :: debug

!allocate( vnam(100), vval(100) ) ! users prolly don't need more than 100 vars
allocate( svar(100) )
!vnam(:) = ""
!vval(:) = ""
svar(:)%n = ""
svar(:)%nv = 0 ! zero means unallocated
Nvars = 0
debug = .false.
lastcall = 1 ! script line number (snum) of last called GOTO with return feature

i = 0
do ! infinite loop cuz the script might never end!
 i = i + 1
 if (i.gt.Nscript) then ! out of bounds or end of script?
   exit
 endif
 line = script(i)
 !!!! Parse line for variables and expressions before commands !!!!{{{
!write(0,*) "s:"//trim(line)
 ! need to substitute variables if they exist
 do while (index(line,"[$").ne.0)  ! while there are brakets with variables
  w(1) = s_get_between("[",line)  ! get the expression
! if (debug) write(0,*) i,"a:"//trim(w(1))
  w(2) = subVars( w(1) )          ! substitute the variables
! if (debug) write(0,*) i,"b:"//trim(w(2))
  line = s_sub( w(1), line, w(2) ) ! redefine the line
! if (debug) write(0,*) i,"c:"//trim(line)
 enddo
 w(1) = subVars( line )  ! sub unnested variables
 line = w(1)
!write(0,*) "v:"//trim(line)
 wc = s_word_count(line)
!write(0,*) "wc:",wc
 ! evaluate the "words" for arithmetic expressions
 do j = 1, wc
   w(j) = s_get_word(j,line)
!write(0,*) j,"w:"//trim(w(j))
   if (j .eq. 1) cycle ! no expressions as first word, only commands
   !write(w(j+1),*) evaluate( w(j) )  ! this makes a number
   w(j+1) = evaluate( w(j) )  
   w(j+2) = s_sub( w(j), line, w(j+1) ) ! substitute the evaluation 
   !w(j+2) = s_sub( w(j), line, evaluate(w(j)) ) ! substitute the evaluation 
   w(j) = trim(w(j+1))  !using w(j+1) as temp variable
   line = trim(w(j+2))  !using w(j+2) as temp variable
 enddo
 if (debug) write(0,*) i,"r:"//trim(line) !}}}
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!60
 ! Let's parse the commands !!! {{{
 select case (trim(w(1)))
  case ("DEBUG"); debug = .true.
  case ("PRINT"); tmp = adjustl(line)
   write(6,'(A)') trim(tmp(6:80))
  case ("EXEC"); tmp = adjustl(line)
   call system(trim(tmp(6:80)))
  case ("LBL"); ! do nothing
  case ("GOTO")
   if (wc.ne.2.and.wc.ne.3) then
    write(0,*) "ERROR: sline:",i,"bad command: "//trim(line)
    write(0,*) "ERROR: the GOTO command requires an argument:"
    write(0,*) " Usage:  GOTO  LABEL [RETURN]"
    write(0,*) "   LABEL  can be an arithmetic expression or string that evaluates to an"
    write(0,*) "          existing label."
    write(0,*) "   [RETURN] optional keyword will return to this line when RETURN is called"
    STOP
   endif 
   do j = 1, Nlabels  ! find the label
    if (trim(w(2)).eq.trim(labels(j))) then
     if (wc.ge.3) then
       if (trim(w(3)).eq."RETURN") lastcall = i ! save this line incase of a RETURN
     endif
     i = snum(j)
!write(0,*) "jumping to sline:",i,"label:"//trim(labels(j))
     exit  ! found. change script line to line after label
    endif
   enddo
  case ("IFGO"); 
   if (wc.ne.3.and.wc.ne.4) then
    write(0,*) "ERROR: sline:",i,"bad command: "//trim(line)
    write(0,*) "ERROR: the IFGO command requires two arguments:"
    write(0,*) " Usage:  IFGO  EXPR  LABEL [RETURN]"
    write(0,*) "   EXPR   is an arithmetic expression, if (EXPR >= 1.0) is true"
    write(0,*) "   LABEL  a label for a target GOTO"
    write(0,*) "   [RETURN] optional keyword will return to this line when RETURN is called"
    STOP
   endif 
   read(w(2),*) rv ! read condition
   if (rv.ge.1.0) then ! condition true
   do j = 1, Nlabels  ! find the label
    if (trim(w(3)).eq.trim(labels(j))) then
     if (wc.ge.4) then
       if (trim(w(4)).eq."RETURN") lastcall = i ! save this line incase of a RETURN
     endif
     i = snum(j)  ! found. change script line to line after label
!write(0,*) "jumping to sline:",i,"label:"//trim(labels(j))
     exit
    endif
   enddo
   endif
  case ("RETURN"); i = lastcall+1 ! actually start on the line after the call 
!!!!##########################################################################80
!!!! Add your user defined commands here !!!!
!  case ("mycommand")
!   read(w(2),*) myargument
!!!!##########################################################################80
  case default !}}}
!!! USER ROUTINE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
  call interpretLine( line, err ) ! err=-1 user called exit, -2 command not found
  if (err.eq.-1) exit ! user defined exit
  if (err.gt.0) then
    write(0,*) "Fatal ERROR: sline:",i,":"//trim(line)
    exit
  endif
  if (err.eq.0) cycle ! successful read of interpret line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!40
!!!! Set variables and resolve variable functions !!!!{{{
   if (wc.gt.1) then
    select case (trim(w(2))) ! check for variable functions
     case ("WORD")
      if (wc.lt.4) then
       write(0,*) "ERROR: in WORD: expecting more arguments, sline",i
       STOP
      endif
       j = index(trim(line),trim(w(3))) + len(trim(w(3)))
       if (index(w(3),".").ne.0) then
         read(w(3),*) rv; ai=nint(rv) ! make integer 
       else
         read(w(3),*) ai ! it's an integer
       endif
       w(2) = trim( s_get_word(ai,line(j:80)) )
     case ("EQUAL")
      if (wc.lt.4) then
       write(0,*) "ERROR: in EQUAL: expecting more arguments, sline",i
       STOP
      endif
      if (trim(w(3)).eq.trim(w(4))) then
       w(2) = "1.0"
      else
       w(2) = "0.0"
      endif
     case ("WC")
      write(w(2),*) wc-2 ! word count in remaining line
!!!!##########################################################################80
!!!! Add your user defined variable functions here !!!!
!   case ("mycommand")
!    w(2) = callfunc( w(3) )
!!!!##########################################################################80
     case default; ! override w(2) with fill line
     ! check for a possibel = sign, don't want to include that
     j = max(index(trim(line),trim(w(1)))+len(trim(w(1))) , index(trim(line),"=")+1)
     w(2) = line(j:80)
     !w(2) = s_get_word_range( 2, wc, line ) ! this removes too much stuff
    end select
    !! set the result, if any, to the variable value
    if (index(trim(line),"=").gt.0) then ! there's an equal sign in here
     tmp = w(1)
     ai = 1 ! default to a single array index
     if (index(w(1),"[").ne.0) then ! there's an array index
       call left_of("[",tmp) ! strip it off nice and slow
       tmp2 = s_get_between("[",w(1)) ! get array index
       if (index(tmp2,".").ne.0) then
         read(tmp2,*) rv; ai=nint(rv) ! integer of array index
       else
         read(tmp2,*) ai ! it's an integer
       endif
     endif
     ! look to see if it exists, if not make it
     n = 0
     do j = 1, Nvars ! loop over existing variable names
      !if (trim(w(1)).eq.trim(svar(j)%n)) then ! found it
      if (trim(tmp).eq.trim(svar(j)%n)) then ! found it
       n = j
      endif
     enddo
     if (n .eq. 0) then ! wasn't found, then make it
      n = Nvars + 1
      Nvars = Nvars + 1
      svar(n)%n = trim(tmp) ! set the variable name
      allocate( svar(n)%v(ai) ) ! if you don't use the ARRAY command
                                ! you could just: x[10] = 0, to make it.
      svar(n)%nv = ai
     endif
     if (.not.allocated( svar(n)%v ).or.ai.gt.svar(n)%nv) then
      write(0,*) "ERROR: trying to set variable value to out of bounds of " &
           & //"variable array:"//trim(svar(n)%n)//" (",ai," of ",svar(n)%nv, &
           & ") on sline:",i
      STOP
     endif
     svar(n)%v(ai) = trim(w(2))  ! set the word to the variable value
  !write(0,*) "variable:"//trim(vnam(n))//" set to:"//trim(vval(n))
     cycle
    endif
   endif
   ! if it's not an equation then...
   write(0,*) "WARNING: sline:",i,"unknown command line: "//trim(line)
 end select !}}}
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!60
enddo

call flush(0); call flush(6)
if (err.ne.0) STOP "done reading script file"

end subroutine runScript
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
function subVars( s ) !{{{
use lineParse
implicit none
character(80) :: s, subVars, vv, tmp
character(40) :: delims
integer :: nv, n, i, ai, k, endpos, startpos, ls
real(4) :: rv

delims = '" +-=/^$~`!@#%*(){}[];:<>,.?|\'
ls = len(trim(s))
nv = 0 ! number of variables on line
do i = 1, ls
 if (s(i:i).eq."$") then
  nv = nv + 1
 endif
enddo

subVars = s
do n = 1, nv
  ! this only makes one substitution
!write(0,*) "subVars",n, trim(subVars)
  startpos = 0
  endpos = 0
  ls = len(trim(subVars))
  do i = 1, ls
   if (subVars(i:i).eq."$") then
    startpos = i
    exit
   endif
  enddo
 
  do i = startpos+1, ls
    if (index(delims,subVars(i:i)).ne.0) then
      endpos = i-1
      ! found variable position, now substitute it
      exit ! done with this variable
    elseif (i.eq.ls) then
      endpos = i
      exit
    endif
  enddo ! substitution of var

  k = 0
  do i = 1, Nvars
   if (subVars(startpos+1:endpos).eq.trim(svar(i)%n)) then
    k = i  ! found the variable name
    exit
   endif
  enddo
  if (k.eq.0) then
    write(0,*) "ERROR: variable:"//subVars(startpos:endpos)//" not found"
    STOP
  endif

  ! need to check for possible array index
  ai = 1
  if (endPos.lt.ls-3) then ! maybe?
    if (subVars(endpos+1:endpos+1) .eq. "[") then ! yup
      tmp = s_get_between( "[", subVars(endpos+1:ls) )
      endpos = min(endpos+len(tmp),ls)
       if (index(tmp,".").ne.0) then
         read(tmp,*) rv; ai=nint(rv) ! integer of array index
       else
         read(tmp,*) ai ! it's an integer
       endif
    endif
  endif

  vv = s_sub( subvars(startpos:endpos), trim(subVars), svar(k)%v(ai) )
  subvars = vv
enddo ! n vars

end function subVars !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
end module scriptMod
