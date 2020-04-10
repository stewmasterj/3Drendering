!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
module controlMod !{{{
logical :: rotate, translate, record, scrot, run, impulseControl, contact
logical :: cameraContact
integer :: contactType, Ofollow
end module controlMod !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
module options !{{{
character(80) :: FILobj, FILscreen
end module options !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80
module inputShareMod !{{{
integer :: col(8)   ! column identifiers from data file
real(4) :: rng(2,8) ! ranges for each column
logical :: h        ! flag for HSV
character(40) :: md ! mode 
end module inputShareMod !}}}
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!80

