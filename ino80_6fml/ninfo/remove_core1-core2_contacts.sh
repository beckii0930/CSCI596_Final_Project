awk '
 BEGIN{
  d=1243
  b1=1243-d
  e1=1252-d
  b2=1705-d
  e2=1856-d
 }
 {
  if($1=="contact") {
   if($7>=b1&&$7<=e1&&$8>=b2&&$8<=e2) tmp=0
   else print $0
  }
  else print $0
 }
' 6fml_223bp_ino80_unit11-11.ninfo > 6fml_223bp_ino80_unit11-11_nocore1-2.ninfo
