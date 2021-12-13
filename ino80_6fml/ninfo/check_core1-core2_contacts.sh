awk '
 /contact/{
  d=1243
  b1=1243-d
  e1=1252-d
  b2=1705-d
  e2=1856-d
  if($7>=b1&&$7<=e1&&$8>=b2&&$8<=e2) print $0
 }
' 6fml_223bp_ino80_unit11-11.ninfo
