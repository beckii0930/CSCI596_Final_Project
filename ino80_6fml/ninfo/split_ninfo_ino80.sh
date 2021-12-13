nunit=11
ininfo=0
fname=6fml_223bp_ino80
for i in `seq 1 $nunit`
do
 for j in `seq $i $nunit`
 do
  awk '{if($3=='$i'&&$4=='$j') print $0}' $fname.ninfo > ${fname}_unit${i}-${j}.ninfo 
  nlines=`wc -l ${fname}_unit${i}-${j}.ninfo | awk '{print $1}'`
  if [ $nlines -eq 0 ]
  then rm ${fname}_unit${i}-${j}.ninfo
  else ininfo=$(awk 'BEGIN{print '$ininfo'+1}') ; echo "NINFO(${i}/${j}) $ininfo" ; echo "$ininfo = ${fname}_unit${i}-${j}.ninfo"
  fi
 done
done
