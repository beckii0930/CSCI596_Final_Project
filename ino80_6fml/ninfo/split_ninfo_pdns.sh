fname=1kx5-3lz0_0bp_pdns
ininfo=31
for i in `seq 3 10`
do
 awk '{if($3=='$i'||$1=="<<<<"||$1==">>>>") print $0}' $fname.ninfo > ${fname}_unit${i}.ninfo
 ininfo=$(awk 'BEGIN{print '$ininfo'+1}')
 echo "NINFO(1-2/$i) $ininfo"
 echo "$ininfo = ${fname}_unit${i}.ninfo"
done
