Nbp=223
Ndna=$(awk 'BEGIN{print '$Nbp'*6-2}') # number of dna beads
Noct=980
fname=6fml_chainG_964-1273_1549-1704
awk '{print "CHANGE_CHARGE",$1+1336+980-963,$2}' $fname.inp > ${fname}_nucl223bp_respac.inp
# awk '{print "CHANGE_CHARGE",$1+"$Ndna"+"$Noct",$2}' $fname.inp > ${fname}_nucl223bp_respac.inp
# awk '{print $1,$2+'$Ndna'+'${Noct}',$3}' $fname.inp > ${fname}_nucl223bp_respac.inp
