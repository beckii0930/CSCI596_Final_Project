Nbp=223
Ndna=$(awk 'BEGIN{print '$Nbp'*6-2}') # number of dna beads
fname=1kx5_octamer_model_respac_A44-eB24-eC15-116D36-e
awk '{print $1,$2+'$Ndna',$3}' $fname.inp > ${fname}_shift${Nbp}bp.inp
