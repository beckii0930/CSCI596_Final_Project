cat 6fml_223bp_ino80_unit11-11_nocore1-2.ninfo | awk '{if($1=="contact"&&($(NF-1)>0.01))$(NF-1)=1.5; print $0}' > 6fml_223bp_ino80_unit11-11_nocore1-2_strong3.ninfo
#cat 6fml_223bp_ino80_unit11-11_nocore1-2.ninfo | awk '{fac=8.5;if($1=="contact"&&($(NF-1)>0.01)){if(fac*$(NF-1)<1.6) $(NF-1)=fac*$(NF-1); else $(NF-1)=1.6}; print $0}' > 6fml_223bp_ino80_unit11-11_nocore1-2_strong3.ninfo
echo "# original contacts"
cat 6fml_223bp_ino80_unit11-11_nocore1-2.ninfo | awk '{if($1=="contact"&&($(NF-1)>0.01))print $(NF-1)}' | average.sh
echo "# contacts fixed at xxx"
cat 6fml_223bp_ino80_unit11-11_nocore1-2_strong3.ninfo | awk '{if($1=="contact"&&($(NF-1)>0.01))print $(NF-1)}' | average.sh
#echo "# contacts rescaled by 3"
#cat 6fml_223bp_ino80_unit11-11_nocore1-2_strong3.ninfo | awk '{if($1=="contact"&&($(NF-1)>0.01))print $(NF-1)}' | average.sh

