# #!/usr/bin/python

import sys
import numpy as np
from numpy import linalg as LA

import MDAnalysis
import MDAnalysis.analysis.hbonds

class Hbond:
    pass

class Pdns_entry:
    pass

# define universe
uarray=[]
uarray.append(MDAnalysis.Universe('1kx5_model_charmm27.pdb'))
uarray.append(MDAnalysis.Universe('3lz0_adjust_resid_charmm27.pdb'))

# analyse hbonds
# tails edges essentially the same as those used in respac
selH3 =(44,134)
selH4 =(24,101)
selH2A=(16,116) # only exception 16 instead of 15, because of the non-symmetric residue
selH2B=(36,124)
harray=[]
for u in uarray:
  harray.append(MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u,selection1='(segid A and resid %d:%d) or (segid B and resid %d:%d) or (segid C and resid %d:%d) or (segid D and resid %d:%d) or (segid E and resid %d:%d) or (segid F and resid %d:%d) or (segid G and resid %d:%d) or (segid H and resid %d:%d)'%(selH3[0],selH3[1],selH4[0],selH4[1],selH2A[0],selH2A[1],selH2B[0],selH2B[1],selH3[0],selH3[1],selH4[0],selH4[1],selH2A[0],selH2A[1],selH2B[0],selH2B[1]),selection2='resname DA DT DC DG',acceptors=('O1P','O2P')))
  harray[-1].run()

chain2num={} # convert chain consistently with the go model, first dna, then proteins
chain2num["A"]=3
chain2num["B"]=4
chain2num["C"]=5
chain2num["D"]=6
chain2num["E"]=7
chain2num["F"]=8
chain2num["G"]=9
chain2num["H"]=10
chain2num["I"]=1
chain2num["J"]=2

# define set
bsets=[]
for h,u in zip(harray,uarray):
  bsets.append(set())
  for i in h.timeseries[0]:
    idonor=i[0] # donor
    iaccep=i[1] # acceptor
    b=Hbond()
    b.dsegid   =chain2num[u.atoms[idonor].segid]
    b.dresname =u.atoms[idonor].resname
    b.dresid   =u.atoms[idonor].resid
    b.asegid   =chain2num[u.atoms[iaccep].segid]
    #b.aresname =u.atoms[iaccep].resname
    b.aresid   =u.atoms[iaccep].resid
    b.true =1
    bsets[-1].add((b.dsegid,b.dresname,b.dresid,b.asegid,b.aresid,b.true))

# symmetrize contacts
for bset in bsets:
  bsetsymm=set()
  for b in bset:
    dsegid_symm=0
    asegid_symm=0
    if(  b[0]>=3 and b[0]<=6):  dsegid_symm=b[0]+4
    elif(b[0]>=7 and b[0]<=10): dsegid_symm=b[0]-4
    else: print('# error: donor is not part of the protein'); sys.exit()
    if(  b[3]==1): asegid_symm=b[3]+1
    elif(b[3]==2): asegid_symm=b[3]-1
    else: print('# error: acceptor is not part of the dna'); sys.exit()
    bsetsymm.add((dsegid_symm,b[1],b[2],asegid_symm,b[4],1))
  for b in bsetsymm: bset.add(b)

# add contacts to the other native structure
bsetsymms=[set(),set()]
for i in range(len(bsets)):
  bset=bsets[i]
  isymm=0
  for b in bset:
    aresid_symm=0
    if(i==0):
      isymm=1
      if(b[4]>-51 and b[4]<51):
        aresid_symm=b[4]
      elif(b[4]>=51):
        aresid_symm=b[4]-1
      elif(b[4]<=-51):
        aresid_symm=b[4]+1
      else:
        print('# error: symmetrization silly mistake'); sys.exit()
    elif(i==1):
      isymm=0
      if(b[4]>-51 and b[4]<51):
        aresid_symm=b[4]
      elif(b[4]>=51):
        aresid_symm=b[4]+1
      elif(b[4]<=-51):
        aresid_symm=b[4]-1
      else:
        print('# error: symmetrization silly mistake'); sys.exit()
    else: print('# error: isymm'); sys.exit()
    bsetsymms[isymm].add((b[0],b[1],b[2],b[3],aresid_symm,1))
for i in range(len(bsets)):
  bset=bsets[i]
  for b in bsetsymms[i]: bset.add(b)
    
# sort
for bset in bsets:
  bset=set(sorted(bset, key=lambda b: b[4])) # acceptor_resid
  bset=set(sorted(bset, key=lambda b: b[3])) # acceptor_segid
  bset=set(sorted(bset, key=lambda b: b[2])) # donor_resid
  bset=set(sorted(bset, key=lambda b: b[0])) # donor_segid

# buit topology from the coarse grained structure
# define universe
uarray=[]
uarray.append(MDAnalysis.Universe('1kx5_A38_AS147_T38_octamer_model_cg.pdb'))
uarray.append(MDAnalysis.Universe('3lz0_A39_601-3lz0-145bp_T39_octamer_model_cg.pdb'))

# define units and intra-unit resid shift
chain2num={}
chain2num["A"]=1
chain2num["B"]=2
chain2num["C"]=3
chain2num["D"]=4
chain2num["E"]=5
chain2num["F"]=6
chain2num["G"]=7
chain2num["H"]=8
chain2num["I"]=9
chain2num["J"]=10
num2chain={}
num2chain[1]="A"
num2chain[2]="B"
num2chain[3]="C"
num2chain[4]="D"
num2chain[5]="E"
num2chain[6]="F"
num2chain[7]="G"
num2chain[8]="H"
num2chain[9]="I"
num2chain[10]="J"
for i in chain2num:
  if(i!=num2chain[chain2num[i]]): print('# error: chain2num not inverse of num2chain'); sys.exit()
ridshift={}
ridshift['A']=0
ridshift['B']=len(uarray[0].select_atoms('segid A'))
ridshift['C']=len(uarray[0].select_atoms('segid A B'))
ridshift['D']=len(uarray[0].select_atoms('segid A B C'))
ridshift['E']=len(uarray[0].select_atoms('segid A B C D'))
ridshift['F']=len(uarray[0].select_atoms('segid A B C D E'))
ridshift['G']=len(uarray[0].select_atoms('segid A B C D E F'))
ridshift['H']=len(uarray[0].select_atoms('segid A B C D E F G'))
ridshift['I']=len(uarray[0].select_atoms('segid A B C D E F G H'))
ridshift['J']=len(uarray[0].select_atoms('segid A B C D E F G H I'))
if(len(set(uarray[0].atoms.segids))!=10): print('# error: nunits!=10'); sys.exit()
for j in range(1,len(uarray)):
  u=uarray[j]
  for i in range(len(u.atoms)):
    if(u.atoms[i].resid!=uarray[0].atoms[i].resid or u.atoms[i].segid!=uarray[0].atoms[i].segid): print('# error: inconsistent cg topologies'); sys.exit()

# create symmetric contacts on the same structure
bsymm_dic={} # dictionary of symmetric bonds
for bset in bsets:
  bsetsymm=set()
  for b in bset:
    # check symmetric contact
    dsegid_symm=0
    if(  b[0]>=3 and b[0]<=6):  dsegid_symm=b[0]+4
    elif(b[0]>=7 and b[0]<=10): dsegid_symm=b[0]-4
    else: print('# error: donor is not part of the protein'); sys.exit()
    asegid_symm=0
    if(  b[3]==1): asegid_symm=b[3]+1
    elif(b[3]==2): asegid_symm=b[3]-1
    else: print('# error: acceptor is not part of the protein'); sys.exit()
    btmp =(dsegid_symm,b[1],b[2],asegid_symm,b[4],1) # use only for test
    if(btmp not in bset): print('# error: bsymm should have been added before'); print(btmp); sys.exit()
    bsymm_dic[(b[0],b[1],b[2],b[3],b[4])]=(dsegid_symm,b[1],b[2],asegid_symm,b[4])
    bsymm_dic[(dsegid_symm,b[1],b[2],asegid_symm,b[4])]=(b[0],b[1],b[2],b[3],b[4])

# start contact analysis
pdns={} # dictionary of pdns contacts
native=0
for bset,u in zip(bsets,uarray):
  native+=1
  Nbp=(len(u.select_atoms('segid A'))+1)/3
  if(Nbp!=(len(u.select_atoms('segid B'))+1)/3): print('# error: inconsistent Nbp'); sys.exit()
  resid_dyad=(Nbp+1)/2
  for b in bset:
    # standard contact
    ca=u.select_atoms('segid %s and resname %s and resid %d and name CA'%(num2chain[b[0]],b[1],b[2]))
    if(len(ca)!=1): print('# error: CA bead not found'); sys.exit()
    ca=ca[0]
    dp=u.select_atoms('segid %s and resname DA DT DC DG and resid %d and name DP'%(num2chain[b[3]],b[4]+resid_dyad))
    if(len(dp)!=1): print('# error: P bead not found'); sys.exit()
    dp=dp[0]
    ds =u.atoms[dp.index+1]
    can=u.atoms[ca.index-1]
    cac=u.atoms[ca.index+1]
    if(ds.name!='DS'):  print('# error: name DS !='), ds.name; sys.exit()
    if(can.name!='CA'): print('# error: name CA !='),can.name; sys.exit()
    if(cac.name!='CA'): print('# error: name CA !='),cac.name; sys.exit()
    p2s  =(ds.position -dp.position )/LA.norm(ds.position -dp.position)
    p2ca =(ca.position -dp.position )/LA.norm(ca.position -dp.position)
    cc2cn=(can.position-cac.position)/LA.norm(can.position-cac.position)
    dist=LA.norm( dp.position-ca.position )
    theta1=np.arccos(np.dot(cc2cn,p2ca))/np.pi*180.0
    theta2=np.arccos(np.dot(p2s,  p2ca))/np.pi*180.0
    # check the properties of the symmetric contacts
    bsymm=bsymm_dic[(b[0],b[1],b[2],b[3],b[4])]
    casymm=u.select_atoms('segid %s and resname %s and resid %d and name CA'%(num2chain[bsymm[0]],bsymm[1],bsymm[2]))
    if(len(casymm)!=1): print('# error: CA bead not found'); sys.exit()
    casymm=casymm[0]
    dpsymm=u.select_atoms('segid %s and resname DA DT DC DG and resid %d and name DP'%(num2chain[bsymm[3]],bsymm[4]+resid_dyad))
    if(len(dpsymm)!=1): print('# error: P bead not found'); sys.exit()
    dpsymm=dpsymm[0]
    dssymm =u.atoms[dpsymm.index+1]
    cansymm=u.atoms[casymm.index-1]
    cacsymm=u.atoms[casymm.index+1]
    if(dssymm.name!='DS'):  print('# error: name DS !=', dssymm.name); sys.exit()
    if(cansymm.name!='CA'): print('# error: name CA !=',cansymm.name); sys.exit()
    if(cacsymm.name!='CA'): print('# error: name CA !=',cacsymm.name); sys.exit()
    p2ssymm  =(dssymm.position -dpsymm.position )/LA.norm(dssymm.position -dpsymm.position)
    p2casymm =(casymm.position -dpsymm.position )/LA.norm(casymm.position -dpsymm.position)
    cc2cnsymm=(cansymm.position-cacsymm.position)/LA.norm(cansymm.position-cacsymm.position)
    distsymm=LA.norm( dpsymm.position-casymm.position )
    theta1symm=np.arccos(np.dot(cc2cnsymm,p2casymm))/np.pi*180.0
    theta2symm=np.arccos(np.dot(p2ssymm,  p2casymm))/np.pi*180.0
    # check difference
    if(np.fabs(dist-distsymm)>1.0 or np.fabs(theta1-theta1symm)>11.0 or np.fabs(theta2-theta2symm)>12.5):
      print('# error: non-symm contact'); print(b,dist,distsymm,theta1,theta1symm,theta2,theta2symm); sys.exit()
    # symmetrize pdns contact
    dist=(dist+distsymm)/2.0
    theta1=(theta1+theta1symm)/2.0
    theta2=(theta2+theta2symm)/2.0
    if(np.fabs(dist-distsymm)>0.5): print('# error: incorrect dist average'); sys.exit()
    if(np.fabs(theta1-theta1symm)>5.5): print('# error: incorrect theta1 average'); sys.exit()
    if(np.fabs(theta2-theta2symm)>6.25): print('# error: incorrect theta2 average'); sys.exit()
    # add contact to the list
    if(ca.id not in pdns.keys()): pdns[ca.id]=[]
    c=Pdns_entry()
    c.chain=chain2num[ca.segid]
    c.id=ca.id
    c.id2=ca.id-ridshift[ca.segid]
    c.dist=dist
    c.theta1=theta1
    c.theta2=theta2
    c.true=b[-1] # was this a true contact, or just the one obtained by symmetry?
    c.native=native # 1 for 1kx5 and 2 for 601
    c.b=b
    pdns[ca.id].append(c)
# sort
for index in pdns:
  pdns[index]=sorted(pdns[index], key=lambda c: c.b[4]) # sort by acceptor_resid

icontact=0
print('<<<< pdns')
keylist = sorted(pdns)
for index in keylist:
  if(len(pdns[index])%2==0):
    # average contacts from different native structures, they should be sorted in the proper way
    for i in range(len(pdns[index])//2):
      c1=pdns[index][i*2]
      c2=pdns[index][i*2+1]
      if(c1.native==c2.native): # error, you should average different native structures
        print('# error: len(pdns[index])==2 and c1.native==c2.native'); sys.exit()
      if(np.fabs(c1.dist-c2.dist)>1.5 or np.fabs(c1.theta1-c2.theta1)>15.0 or np.fabs(c1.theta2-c2.theta2)>27.0):
        print('# error: dist_1kx5-dist_601>1.5 A or dtheta1>15 deg or dtheta2>27')
        print(c1.b,'native=',c1.native,c1.dist,c1.theta1,c1.theta2)
        print(c2.b,'native=',c2.native,c2.dist,c2.theta1,c2.theta2)
        sys.exit()
      c=pdns[index][0]
      dist=(c1.dist+c2.dist)/2.0
      theta1=(c1.theta1+c2.theta1)/2.0
      theta2=(c1.theta2+c2.theta2)/2.0
      if(np.fabs(dist-c2.dist)>0.75 or np.fabs(dist-c1.dist)>0.75): print('# error: incorrect dist average'); sys.exit()
      if(np.fabs(theta1-c2.theta1)>7.5  or np.fabs(theta1-c1.theta1)>7.5 ): print('# error: incorrect theta1 average'); sys.exit()
      if(np.fabs(theta2-c2.theta2)>13.5 or np.fabs(theta2-c1.theta2)>13.5): print('# error: incorrect theta2 average'); sys.exit()
      icontact+=1; print("pdns",icontact,c.chain,c.id2,c.id2,dist,theta1,theta2,"-1.0")
  if(len(pdns[index])%2!=0 or len(pdns[index])>4): print('# error: len(pdns[index])%2!=0 or len(pdns[index])>4'); sys.exit()
print('>>>>')
