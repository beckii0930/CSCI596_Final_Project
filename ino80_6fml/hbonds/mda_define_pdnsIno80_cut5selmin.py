#!/usr/bin/env python

import sys
import numpy as np
from numpy import linalg as LA

import MDAnalysis
import MDAnalysis.analysis.hbonds
# from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

class Hbond:
    pass

class Pdns_entry:
    pass

# define universe
uarray=[]
uarray.append(MDAnalysis.Universe('6fml_extend53.charmm27.pdb'))
# uarray.append(MDAnalysis.Universe('6fml.gro'))
# ridshift=665
# consistent with ninfo/define_groups.sh
ridshift=963
# for u in uarray:
#   for a in u.select_atoms('segid K'): a.residue.resid+=-ridshift  # residue numbering consistent with the cg topology

# analyse =( hbonds
# selO743-ridshift,940-ridshift,1046-ridshift,1212-ridshift)
# selK =(964-ridshift,1273-ridshift,1549-ridshift,1705-ridshift)
selK =(964,1273,1549,1705)
# print(selK)
harray=[]
for u in uarray:
  hbonds=MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u,selection1='(segid K and (resid %d:%d or resid %d:%d))'%(selK),selection2='resname DA DT DC DG and (resid 48:58 or resid 165:175)',acceptors=('O1P','O2P'),distance=5.0,angle=0.0)
  # hbonds=MDAnalysis.analysis.hbonds.HydrogenBondAnalysis(u,selection1='segid K',selection2='resname DA DT DC DG and (resid 48:58 or resid 165:175)',acceptors=('O1P','O2P'),distance=5.0, angle=0.0)
  # print(hbonds)
  harray.append(hbonds)
  # selection1='(segid K and (resid %d:%d or resid %d:%d))'%(selK)
  # # segid K and (resid 1:310 or resid 586:742))
  # selection2='resname DA DT DC DG and (resid 48:58 or resid 165:175)'
  # hbonds = HydrogenBondAnalysis(
  #   universe=u,
  #   acceptors_sel='name O1P O2P',
  #   between=[
  #       [selection1, selection2] # for atpase-DNA hbonds
  #   ],
  #   d_a_cutoff=5.0,d_h_a_angle_cutoff=0.0
  # )

  # harray.append(hbonds)
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
chain2num["K"]=11

# define set
bsets=[]
for h,u in zip(harray,uarray):
  bsets.append(set())
  for i in h.timeseries[0]:
    # print("-----------------------in for loops")
    idonor=i[0] # donor
    iaccep=i[1] # acceptor
    # print(i)
    b=Hbond()
    b.dsegid   =chain2num[u.atoms[idonor].segid]
    b.dresname =u.atoms[idonor].resname
    b.dresid   =u.atoms[idonor].resid-ridshift
    b.asegid   =chain2num[u.atoms[iaccep].segid]
    #b.aresname =u.atoms[iaccep].resname
    b.aresid   =u.atoms[iaccep].resid
    b.true =1
    bsets[-1].add((b.dsegid,b.dresname,b.dresid,b.asegid,b.aresid,b.true))

# sort
for bset in bsets:
  bset=set(sorted(bset, key=lambda b: b[4])) # acceptor_resid
  bset=set(sorted(bset, key=lambda b: b[3])) # acceptor_segid
  bset=set(sorted(bset, key=lambda b: b[2])) # donor_resid
  bset=set(sorted(bset, key=lambda b: b[0])) # donor_segid
# print(bset)
# print("-----------------------building topology")

# build topology from the coarse grained structure
# define universe
uarray=[]
uarray.append(MDAnalysis.Universe('6fml_223bp_1kx5oct_ino80_minim_model2.pdb'))

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
chain2num["K"]=11
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
num2chain[11]="K"
for i in chain2num:
  if(i!=num2chain[chain2num[i]]): print('# error: chain2num not inverse of num2chain'); sys.exit()
ridshift={}
ridshift['K']=len(uarray[0].select_atoms('segid A B C D E F G H I J'))
if(len(set(uarray[0].atoms.segids))!=11): print('# error: nunits!=11'); sys.exit()
for j in range(1,len(uarray)):
  u=uarray[j]
  for i in range(len(u.atoms)):
    if(u.atoms[i].resid!=uarray[0].atoms[i].resid or u.atoms[i].segid!=uarray[0].atoms[i].segid): print('# error: inconsistent cg topologies'); sys.exit()

# print("-----------------------start contact analysis")
# start contact analysis
pdns={} # dictionary of pdns contacts
native=0
for bset,u in zip(bsets,uarray):
  native+=1
  Nbp=(len(u.select_atoms('segid A'))+1)/3
  if(Nbp!=(len(u.select_atoms('segid B'))+1)/3): print('# error: inconsistent Nbp'); sys.exit()
  for b in bset:
    # print("-----------------------standard contact")
    # standard contact
    # print(num2chain[b[0]],b[1],b[2])
    # ca=u.select_atoms('segid %s and resname %s and resid %d'%(num2chain[b[0]],b[1],b[2]))
    ca=u.select_atoms('segid %s and resname %s and resid %d and name CA'%(num2chain[b[0]],b[1],b[2]))
    if(len(ca)!=1): print('# error: CA bead not found'); sys.exit()
    ca=ca[0]
    # print(ca)
    dp=u.select_atoms('segid %s and resname DA DT DC DG and resid %d and name DP'%(num2chain[b[3]],b[4]))
    if(len(dp)!=1): print('# error: P bead not found'); sys.exit()
    dp=dp[0]
    ds =u.atoms[dp.index+1]
    can=u.atoms[ca.index-1]
    cac=u.atoms[ca.index+1]
    if(ds.name!='DS'):  print('# error: name DS !=', ds.name); sys.exit()
    if(can.name!='CA'): print('# error: name CA !=',can.name); sys.exit()
    if(cac.name!='CA'): print('# error: name CA !=',cac.name); sys.exit()
    p2s  =(ds.position -dp.position )/LA.norm(ds.position -dp.position)
    p2ca =(ca.position -dp.position )/LA.norm(ca.position -dp.position)
    cc2cn=(can.position-cac.position)/LA.norm(can.position-cac.position)
    dist=LA.norm( dp.position-ca.position )
    theta1=np.arccos(np.dot(cc2cn,p2ca))/np.pi*180.0
    theta2=np.arccos(np.dot(p2s,  p2ca))/np.pi*180.0
    # add contact to the list
    if(ca.id not in pdns.keys()): pdns[ca.id]=[]
    # print("changed")
    c=Pdns_entry()
    c.chain=chain2num[ca.segid]
    c.id=ca.id
    c.id2=ca.id-ridshift[ca.segid]
    c.dist=dist
    c.theta1=theta1
    c.theta2=theta2
    c.true=b[-1] # was this a true contact, or just the one obtained by symmetry?
    c.native=native
    c.b=b
    pdns[ca.id].append(c)
# sort
for resid in pdns:
  pdns[resid]=sorted(pdns[resid], key=lambda c: c.b[4]) # sort by acceptor_resid

# print out
icontact=0+78 # shift by the number of hydrogen bonds in the histone octamer
print('<<<< pdns')
keylist = pdns.keys()
keylist=sorted(keylist)
if(len(uarray)!=1): print("error: len(uarray)!=1"); sys.exit()
for resid in keylist:
  for i in range(len(pdns[resid])):
    c=pdns[resid][i]
    dist=c.dist
    theta1=c.theta1
    theta2=c.theta2
    # contact+=1; print("pdns",icontact,c.chain,c.id,c.id2,dist,theta1,theta2,"-1.0")
    icontact+=1; print("pdns",icontact,c.chain,c.id,c.id2,dist,theta1,theta2,"-1.0")
  if(len(pdns[resid])>2): print('# error: len(pdns[resid])>2'); sys.exit()
print('>>>>')
