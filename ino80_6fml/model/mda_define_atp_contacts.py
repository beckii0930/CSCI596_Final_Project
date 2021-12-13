#!/usr/bin/python

import sys
# from sets import Set
import numpy as np
from numpy import linalg as LA

import MDAnalysis
import MDAnalysis.analysis.hbonds

class Hbond:
    pass

class Pdns_entry:
    pass

# define universe
u1=MDAnalysis.Universe('6fml_chainG_resid964-1273_fit.pdb')
u2=MDAnalysis.Universe('6fml_chainG_resid1549-1701_fit.pdb')
ridshift=665

# select CA beads
sel1=u1.select_atoms('(name CA)')
sel2=u2.select_atoms('(name CA)')

# example of ninfo entry
# contact   7485     11     11   2996   3000    680    684      6.4815      1.0000      1      1.6000 p-p
i=10000
unit1=11
unit2=11
for a1 in sel1:
  for a2 in sel2:
    p1=a1.position
    p2=a2.position
    d=LA.norm(p1-p2)
    if(d<7.5):
      i=i+1
      resid1= a1.resid-ridshift
      resid2= a2.resid-ridshift
      serial1=resid1
      serial2=resid2
      print('contact',i,unit1,unit2,serial1,serial2,resid1,resid2,d,'1.0000 1 XXX p-p')

