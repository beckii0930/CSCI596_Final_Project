# CSCI596 Final Project: Coarse-Grain MD Simulation for Chromatin Remodeler INO80
#### Yibei Jiang


## I. Introduction

### 1. Chromatin Remodeler IN080 (PDB ID: 6fml)
<p align="center">
  <img src="./ino80.png" width="800">
</p>

<p align="center">
  <b>Figure 1. INO80 Chromatin Remodeler Cryo-EM Structure</b>
</p>

Chromatin Remodelers are protein complexes that can interact directly to consensed chromatin and displace nucleosomes to allow downstream Transcription Factors to access the underlying genomic DNA information. It does so using ATP hydrolysis. The recent influx of remodeler cystal structures shed lights on potential modes of nucleosome sliding but the detailed mechanism is unclear. Here, I investigated the mechanism of chromatin remodeler, INO80 using Coarse Grained Molecular Dynamics Simulations.

```
The Biological Question: 
What is the moleclar mechanism INO80 uses to slide the nucleosome using ATP hydrolysis?
```
## II. Methods
### 1. Coarse Grain Molecular Dynamics
Coarse Grain MD simulations map atoms into virtual particles called beads and, as a result, reduces the number of particles in the system and simplifies calculation. Using CG MD, we can simulate larger system for a longer time. The mapping scheme of each CG situation varies and application specific. Here, I used a 3SPN.2C model by [Hinckley et al.,2013](https://doi.org/10.1063/1.4822042) that is well suited for DNA and AICG2+ model by [Li et al.,2014](https://doi.org/10.1073/pnas.1201807109) for protein.

The general workflow to set up CG MD simulation for the Chromatin Remodeler is described below:
<p align="center">
  <img src="./workflow.png" width="1000">
</p>

<p align="center">
  <b>Figure 2. CG MD Workflow for INO80</b>
</p>

A "**3 Sites per Nucleotide**" strategy is used to coarse grain DNA. 
Below is a sample CG scheme for a nucleotide Cytosine.
<p align="center">
  <img src="./Ccgscheme.png" width="300">
</p>

The potential energy calculations for DNA includes bonded potentials and non-bonded potentials:
<p align="center">
  <img src="./3snp2c.png" width="500">
</p>

<p align="center">
  Equation 1. 3SPN.2C Potential for DNA
</p>
where the detailed calculations are described in [Hinckley et al.,2013](https://doi.org/10.1063/1.4822042).


<p align="center">
  <b>Figure 3. DNA CG Scheme</b>
</p>

A "**1 Site per Amino Acid**" strategy is used to coarse grain proteins. 
Below is a sample CG scheme for a protein with amino acid sequence GLN-GLU-ASP-ASP-ALA.

<p align="center">
  <img src="./pro_cg_scheme.png" width="300">
</p>

<p align="center">
  <b>Figure 4. Protein CG Scheme</b>
</p>

<p align="center">
  <img src="./aicg2+.png" width="500">
</p>

<p align="center">
  Equation 2. AICG2+ Potential for Protein
</p>


### 2. Chromatin Remodeler Slides Nucleosome Using ATP Hydrolysis
Chromatin Remodelers are SNF2 family ATPases. They contain 2 lobes in its ATPase domain. 
<p align="center">
  <img src="./snf2domain.png" width="500">
</p>

<p align="center">
  <b>Figure 5. SNF2 ATPase Domain Structure Has 2 Catalytic Lobes</b>
</p>

Previous CG simulation by [Brandani and Takada,2018](https://doi.org/10.1101/297762) has simulated the ATP Hydrolysis process (apo->ATP->ADP) through changing the **Hydrogen bond strength** and **lobe1-lobe2 Go interaction strength** to facilitate conformational change that allow the Lobe 2 to close up to Lobe 1. The details for Energy Calculation Functions of Go Interaction can be found [here](https://www.charmm.org/wiki//index.php/Coarse_Grained_Go_Models#The_functional_form_of_the_coarse-grained_Go_model). Inspired by their research, I adopted similar strategy to simulate DNA translocation by INO80 using ATP hydrolysis.
<p align="center">
  <img src="./atpaseModel.png" width="300">
</p>

<p align="center">
  <b>Figure 6. ATP Hydrolysis Cycle</b>
</p>

### 3. Define Translocase Lobe GO Contacts
The residues belong to lobe1 and lobe2 were identified using the [Cryo-EM structure 6fml](https://www.nature.com/articles/s41586-018-0029-y). This gives us the atpase in closed conformation. A simple homology model was constructed with [6hts](https://www.rcsb.org/structure/6HTS) to obtain the lobes in open conformation by aligning each of the two lobes with the respective ones from the [6fml](https://www.rcsb.org/structure/6fml) . A list of go interactions between lobes 1 and 2 are then generated (script: "mda_define_atp_contacts.py").
The strength of these contacts in each translocate state is set according to [Brandani and Takada,2018](https://doi.org/10.1101/297762):

|state|Interaction Strength|resulting state    |
|-----|--------------------|-------------------|
|apo  |0.0 kcal/mol        |preferably open    |
|ATP  |0.6 kcal/mol        |preferably closed  |
|ADP  |0.1 kcal/mol        |preferably open    |

**Table 1. Translocase lobes Go Interaction Strength per State**

### 4. Define DNA-Translocase Hydrogen Bonds
The hydrogen bonds are based on the all-atom (charmm27) conformation of the translocase
in complex with nucleosomal DNA as found in PDB 6fml with pdb2gmx. Hydrogen bonds with a cutoff distance of 5Å are identified (script: "mda_define_pdnsSnf2_cut5selmin.py").
The strength of hydrogen bonds in each translocate state is set according to [Brandani and Takada,2018](https://doi.org/10.1101/297762)

|state     |H Bond Strength|resulting state                            |
|----------|---------------|-------------------------------------------|
|apo       |1.8 kBT        |preferably open                            |
|ATP       |1.8 kBT        |preferably closed                          |
|ADP lobe1 |1.8 kBT        |preferably open                            |
|ADP lobe2 |1.44 kBT       |preferably unattach and translocate forward|

**Table 2. DNA-Translocase Hydrogen Bond Strength per State**

The analysis scripts are written in Python v3.6 with MDAnalysis v1.0.0.

### 5. Simulation
First, the CG initial state for INO80 is mimimized using Steepest Descend algorithm for 2X10<sup>7</sup> time steps
Then, I began production run for the DNA sliding activity of INO80 by running a simulation of an ATPase cycle with 3 states (apo-ATP-ADP), one after another. Each state was ran 2X10<sup>5</sup> timesteps with interactions strength modified as described in section 3 an 4 above. This CG MD simulation is performed using [CafeMOl3.0](https://doi.org/10.1021/ct2001045) with equation of motion via Constant temperature Langevin dynamics at 300K.

The last frame of the .dcd trajectory obtained from simulation is used as the intial structure to input file for the next states with changes in parameters in the bond strength files to simulate a whole ATP Hydrolysis cycle. 

## III. Results

First, we obtain the entire Initial CG structure for this system using the CG scheme described in the Methods Section.
![image](https://user-images.githubusercontent.com/25398675/143984154-7b7f0b93-97b7-4076-8595-bdf312867ebc.png)

<p align="center">
  <b>Figure 7. INO80 All Atom to CG Structure</b>
</p>

Next, we started simulation and below is the 2 second snapshots of the 3 state ATP Hydrolysis Cycle. 
In the videos below, 
- ![](https://via.placeholder.com/15/00FFFF/000000?text=+) `Lobe 1` is colored Cyan
- ![](https://via.placeholder.com/15/00FF00/000000?text=+) `Lobe 2` is colored Green
- ![](https://via.placeholder.com/15/FFFF00/000000?text=+) `Insertion Domain` is colored Yellow
- ![](https://via.placeholder.com/15/FFA500/000000?text=+) `DNA` is colored Orange
- ![](https://via.placeholder.com/15/808080/000000?text=+) `Histones` are colored Grey

### 1. The apo state

https://user-images.githubusercontent.com/25398675/149243670-902e0e03-bc89-4086-abda-1d688983e5db.mov

### 2. The ATP state

https://user-images.githubusercontent.com/25398675/149243689-e843434b-e2ba-4787-adf1-7864d02e1c97.mov

### 3. The ADP state

https://user-images.githubusercontent.com/25398675/149243712-20a53ccf-aac9-4f33-95a8-a09eee4d7166.mov

Last, we performed RMSD and RMSF analysis for the trajectory:
#### 1. RMSD: 
RMSD analysis aims to understand the root mean square deviation of the structure compared to the initial structue over time. Below shows the RMSD of the ATPase Lobe1, ATPase Lobe2, Histone, and DNA. 

We see that DNA has the greatest movement due to its 3' end and 5' end linkers. 

The Lobe 2 has greater movements compared to Lobe 1. Lobe 2 has stabilized more in ATP and ADP states compared to apo state. In ATP state, we increased the contact between Lobe 1 and Lobe 2, which can explain the lower RMSD of Lobe2. In ADP state, we lowered contact strength between Lobe 2 and DNA, which is reflected by increased RMSD of Lobe 2. 

Last, the histone movements are similar. 
<p align="center">
  <img src="./RMSD.png" width="140%">
</p>

<p align="center">
  <b>Figure 8. RMSD plot</b>
</p>

#### 2. RMSF: 
RMSF analysis is used to further investigate to pinpoint which parts of the protein or DNA have the greated root mean square fluctuation compare to a computed "middle" structure. 

First, DNA and ATPase RMSF is plotted. The areas close to the domain ends show the greatest flucuation. The RMSF of Lobe2 shows greater fluctuation compared that of Lobe1 as we have seen in RMSD. 

For DNA, the nucleosome DNA dyad position is indicated by a verticle line and we can see the DNA exhibits lots of flucuations in the dyad region and demonstrates a periodic trend. It is well known that the nucleosome DNA has a 10bp periodicity. This flucuation may explain the DNA sliding mechanims by remodeler couples rotation rather than simply sliding on the nucleosome. In the latter case, the RMSF plot will not have such periodicity. 

<p align="center">
  <img src="./RMSF_DNA_ATPase.png" width="120%">
</p>

<p align="center">
  <b>Figure 9. RMSF of 2 Catalytic Lobes and DNA</b>
</p>

Next, histone RMSF is investigated. Solid lines indicate one copy of each histone octamer, whereas the dotted lines indicate the second copy of that. In general, the movements of the 2 copies fluctuates in a similar manner.

<p align="center">
  <img src="./RMSF_Histone_Per_State.png" width="1500">
</p>

<p align="center">
  <b>Figure 10. RMSF of Histone Octamers in Each State</b>
</p>

Last, we investigated histone octamer fluctuations across 3 states: apo, ADP, ATP. H2A fluctuation is the smallest which can point to stabilizing forces of the ATPase with Histone H2A. This can be explained by remodeler ATPase interactions with histone acidic patch formed by H2A-H2B.

<p align="center">
  <img src="./RMSF_Histone_Per_Domain.png" width="800">
</p>

<p align="center">
  <b>Figure 11. RMSF of Each Histone Octamer in 3 States </b>
</p>

## IV. Running Simulations
One will need to download and configure Cafemol software to run this simulation. All files needed to run simulation is in the ino80_6fml folder.

|Folder            |Purpose                                         |
|------------------|------------------------------------------------|
|*.inp             | input files to run each part of the simulation |
|DNA_3spn2c        | DNA parameter files                            |
|hbonds            | Code to genenrate Hydrogen Bond parameter files|
|model             | Code to genenrate Translocase parameter files  |
|ninfo             | All the parameter files needed for simulation  |
|para3_exv.1ssns   | Force Fields for CG models               |
|respac_histones   | Code to genenrate Histone Surface charges      |
|respac_translocase| Code to genenrate Translocase Surface charges  |


**Table 3. Simulation Folder Architecture**

To identify GO contacts:
```
python3 ./ino80_6fml/model/mda_define_atp_contacts.py
```

To identify hydrogen bonds between histone and DNA:
```
python3 ./ino80_6fml/ninfo/mda_define_nsgohb.py
```

To identify hydrogen bonds between translocase and DNA:
```
python3 ./ino80_6fml/hbonds/mda_define_pdnsIno80_cut5selmin.py
```

To run the simulation with input file pt1.inp, type:
```
cafemol ./ino80_6fml/pt1.inp
```
## IV. Conclusion

In this repo, I completed simulations for the INO80 chromatin remodeler and demonstrated the general workflow to simulate ATP Hydrolysis Cycle of this system. The RMSD and RMSF plots are used to examine the simulated trajactory to see the behaviour of the system. In the next phase, we need to perform more analysis, such as free energy profile in order to see whether we have achieved expected INO80 translocation activity.
