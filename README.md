# CSCI596 Final Project: Coarse-Grain MD Simulation for Chromatin Remodeler INO80
#### Yibei Jiang


## I. Introduction

### 1. Chromatin Remodeler IN080 (PDB ID: 6fml)
<p align="center">
  <img src="./ino80.png" width="800">
</p>
**Figure 1. INO80 Chromatin Remodeler Cryo-EM Structure**

Chromatin Remodelers are protein complexes that can interact directly to consensed chromatin and displace nucleosomes to allow downstream Transcription Factors to access the underlying genomic DNA information. It does so using ATP hydrolysis. The recent influx of remodeler cystal structures shed lights on potential modes of nucleosome sliding but the detailed mechanism is unclear. Here, I investigated the mechanism of chromatin remodeler, INO80 using Coarse Grained Molecular Dynamics Simulations.

```
The Biological Question: 
What is the moleclar mechanism INO80 uses to slide the nucleosome using ATP hydrolysis?
```
## II. Methods
### 0. Coarse Grain Molecular Dynamics
Coarse Grain MD simulations map atoms into virtual particles called beads and, as a result, reduces the number of particles in the system and simplifies calculation. Using CG MD, we can simulate larger system for a longer time. The mapping scheme of each CG situation varies and application specific. Here, I used a 3SPN.2C model by [Hinckley et al.,2013](https://doi.org/10.1063/1.4822042) that is best suited for DNA and AICG2+ model by [Jones et al.](https://doi.org/10.1021/jacs.1c05219) for protein.

The general workflow to set up CG MD simulation for the Chromatin Remodeler is described below:
<p align="center">
  <img src="./workflow.png" width="1000">
</p>
**Figure 2. CG MD Workflow for INO80**


A "**3 Sites per Nucleotide**" strategy is used to coarse grain DNA. Below is a sample CG scheme for a nucleotide C.
<p align="center">
  <img src="./Ccgscheme.png" width="300">
</p>
*Figure 3. DNA CG Scheme*

A "**1 Site per Amino Acid**" strategy is used to coarse grain protein. Below is a sample CG scheme for a protein
sequence GLN-GLU-ASP-ASP-ALA.
<p align="center">
  <img src="./pro_cg_scheme.png" width="300">
</p>
**Figure 4. Protein CG Scheme**

Now we obtain the entire Initial CG structure for this system.
![image](https://user-images.githubusercontent.com/25398675/143984154-7b7f0b93-97b7-4076-8595-bdf312867ebc.png)
**Figure 5. INO80 All Atom to CG Structure**

### 2. Chromatin Remodeler Slides Nucleosome Using ATP Hydrolysis
Chromatin Remodelers are SNF2 family ATPases. They contain 2 lobes in its ATPase domain. 
<p align="center">
  <img src="./snf2domain.png" width="1000">
</p>
**Figure 6. SNF2 ATPase Domain Structure Has 2 Catalytic Lobes**

Previous CG simulation by [Brandani and Takada,2018](https://doi.org/10.1101/297762) has simulated the ATP Hydrolysis process (apo->ATP->ADP) through changing the **Hydrogen bond strength** and **lobe1-lobe2 Go interaction strength** to facilitate conformational change that allow the lobe2 to close up to lobe1. Details such as the Energy Calculation Functions of Go Interaction can be found [here](https://www.charmm.org/wiki//index.php/Coarse_Grained_Go_Models#The_functional_form_of_the_coarse-grained_Go_model).
Inspired by their research, I adopted similar strategy to simulate DNA translocation by INO80 using ATP hydrolysis.
<p align="center">
  <img src="./atpaseModel.png" width="300">
</p>
**Figure 7. ATP Hydrolysis Cycle**


### 3. Define Translocase Lobe GO Contacts
The residues belong to lobe1 and lobe2 were identified using the [Cryo-EM structure 6fml](https://www.nature.com/articles/s41586-018-0029-y). This gives us the atpase in closed conformation. A simple homology model was constructed with 6hts to obtain the lobes in open conformation by aligning each of the two lobes with the respective ones from the 6fml. A list of go interactions between lobes 1 and 2 are then generated (script: "mda_define_atp_contacts.py").
The strength of these contacts in each translocate state is set according to [Brandani and Takada,2018](https://doi.org/10.1101/297762):
- apo: 0.0 kcal/mol (preferably open)
- ATP: 0.6 kcal/mol (preferably closed)
- ADP: 0.1 kcal/mol (preferably open)

### 4. Define DNA-Translocase Hydrogen Bonds
The hydrogen bonds are based on the all-atom (charmm27) conformation of the translocase
in complex with nucleosomal DNA as found in PDB 6fml with pdb2gmx. Hydrogen bonds with a cutoff distance of 5Å are identified (script: "mda_define_pdnsSnf2_cut5selmin.py").
The strength of hydrogen bonds in each translocate state is set according to [Brandani and Takada,2018]
- apo: 1.8 kBT (preferably attached to DNA)
- ATP: 1.8 kBT (preferably attached to DNA)
- ADP lobe1: 1.8 kBT  (preferably attached to DNA)
- ADP lobe2: 1.44 kBT  (preferably unattach and translocate forward)

### 5. Simulation
I simulated the DNA sliding activity of INO80 by running a simulation of an ATPase cycle with 3 states (apo-ATP-ADP), one after another. Each state was ran 10<sup>7</sup> timesteps with interactions strength modified as described in section 3 an 4 above. This CGMD simulation is performed using [CafeMOl3.0](https://doi.org/10.1021/ct2001045) with equation of motion via Constant temperature Langevin dynamics at 300K.

## IV. Results

## V. Running the simulation
One will need to download and configure Cafemol software to run this simulation.
To run the simulation, type:
cafemol p1.inp
