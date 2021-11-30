# CSCI596 Final Project: Coarse-Grain MD Simulation for Chromatin Remodeler INO80
#### Yibei Jiang


## I. Background

### 1. Chromatin Remodeler IN080 (PDB ID: 6fml)
<p align="center">
  <img src="./ino80.png" width="800">
</p>

```
The Biological Question: 
What is the moleclar mechanism INO80 uses to slide the nucleosome using ATP hydrolysis?
```
### 2. Coarse Grain Molecular Dynamics

A "**3 Sites per Nucleotide**" strategy is used to coarse grain DNA. Below is a sample CG scheme for a nucleotide C.
<p align="center">
  <img src="./Ccgscheme.png" width="400">
</p>

A "**1 Site per Amino Acid**" strategy is used to coarse grain protein. Below is a sample CG scheme for a protein
sequence GLN-GLU-ASP-ASP-ALA.
<p align="center">
  <img src="./pro_cg_scheme.png" width="400">
</p>

Now we obtain the entire Initial CG structure for this system.
![image](https://user-images.githubusercontent.com/25398675/143984154-7b7f0b93-97b7-4076-8595-bdf312867ebc.png)


## II. Introduction
If the C compiler on your computer is cc (also common is gcc for Gnu C
compiler), type:
cc -O -o md md.c -lm
This will create an executable named md. To run the executable, type:
./md < md.in
## III. Methods
The following files are included in this folder, in addition to this readme
file, readme.md.
<ul>
<li>md.c: Main C program</li>
<li>md.h: Header file for md.c</li>
<li>md.in: Input parameter file (to be redirected to the standard input)</li>
</ul>
![Screen shot of MD simulation](ScreenShot.png)

## IV. Results
