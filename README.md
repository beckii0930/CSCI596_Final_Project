# CSCI596 Final Project: Coarse-Grain MD Simulation for Chromatin Remodeler INO80



## 0. Background

### Chromatin Remodelers
<p align="center">
  <img src="./ino80.png" width="800">
</p>


### Coarse Graining

A "3 Sites per Nucleotide" strategy is used to coarse grain DNA. Below is a sample CG scheme for a nucleotide C.
<p align="center">
  <img src="./Ccgscheme.png" width="500">
</p>

A "1 Sites per Amino Acid" strategy is used to coarse grain protein.Below is a sample CG scheme for a protein sequence ALA-ASP-ASP-GLU-GLN.
<p align="center">
  <img src="./pro_cg_scheme.png" width="500">
</p>

Now we obtain the entire Initial CG structure for this system.
![image](https://user-images.githubusercontent.com/25398675/143984154-7b7f0b93-97b7-4076-8595-bdf312867ebc.png)


## 1. Introduction
If the C compiler on your computer is cc (also common is gcc for Gnu C
compiler), type:
cc -O -o md md.c -lm
This will create an executable named md. To run the executable, type:
./md < md.in
## 2. Methods
The following files are included in this folder, in addition to this readme
file, readme.md.
<ul>
<li>md.c: Main C program</li>
<li>md.h: Header file for md.c</li>
<li>md.in: Input parameter file (to be redirected to the standard input)</li>
</ul>
![Screen shot of MD simulation](ScreenShot.png)

## 2. Results
