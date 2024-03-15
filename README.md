# MLMM-embeddings-assessment-paper
This repository contains data and software regarding the paper submited to JCIM, entitled "Assessment of embedding schemes in a hybrid machine learning/classical potentials (ML/MM) approach". 

## Content
- The python scripts necessary to reproduce all the Figures of this manuscript, as well as a Jupiter notebook to do so.
- A stand-alone python script to compute the electrostatic ML/MM coupling energy with the Thole model.

## Usage
The pythole code and the notebook require the "ANI-aa-qmmm.h5" data set to be present (https://doi.org/10.6084/m9.figshare.24582738.v2).

## Notebook
This notebook requires the "ANI-aa-qmmm.h5" data set to be present in the same directory it is executed.
It contains the necesary code to reproduce all the Figures presented in the paper.

## PyThole
This is python script that calculates the Thole energy of all the structures in the "ANI-aa-qmmm.h5" data set. For each formula entry in the data set, it creates a .dat file containing the Thole energy for every structure.

### References
We refer the users to the manuscript.


