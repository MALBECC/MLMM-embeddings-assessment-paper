# MLMM-embeddings-assessment-paper
This repository contains data and software regarding the paper submited to JCIM, entitled "Assessment of embedding schemes in a hybrid machine learning/classical potentials (ML/MM) approach". 

## Content
- The python scripts necessary to reproduce all the Figures of this manuscript, as well as a Jupiter notebook to do so.
- A stand-alone python script to compute the electrostatic ML/MM coupling energy with the Thole model.
- The AMBER patch developed to perform the ML/MM simulations presented here as a proof of concept.

## Usage
The pythole code and the notebook require the "ANI-aa-qmmm.h5" data set to be present.

## Notebook
This notebook requires the "ANI-aa-qmmm.h5" data set to be present in the same directory it is executed (https://doi.org/10.6084/m9.figshare.24582738.v2).
It contains the necesary code to reproduce all the Figures presented in the paper.

## PyThole
This is python script that calculates the Thole energy of all the structures in the "ANI-aa-qmmm.h5" data set. For each formula entry in the data set, it creates a .dat file containing the Thole energy for every structure.

## AMBER ML/MM Patch
A patch for a proof of concept on ML/MM simulations in AMBER. The implementation works with a File-Based-Interface with Torchani and Psi4. The code works by combining the ANI-2x potential and Psi4-predicted atomic charges in each forces calculation, allowing to run molecular dynamics and minimizations with Sander, using a ML/MM scheme. While this implementation is not computationally efficient ficient (since it uses an actual QM method to compute the atomic charges), as this property would be ML-predicted in an actual implementation of our model. 

### Requirements
- Python 3
- Torch, Torchani, NumPy
- Psi4

### Installation

```
./patch.sh /path/to/sander_folder
```

### Example 

See below an example `mdin` input file run an NVT MD with this interface as a ML/MM engine in amber:

```
 &cntrl
    imin=0, ntx=5, ntwr=100,ntpr=1,
    ntwx=10,ioutfm=1,ntxo=1,
    nstlim=5000,dt=0.001,
    ntt=3,tempi=300.0,temp0=300.0,gamma_ln=5.0,
    ntp=0,ntb=1,ntf=1,ntc=2,
    cut=10.0,
    ifqnt=1, !turn on QMMM/MLMM calculations
  &end
/
 &qmmm
  qm_theory='EXTERN',
  qmmask=':1',
  qmmm_int =2, !set to 2 to turn on the ML/MM coumpling
  qmshake=0,
  qm_ewald=0,
  qmcut=15.0,
  writepdb=1,
  verbosity=2,
/
 &ml
  mlmm_patch_path='/path/to/amber-mlmm-path_folder'
  python_executable='python' !python excecutable
  psi4_path='/path/to/psi4_executable' !should be psi4_executable
  psi4_nthreads=16 !number of threads used for the Psi4 calculation of MBIS_CHARGES  
  charges_functional='wb97x' !level of theory for the Psi calculation (wb97x recommended)
  charges_basis='6-31G(d)'   !level of theory for the Psi calculation (6-31G(d) recommended)

```

### References
We refer the users to the manuscript.


