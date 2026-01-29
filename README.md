# SWIS: Spin Wave and Inelastic Scattering

**SWIS** is a code package for computing **spin-wave dispersion** and **spin-polarized inelastic scattering spectra** of **noncollinear magnetic systems**.  
It was written and developed by **Flaviano José dos Santos**.

---

## Features
- Calculation of spin-wave dispersion for general noncollinear magnetic states.
- Computation of spin-polarized inelastic neutron scattering spectra.
- Flexible setup for different lattice geometries and interaction parameters.
- Designed for research on frustrated magnets, complex spin textures, and more.

---

## Compilation

### Using the Intel fortran compiler (ifx)

Download SWIS and in the terminal, go inside of the code folder.

#### Option 1 - Using cmake:

```
mkdir build
cd build
cmake ../
make
```

If configurations is needed, you should modify "CMakeLists.txt".

#### Option 2 - Using make:

```
make clean
make
```

The executable will created inside of a "build" folder.

#### Option 3 - Manual compilation:

```
ifx -c mod_global.f90 -qmkl -fopenmp -llapack -lblas
ifx -c bogoliubov_transf.f90 -qmkl -fopenmp -llapack -lblas
ifx -c spinwavesNonCol.f90 -qmkl -fopenmp -llapack -lblas
ifx -o main.exe mod_global.o bogoliubov_transf.o spinwavesNonCol.o -qmkl -fopenmp
```

---

## Citation

If SWIS was useful in any way for your research, please cite:

dos Santos, F. J. et al., Phys. Rev. B 97, 024431 (2018)\
https://doi.org/10.1103/PhysRevB.97.024431

dos Santos, F. J., PhD Thesis, RWTH Aachen University & Forschungszentrum Jülich GmbH (2019)\
DOI: 10.18154/RWTH-2020-01879\
https://publications.rwth-aachen.de/record/782441


## Contact
**Flaviano José Marchiori dos Santos, Dr. rer. nat.**\
Researcher at the Department of Theoretical Physics (COTEO)\
Brazilian Center for Research in Physics (CBPF)\
Ministry of Science, Technology, and Innovation (MCTI)\
Rua Dr.Xavier Sigaud, 150 - Urca, 22290-180, Rio de Janeiro - RJ, Brazil\
E-mail: flaviano@cbpf.br


## Installing Intel Fortran Compiler (ifx):

If this does not work, try to find instructions in the Intel website.

```
sudo apt update
sudo apt install curl
sudo cp /usr/share/keyrings/oneapi-archive-keyring.gpg /usr/share/keyrings/oneapi-archive-keyring.gpg_bkp
curl -Lo- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | sudo gpg --dearmor -o /usr/share/keyrings/oneapi-archive-keyring.gpg
sudo tee /etc/apt/sources.list.d/oneAPI.list <<< "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main"
sudo apt update
sudo apt install intel-oneapi-compiler-fortran
source /opt/intel/oneapi/setvars.sh (also add that line on ~/.bashrc)
ifx --version
```

