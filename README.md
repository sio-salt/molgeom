# molgeom ⌬
Easy to use molecular geometry manipulation Python library.

This is the version that requires numpy.
If you want to use the version that does not require numpy, please check the `main` branch.

## molgeom currently can do:
1. Reading and writing of molecular geometries of various formats:
    - XYZ, Gaussian input, GAMESS input, CIF, POSCAR
2. Calculations of molecular properties:
    - Nuclear repulsion energy
    - Center of mass
    - Bond lengths
3. Geometry manipulations:
    - Translation
    - Rotation (supports matrix and Rodrigues' rotation formula)
    - Replication of cells using lattice vectors
    - Replication of cells using symmetry operations (e.g. '-x, y + 1/2, -z + 1/2', ‘-2y+1/2, 3x+1/2, z-y+1/2’) 
    - Merging molecule geometries
4. Other useful functionalities:
    - Filtering atoms by element
    - Fragments clustering
    - Get bond cycles


## Installation
```bash
pip install -U git+https://github.com/sio-salt/molgeom@numpyfy
```
The minimul Python version is 3.11.


## Future plans
- Add Z-matrix support
- Add multi-molecule support
- much more...

