# molgeom ⌬
Easy to use molecular geometry manipulation Python library.

There is the numpyfy version of molgeom.
If you want to use the numpyfy version, please check the [numpyfy branch](https://github.com/sio-salt/molgeom/tree/numpyfy).

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

More functionalities are coming soon!


## Installation
```bash
pip install -U git+https://github.com/sio-salt/molgeom/tree/main
```
or
```bash
git clone https://github.com/sio-salt/molgeom.git
cd molgeom
pip install -e .
```
The minimul Python version is 3.9.


## Future plans
- Add example codes
- Add Z-matrix support
- Add multi-molecule support
- much more...


