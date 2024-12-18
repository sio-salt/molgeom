# molgeom
Easy to use molecular geometry manipulation Python library.
Makes your quantum chemistry life easier.


### molgeom currently can do:
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


### Installation
```bash
pip install -U git+https://github.com/sio-salt/molgeom
```

