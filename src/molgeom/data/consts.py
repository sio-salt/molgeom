"""
future plan:
    https://github.com/materialsproject/pymatgen/blob/master/dev_scripts/periodic_table.yaml
    more detailed info
"""

# from Gaussian homepage
ANGST2BOHR_GAU16 = 1.0 / 0.52917721092
ANGST2BOHR_GAU09 = 1.0 / 0.52917720860
ANGST2BOHR_GMS19 = 1.0 / 0.52917724924

HARTREE2KCAL_SIMPLE = 627.51
HARTREE2KCAL_GAU09 = 627.509541
HARTREE2KCAL_GMS19 = 627.509469  # from GAMESS
HARTREE2KCAL_WIKI = 627.5094740631  # from Wikipedia

# fmt: off
ATOMIC_NUMBERS = {
    "H" :  1, "He":  2, "Li":  3, "Be":  4, "B" :  5, "C" :  6,
    "N" :  7, "O" :  8, "F" :  9, "Ne": 10, "Na": 11, "Mg": 12,
    "Al": 13, "Si": 14, "P" : 15, "S" : 16, "Cl": 17, "Ar": 18,
    "K" : 19, "Ca": 20, "Sc": 21, "Ti": 22, "V" : 23, "Cr": 24,
    "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y" : 39, "Zr": 40, "Nb": 41, "Mo": 42,
    "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48,
    "In": 49, "Sn": 50, "Sb": 51, "Te": 52, "I" : 53, "Xe": 54,
    "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
    "Ho": 67, "Er": 68, "Tm": 69#, "Yb": 70, "Lu": 71, "Hf": 72,
    # "Ta": 73, "W" : 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
    # "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82,
}

NUMBER_TO_SYMBOL = {v: k for k, v in ATOMIC_NUMBERS.items()}

SPECIAL_ELEMENTS = {
    "X", "Tv" 
}

ATOMIC_MASSES = {
    "H" :  1.00790, "He":  4.00260, "Li":  6.94000, "Be":  9.01218, "B" : 10.81000, "C" : 12.01100,
    "N" : 14.00670, "O" : 15.99940, "F" : 18.99840, "Ne": 20.17900, "Na": 22.98977, "Mg": 24.30500,
    "Al": 26.98154, "Si": 28.08550, "P" : 30.97376, "S" : 32.06000, "Cl": 35.45300, "Ar": 39.94800,
    "K" : 39.09830, "Ca": 40.08000, "Sc": 44.95590, "Ti": 47.90000, "V" : 50.94150, "Cr": 51.99600,
    "Mn": 54.93800, "Fe": 55.84700, "Co": 58.93320, "Ni": 58.71000, "Cu": 63.54600, "Zn": 65.38000,
    "Ga": 69.73500, "Ge": 72.59000, "As": 74.92160, "Se": 78.96000, "Br": 79.90400, "Kr": 83.80000,
    "Rb": 85.46780, "Sr": 87.62000, "Y" : 88.90590, "Zr": 91.22000, "Nb": 92.90640, "Mo": 95.94000,
    "Tc": 98.90620, "Ru": 101.0700, "Rh": 102.9055, "Pd": 106.4000, "Ag": 107.8680, "Cd": 112.4100,
    "In": 114.8200, "Sn": 118.6900, "Sb": 121.7500, "Te": 127.6000, "I" : 126.9045, "Xe": 131.3000,
    "Cs": 132.9054, "Ba": 137.3300, "La": 150.0000, "Ce": 178.4900, "Pr": 180.9479, "Nd": 183.8500,
    "Pm": 186.2070, "Sm": 190.2000, "Eu": 192.2200, "Gd": 195.0900, "Tb": 196.9665, "Dy": 200.5900,
    "Ho": 204.3700, "Er": 207.2000, "Tm": 208.9804
}
