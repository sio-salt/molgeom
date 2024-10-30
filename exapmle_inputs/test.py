from molgeom import parse_file

# Read the xyz file
mole = parse_file("1.input_example.xyz")

mole.display_bonds(1, 2)
