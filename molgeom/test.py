from molgeom import parse_file

file_path = "../exapmle_inputs/1.input_example.xyz"
mole = parse_file(file_path, 'r')

print(mole[3])
print(mole[3:5])
