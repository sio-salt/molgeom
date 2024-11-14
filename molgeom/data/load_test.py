import yaml


p_yaml = "periodic_table.yaml"
with open(p_yaml, "r") as file:
    data = yaml.safe_load(file)

ag = data["Ag"]["Atomic mass"]
print(ag)
