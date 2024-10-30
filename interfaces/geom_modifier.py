import sys
import importlib
from easyvec import Vec3
import molgeom
from molgeom import parse_file

importlib.reload(molgeom)
if importlib.util.find_spec("readline"):
    import readline  # noqa


def translate_molecule(mole):
    print()
    print("---------------")
    print("  translation  ")
    print("---------------")
    print("enter coordinate A and B  (translate A -> B)")
    while True:
        try:
            pA = input("coordinate A:  ( e.g.  1.2  0.55  3.0 )\n > ")
            pA = Vec3(*list(map(float, pA.strip().split())))
            pB = input("coordinate B:  ( e.g.  1.2  0.55  3.0 )\n > ")
            pB = Vec3(*list(map(float, pB.strip().split())))
            break
        except ValueError:
            print("Please enter valid coordinates.")
            continue
        except TypeError as e:
            print(f"Parsing error: {e}")
            continue
    trans_vec = pB - pA
    print(f"trans vector : \n {trans_vec}")
    mole.translate(*trans_vec)

    return mole


def mirror_molecule(mole):
    print()
    print("---------------")
    print("  reflection   ")
    print("---------------")
    print("define mirror plane with 3 points")
    while True:
        try:
            inp1 = input("point 1:  ( e.g.  1.2  0.55  3.0 )\n > ")
            p1 = Vec3(*list(map(float, inp1.strip().split())))
            inp2 = input("point 2:  ( e.g.  1.2  0.55  3.0 )\n > ")
            p2 = Vec3(*list(map(float, inp2.strip().split())))
            inp3 = input("point 3:  ( e.g.  1.2  0.55  3.0 )\n > ")
            p3 = Vec3(*list(map(float, inp3.strip().split())))
            break
        except ValueError:
            print("Please enter valid points.")
            continue
        except TypeError as e:
            print(f"Parsing error: {e}")
            continue

    mole.mirror_by_plane(p1, p2, p3)

    # sign_xyz = input("enter mirror signs  ( e.g.  1  -1  1 )\n > ")
    # sign_xyz = list(map(int, sign_xyz.strip().split()))
    # mole.mirror(*sign_xyz)

    return mole


def rotate_molecule(mole):
    print()
    print("---------------")
    print("   rotation    ")
    print("---------------")
    while True:
        try:
            axis_point1 = input("enter axis point1  ( e.g.  1.2  0.55  3.0 )\n > ")
            axis_point1 = Vec3(*list(map(float, axis_point1.strip().split())))
            axis_point2 = input("enter axis point2  ( e.g.  1.2  0.55  3.0 )\n > ")
            axis_point2 = Vec3(*list(map(float, axis_point2.strip().split())))
            break
        except ValueError:
            print("Please enter valid points.")
            continue
        except TypeError as e:
            print(f"Parsing error: {e}")
            continue

    rot_axis = axis_point2 - axis_point1
    print("rotation axis : \n", rot_axis)

    angle_degrees = float(input("enter rotation angle (degrees)  ( e.g.  90 )\n > "))
    mole.rotate_by_axis(axis_point1, axis_point2, angle_degrees)

    return mole


operation_functions = {
    "translate": translate_molecule,
    "reflect": mirror_molecule,
    "rotate": rotate_molecule,
}


def get_operation_order():
    operations = ["translate", "reflect", "rotate"]
    while True:
        try:
            input_str = input(
                "---------------------------------------------------\n"
                " Enter the order of operations separated by spaces \n"
                "     0 if you don't want to apply the operation    \n"
                "   translate  reflect  rotate   ( e.g.  2  0  1 )  \n"
                "---------------------------------------------------\n"
                " > "
            )
            values = [int(v.strip()) for v in input_str.strip().split()]
            if len(values) != 3:
                print("Please enter 3 numbers.")
                continue
            if any(v < 0 or v > 3 for v in values):
                print("Please enter numbers between 0 and 3.")
                continue
            if len(set(v for v in values if v != 0)) != len(
                [v for v in values if v != 0]
            ):
                print("Please enter non-duplicate numbers.")
                continue
            non_zero_values = [v for v in values if v != 0]
            if sorted(non_zero_values) != list(range(1, len(non_zero_values) + 1)):
                print(
                    "The order of operations is not consecutive. Please enter consecutive numbers starting from 1."
                )
                continue
            order = dict(zip(operations, values))
            break
        except ValueError:
            print("Please enter integers.")
    return order


def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <xyz_file_path>")
        sys.exit(1)

    filepath = sys.argv[1]
    print(f"\n{filepath}\n")
    mole = parse_file(filepath, "r")

    operation_order = get_operation_order()
    print()
    if all(v == 0 for v in operation_order.values()):
        print("bye")
        print()
        sys.exit(0)

    print(f"order of operations: {operation_order}")

    for i in range(1, 4):
        for op, order in operation_order.items():
            if order == i:
                mole = operation_functions[op](mole)

    print()
    print(" final molecule geometry ")
    print()
    print(mole.to_xyz())
    print()


if __name__ == "__main__":
    main()
