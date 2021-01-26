import sys
import os


def from_args():
    vert_path = sys.argv[1]
    tri_path = sys.argv[2]

    print(vert_path)
    print(tri_path)

    if not vert_path.endswith("vert") or not tri_path.endswith("tri"):
        print("USAGE: python convert_to_obj.py mesh.vert mesh.tri (respect arg order)")
        exit(1)

    name = vert_path.split(".")[0]

    vert_file = open(vert_path)
    vert_lines = vert_file.readlines()
    vert_file.close()

    tri_file = open(tri_path)
    tri_lines = tri_file.readlines()
    tri_file.close()

    obj_file = open(name + ".obj", "w")

    for line in vert_lines:
        obj_file.write("v " + line)

    for line in tri_lines:
        obj_file.write("f " + line)

    obj_file.close()

def from_cur_dir():
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    names = set()
    for f in files:
        if f.endswith(".py"): continue
        names.add(f.split(".")[0])
    for name in names:
        print("converting "+name)
        convert_file(name)

def convert_file(name):
    vert_path = name + ".vert"
    tri_path = name + ".tri"

    vert_file = open(vert_path)
    vert_lines = vert_file.readlines()
    vert_file.close()

    tri_file = open(tri_path)
    tri_lines = tri_file.readlines()
    tri_file.close()

    obj_file = open(name + ".obj", "w")

    for line in vert_lines:
        obj_file.write("v " + line)

    for line in tri_lines:
        obj_file.write("f " + line)

    obj_file.close()

from_cur_dir()