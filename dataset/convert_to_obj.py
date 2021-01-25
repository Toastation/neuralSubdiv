import sys

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