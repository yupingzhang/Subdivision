# Subdivision

Loop subdivision code
- add boundary case for open mesh
- comment out the smoothing part for edge vertices

Compile:
clang++ subdivide.cpp slVector.cpp -o subdivide

Run:
./subdivide mesh.obj output.obj subdivide_num
