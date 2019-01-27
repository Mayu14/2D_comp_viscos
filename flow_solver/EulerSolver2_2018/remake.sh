rm *__genmod.*
rm *.mod
rm -r obj/
python3 convert_makefile_for_linux.py
make
./bin/Debug/EulerSolver2
