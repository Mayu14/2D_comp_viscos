# -- coding: utf-8 --
import os
import re

def main():
    CFLAGS_ifort = "CFLAGS = -Wall -check uninit -check pointers -check bounds -check all -inline-level=0 -O0 -traceback -warn all -ftrapuv -debug full -zero -check uninit -check pointers -check bounds -check all -O0 -traceback -warn interfaces -warn all -ftrapuv -debug full  -module $(OBJS_DIR) $(IDIR)\n"
    LFLAGS_ifort = "LFLAGS = -s\n"

    os.replace("Makefile", "Makefile_win")
    flag = 0
    with open(file="Makefile_win", mode="r") as fw:
        with open(file="Makefile", mode="w") as fl:
            for row in fw:
                line = row

                if line.startswith("EXE"):
                    line = line.replace(".exe", "")
                    flag = 1

                if flag == 1:
                    if line.startswith("FC") or line.startswith("LD"):
                        line = line.replace("gfortran.exe", "ifort")

                    if line.startswith("CFLAGS"):
                        line = CFLAGS_ifort

                    if line.startswith("LFLAGS"):
                        line = LFLAGS_ifort
                        flag = 0
                fl.write(line)


if __name__ == '__main__':
    main()