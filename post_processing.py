# coding: utf-8
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
import os

def main():
    path = "D:\\Toyota\\github\\2D_comp_viscos\\flow_solver\\EulerSolver2_2018\\ResultC\\test_NACA0012\\"
    for i in range(1800):
        i12digit = 0 #i % 100
        i34digit = 12 #5 * int(i / 100.0) + 11
        # char_angle = np.zeros((14, 2))
        char_angle = np.zeros((15, 2))
        # for angle in range(0, 42, 3):
        for angle in range(-3, 42, 3):
            deg = int(angle / 3) + 1

            # fname = "NACA" + str(i12digit).zfill(2) + str(i34digit).zfill(2) + "_" + str(angle).zfill(2) + "_AC.dat"
            fname = "NACA0012_fine_rev2_" + str(angle).zfill(2) + "_AC.dat"
            
            if os.path.exists(path + fname):
                # print(fname + " " + str(os.path.exists(path + fname)))
                # drag = np.zeros(1, dtype = float)
                # lift = np.zeros(1, dtype = float)
                number = 0
                with open(path + fname, "r") as f:
                    for line in f:
                        drag_and_lift = line.split()
                        attack = - angle / 180.0 * np.pi * 2.0  # 0.0
                        # drag[number] = float(drag_and_lift[0])
                        # lift[number] = float(drag_and_lift[1])
                        # number += 1
                        drag = float(drag_and_lift[0])
                        lift = float(drag_and_lift[1])
            
            # char_angle[deg, 0] = - np.average(drag)
            # char_angle[deg, 1] = np.average(lift)
                char_angle[deg, 0] = - drag
                char_angle[deg, 1] = lift
                print(-drag)
                
        
        # t = np.linspace(0, 39, 14)
        t = np.linspace(-3, 39, 15)
        plt.plot(t, char_angle[:, 0], "x")
        plt.plot(t, char_angle[:, 1], "o")
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title("NACA0012 Fine-Grid " + r'$C_D$' + " & " + r'$C_L$' + " vs Angle of Attack")
        ax.plot(t, char_angle[:, 0], "x", label="Drag Coefficient: Re=5 mill. Ma=0.7")
        ax.plot(t, char_angle[:, 1], "o", label="Lift Coefficient: Re=5 mill. Ma=0.7")
        ax.legend()
        ax.set_xlabel(r"Angle of Attack [deg]")
        ax.set_ylabel(r"$C_D$" + " & " + r"$C_L$")
        ax.grid()
        ax.set_xlim(-3.0, 39.0)
        ax.set_ylim(-0.5, 1.5)
        plt.show()

        plt.show()
        exit()
        

def cp_graph():
    path = "D:\\Toyota\\github\\2D_comp_viscos\\flow_solver\\EulerSolver2_2018\\ResultC\\test_NACA0012\\"
    char_angle = np.zeros((9, 2))
    i12digit = 00
    i34digit = 12
    # for angle in range(-3, 42, 3):
    angle = 10
    # fname = "CP_NACA" + str(i12digit).zfill(2) + str(i34digit).zfill(2) + "_" + str(angle).zfill(2) + ".dat"
    fname = "CP_NACA0012_fine_rev2_" + str(angle).zfill(2) + ".dat"
    coord = np.zeros(5000, dtype = float)
    pressure = np.zeros(5000, dtype = float)
    print(path + fname)
    if os.path.exists(path + fname):
        # print(fname + " " + str(os.path.exists(path + fname)))
        number = 0
        with open(path + fname, "r") as f:
            for line in f:
                coord_and_pressure = line.split()
                if(not coord_and_pressure):
                    break
                attack = - angle / 180.0 * np.pi * 2.0  # 0.0
                coord[number] = float(coord_and_pressure[0])
                pressure[number] = float(coord_and_pressure[1])
                # print(coord_and_pressure)
                # print(angle)
                number += 1

        alpha = r'$\alpha$'
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title("NACA0012 Fine-Grid " + alpha + " = 10 deg")
        ax.plot(coord[:number], pressure[:number], ".", label="Toyota. Re=5 mill. Ma=0.7")
        ax.legend()
        ax.set_xlabel(r"$\frac{x}{c}$")
        ax.set_ylabel(r"$C_P$")
        ax.grid()
        ax.set_ylim(2.0, -6.0)
        plt.show()
            # exit()
    
def test():
    i12list = [1,2,3,4,5,6,7,8,9]
    for i in range(0, 1620):
        i12 = int(float(i) / 20.0) + 11 + int(float(i) / 180.0)
        # print(i, i1, int(float(i) / 180.0))
        i34 = 4 * (i % 20) + 12
        print(str(i12).zfill(2) + str(i34).zfill(2))
        

if __name__ == '__main__':
    # rc('text', usetex=True)
    # main()
    # cp_graph()
    test()