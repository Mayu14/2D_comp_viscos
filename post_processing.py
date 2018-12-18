# coding: utf-8
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def main():
    path = "D:\\Toyota\\github\\2D_comp_viscos\\flow_solver\\EulerSolver2_2018\\ResultC\\"
    for i in range(1800):
        i12digit = 4#i % 100
        i34digit = 11#5 * int(i / 100.0) + 11
        # char_angle = np.zeros((14, 2))
        char_angle = np.zeros((9, 2))
        # for angle in range(0, 42, 3):
        for angle in range(8, 28, 2):
            deg = int(angle / 3)
            # fname = "NACA" + str(i12digit).zfill(2) + str(i34digit).zfill(2) + "_" + str(angle).zfill(2) + "_AC.dat"
            fname = "NACA0012_course_" + str(angle).zfill(2) + "_AC.dat"
            if os.path.exists(path + fname):
                # print(fname + " " + str(os.path.exists(path + fname)))
                drag = np.zeros(100, dtype = float)
                lift = np.zeros(100, dtype = float)
                number = 0
                with open(path + fname, "r") as f:
                    for line in f:
                        drag_and_lift = line.split()
                        attack = - angle / 180.0 * np.pi * 2.0  # 0.0
                        drag[number] = float(drag_and_lift[0])
                        lift[number] = float(drag_and_lift[1])
                        number += 1
            
            char_angle[deg, 0] = - np.average(drag)
            char_angle[deg, 1] = np.average(lift)
        
        # t = np.linspace(0, 39, 14)
        t = np.linspace(8, 26, 9)
        plt.plot(t, char_angle[:, 0])
        plt.plot(t, char_angle[:, 1])
        plt.show()
        exit()
        

def cp_graph():
    path = "D:\\Toyota\\github\\2D_comp_viscos\\flow_solver\\EulerSolver2_2018\\ResultC\\"
    char_angle = np.zeros((9, 2))
    for angle in range(8, 28, 2):
        deg = int(angle / 3)
        # fname = "NACA" + str(i12digit).zfill(2) + str(i34digit).zfill(2) + "_" + str(angle).zfill(2) + "_AC.dat"
        fname = "CP_NACA0012_course_" + str(angle).zfill(2) + ".dat"
        coord = np.zeros(500, dtype = float)
        pressure = np.zeros(500, dtype = float)
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
                    pressure[number] = -float(coord_and_pressure[1])
                    number += 1
            
            plt.plot(coord[:number], pressure[:number], "x")
            plt.show()
            exit()
        

if __name__ == '__main__':
    # main()
    cp_graph()
