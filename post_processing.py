# coding: utf-8
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

def main():
    path = "D:\\Toyota\\github\\2D_comp_viscos\\flow_solver\\EulerSolver2_2018\\ResultC\\"
    for i in range(1800):
        i12digit = i % 100
        i34digit = 5 * int(i / 100.0) + 11
        char_angle = np.zeros((14, 2))
        for angle in range(0, 42, 3):
            deg = int(angle / 3)
            fname = "NACA" + str(i12digit).zfill(2) + str(i34digit).zfill(2) + "_" + str(angle).zfill(2) + "_AC.dat"
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
        
        t = np.linspace(0, 39, 14)
        plt.plot(t, char_angle[:, 0])
        plt.plot(t, char_angle[:, 1])
        plt.show()
        exit()
        

                
if __name__ == '__main__':
    main()
