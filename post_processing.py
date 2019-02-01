# coding: utf-8
import numpy as np
import paramiko
import matplotlib.pyplot as plt
from matplotlib import rc
import os

def main():
    password = input("input password of FMa037@afifep.ifs.tohoku.ac.jp")
    path = "/work/A/FMa/FMa037/Case2/ResultC/"
    save_path = "G:\\Toyota\\Data\\Compressible_Viscos\\training_data\\NACA4\\"
    number = 6  # data number per wing
    
    max_angle = 39  # [deg]
    min_angle = 0  # [deg]
    delta_angle = 3  # [deg]
    angle_variation = int((max_angle - min_angle) / delta_angle) + 1
    
    wing_valiation = 1620
    total_data = wing_valiation * angle_variation
    
    save_data = np.zeros((total_data, number), dtype = float)
    
    type = 3.0
    save_data[:, 0] = type
    
    with paramiko.SSHClient() as ssh:
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(hostname = "afifep.ifs.tohoku.ac.jp", port = 22, username = "FMa037", password = password)
        
        count = 0
        for i1 in range(1, 10):
            for i2 in range(1, 10):
                for i34 in range(12, 92, 4):
                    naca4 = str(i1).zfill(1) + str(i2).zfill(1) + str(i34).zfill(2)
                    for deg in range(min_angle, max_angle + delta_angle, delta_angle):
                        save_data[count, 1] = float(naca4)
                        save_data[count, 3] = float(deg)
                        
                        fname = "NACA" + naca4 + "_" + str(deg).zfill(2) + "_AC.dat"
                        stdin, stdout, stderr = ssh.exec_command('cat ' + path + fname)
                        cd_cl = stdout.readlines()[0].split()
                        save_data[count, 4] = cd_cl[1]  # CL
                        save_data[count, 5] = cd_cl[0]  # CD
                        count += 1
    
    save_fname = save_path + "NACA4\\s1122_e9988_s4_a" + str(angle_variation).zfill(3) + ".csv"
    np.savetxt(save_fname, save_data, delimiter = ",")
    
def cdcl_plot_test():
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
        

def cp_plot_test():
    path = "D:\\Toyota\\github\\2D_comp_viscos\\flow_solver\\EulerSolver2_2018\\ResultC\\test_NACA0012\\"
    # fname = "CP_NACA" + str(i12digit).zfill(2) + str(i34digit).zfill(2) + "_" + str(angle).zfill(2) + ".dat"
    fname = "CP_NACA0012_fine_rev2_" + str(angle).zfill(2) + ".dat"
    char_angle = np.zeros((9, 2))
    i12digit = 00
    i34digit = 12
    # for angle in range(-3, 42, 3):
    angle = 10
    
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


if __name__ == '__main__':
    main()
    # rc('text', usetex=True)
    # cp_plot_test()
    # cdcl_plot_test()
    