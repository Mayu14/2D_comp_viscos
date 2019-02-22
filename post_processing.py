# coding: utf-8
import numpy as np
import paramiko
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import os
import glob
from itertools import cycle, islice

def main(online=False):
    def main_process(online):
        count = 0
        for i1 in range(1, 10):
            for i2 in range(1, 10):
                for i34 in range(12, 92, 4):
                    naca4 = str(i1).zfill(1) + str(i2).zfill(1) + str(i34).zfill(2)
                    for deg in range(min_angle, max_angle + delta_angle, delta_angle):
                        if(deg != 24):
                            save_data[count, 1] = float(naca4)
                            save_data[count, 3] = float(deg)
                            
                            fname = "NACA" + naca4 + "_" + str(deg).zfill(2) + "_AC.dat"
                            if online:
                                stdin, stdout, stderr = ssh.exec_command('cat ' + path + fname)
                                cd_cl = stdout.readlines()[0].split()
                            else:
                                with open(path + fname, "r") as f:
                                    cd_cl = f.readline().split()
                            
                            save_data[count, 4] = float(cd_cl[1])  # CL
                            save_data[count, 5] = -float(cd_cl[0])  # CD
                            count += 1
    
    save_path = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\NACA4\\"
    number = 6  # data number per wing

    if online:
        password = input("input password of FMa037@afifep.ifs.tohoku.ac.jp")
        path = "/work/A/FMa/FMa037/Case3/ResultC/"
    else:
        path = "G:\\Toyota\\Data\\Case3\\"
        
    max_angle = 39  # [deg]
    min_angle = 0  # [deg]
    delta_angle = 3  # [deg]
    angle_variation = int((max_angle - min_angle) / delta_angle) + 1 - 1    # 24degは除く
    
    wing_valiation = 1620
    total_data = wing_valiation * angle_variation
    
    save_data = np.zeros((total_data, number), dtype = float)
    
    type = 3.0
    save_data[:, 0] = type
    
    if online:
        with paramiko.SSHClient() as ssh:
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            ssh.connect(hostname = "afifep.ifs.tohoku.ac.jp", port = 22, username = "FMa037", password = password)
            main_process(online)
    else:
        main_process(online)
        
    
    save_fname = save_path + "s1122_e9988_s4_a" + str(angle_variation).zfill(3) + ".csv"
    np.savetxt(save_fname, save_data, delimiter = ",")
    
def cdcl_plot_test():
    path = "D:\\Toyota\\github\\2D_comp_viscos\\flow_solver\\EulerSolver2_2018\\ResultC\\test_NACA0012\\"
    path = "G:\\Toyota\\Data\\Case3\\"

    for i in range(1800):
        i12digit = 0 #i % 100
        i34digit = 12 #5 * int(i / 100.0) + 11
        # char_angle = np.zeros((14, 2))
        char_angle = np.zeros((15, 2))
        # for angle in range(0, 42, 3):
        for angle in range(-3, 42, 3):
            deg = int(angle / 3) + 1

            # fname = "NACA" + str(i12digit).zfill(2) + str(i34digit).zfill(2) + "_" + str(angle).zfill(2) + "_AC.dat"
            fname = "NACA1112_" + str(angle).zfill(2) + "_AC.dat"
            
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
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title("NACA1112 Course-Grid "+r'$C_D$'+" and "+r'$C_L$'+" vs Angle of Attack")
        ax.plot(t, char_angle[:, 0], "x", label="Drag Coefficient: Ma=0.7")
        ax.plot(t, char_angle[:, 1], "o", label="Lift Coefficient: Ma=0.7")
        ax.legend()
        ax.set_xlabel(r"Angle of Attack [deg]")
        ax.set_ylabel(r"$C_D$"+" and "+r"$C_L$")
        ax.grid()
        ax.legend()
        ax.set_xlim(-3.0, 39.0)
        ax.set_ylim(-0.25, 1.75)
        plt.show()

        exit()
        

def cp_plot_test():
    path = "D:\\Toyota\\github\\2D_comp_viscos\\flow_solver\\EulerSolver2_2018\\ResultC\\"
    # path = "G:\\Toyota\\Data\\Case3\\"
    # fname = "CP_NACA" + str(i12digit).zfill(2) + str(i34digit).zfill(2) + "_" + str(angle).zfill(2) + ".dat"
    angle = 15
    case_name = "_M015_Roe_LTS_CFL95"
    fname = "CP_NACA0012_" + str(angle).zfill(2) + case_name + ".dat"
    char_angle = np.zeros((9, 2))
    i12digit = 00
    i34digit = 12
    # for angle in range(-3, 42, 3):
    
    
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
        ax.set_title("NACA0012 Fine-Grid " + alpha + " = " + str(angle) +" deg Roe + LTS")
        ax.plot(coord[:number], pressure[:number], ".", label="Toyota. Ma=0.15")
        ax.legend()
        ax.set_xlabel(r"$\frac{x}{c}$")
        ax.set_ylabel(r"$C_P$")
        ax.grid()
        ax.set_ylim(2.0, -6.0)
        plt.show()
            # exit()

def plot_residual_graph():
    path = "G:\\Toyota\\Data\\Compressible_Invicid\\solver_validation\\NACA0012\\M015_15deg_compare\\ResultR\\"
    fname = "RES_NACA0012_15_M015_Roe_LTS_CFL95.csv"

    name = ["TimeStep", "Ave:Density", "Ave:Momentum_X", "Ave:Momentum_Y", "Ave:Momentum_Z", "Ave:Energy", "Ave:Mix",
            "Max:Density", "Max:Momentum_X","Max:Momentum_Y", "Max:Momentum_Z", "Max:Energy", "Max:Entire"]

    d_type = {"TimeStep":"int8", "Ave:Density":"float16", "Ave:Momentum_X":"float16", "Ave:Momentum_Y":"float16",
             "Ave:Momentum_Z":"float16", "Ave:Energy":"float16", "Ave:Mix":"float16", "Max:Density":"float16",
              "Max:Momentum_X":"float16", "Max:Momentum_Y":"float16", "Max:Momentum_Z":"float16",
              "Max:Energy":"float16", "Max:Entire":"float16"}

    df = pd.read_csv(path + fname, names=name, dtype = d_type)

    data_total = df.shape[0]


    angle = 15
    case_name = "Residual History of NACA0012 Fine-Grid\n" + r'$\alpha$' + " = " + str(angle) + " deg, Ma = 0.15 with Roe + LTS"
    timestep = np.arange(data_total)
    
    def make_subplot(locate=1, data_label="Max:Density"):
        ax = fig.add_subplot(3, 2, locate)
        if(data_label.find("Momentum_Z") == -1):
            ax.set_yscale("log")
        ax.set_title(data_label)
        ax.plot(timestep, df[data_label], ".", label = data_label)
        ax.legend()
        ax.set_xlabel(r"$Time Step$")
        ax.set_ylabel(r"$Error$")
        ax.grid()
        

    fig = plt.figure()
    fig.suptitle(case_name)
    
    ax = fig.add_subplot(1,1,1)
    data_label = "Max:Entire"
    #ax.set_yscale("log")
    ax.set_title(data_label)
    ax.plot(timestep[75000:], df[data_label][75000:], ".", label = data_label)
    ax.legend()
    ax.set_xlabel(r"$Time Step$")
    ax.set_ylabel(r"$Error$")
    ax.grid()
    plt.tight_layout()
    plt.show()
    
    make_subplot()
    make_subplot(2, "Max:Momentum_X")
    make_subplot(3, "Max:Momentum_Y")
    make_subplot(4, "Max:Momentum_Z")
    make_subplot(5, "Max:Energy")
    make_subplot(6, "Max:Entire")
    plt.tight_layout()
    plt.show()

    fig = plt.figure()
    fig.suptitle(case_name)
    make_subplot(1, "Ave:Density")
    make_subplot(2, "Ave:Momentum_X")
    make_subplot(3, "Ave:Momentum_Y")
    make_subplot(4, "Ave:Momentum_Z")
    make_subplot(5, "Ave:Energy")
    make_subplot(6, "Ave:Mix")
    plt.tight_layout()
    plt.show()

def make_restart_list(digit, deg_list):

    path = "D:\\Toyota\\Dropbox\\shareTH\\program\\"
    date = "02221855"
    fname = "Case4_C_" + date + ".txt"
    mid_fname = "Case4_U_" + date + ".txt"
    restart_fname = "re_cal_namelist.dat"
    restart_list = []
    restart_list_u = []
    # 再計算しなければならないケースのリスト作成
    for degree in deg_list:
        deg = str(degree).zfill(2)
        name_list = []
        name_list_u = []
        if digit == 4:
            for i1 in range(1,10):
                for i2 in range(1, 10):
                    for i34 in range(12, 92, 4):
                        naca4 = str(i1) + str(i2) + str(i34)
                        name = "CP_NACA" + naca4 + "_" + deg + "_M015.dat"
                        name_u = "NACA" + naca4 + "_" + deg + "_M015_"
                        name_list.append(name)
                        name_list_u.append(name_u)
        else:
            head_list = [210, 220, 230, 240, 241, 250, 251]
            for i123 in head_list:
                for i45 in range(11, 91):
                    naca5 = str(i123).zfill(3) + str(i45).zfill(2)
                    name = "CP_NACA" + naca5 + "_" + deg + "_M015.dat"
                    name_u = "NACA" + naca5 + "_" + deg + "_M015_"
                    name_list.append(name)
                    name_list_u.append(name_u)

        with open(path + fname, "r", encoding = "utf-8") as f:
            for line in f:
                if(line[8+digit:10+digit] == deg):
                    # print(line)
                    if (line[:19+digit] in name_list):
                        name_list.remove(line[:19+digit])


        restart_list.extend(name_list)

    total = len(restart_list)

    # 中間出力の名前に変換
    for i in range(total):
        restart_list[i] = restart_list[i][3:15+digit] + "_"

    # 中間出力が存在すればファイル名，存在しなければNoneを返す
    mid_file = [None]*len(restart_list)
    for degree in deg_list:
        deg = str(degree).zfill(2)

        with open(path + mid_fname, encoding = "utf-8") as f:
            for line in f:
                if(line[5+digit:7+digit] == deg):
                    if(line[:13+digit] in restart_list):
                        mid_file[restart_list.index(line[:13+digit])] = line

    with open(path + restart_fname, "w") as f:
        f.write(str(total) + " case Required.\n")
        for i in range(total):
            name = restart_list.pop(0)
            naca4 = name[4:4+digit]
            deg = name[5+digit:7+digit]
            f.write(naca4 + " " + deg + " " + str(mid_file.pop(0)).replace("\n", "") + "\n")
    exit()

def del_duplication(cp=False):
    def rename_on_win(f_n, new_name):
        if f_n != new_name:
            if (os.path.exists(new_name)):
                os.remove(new_name)
            
            os.rename(f_n, new_name)
            
    # path = "/work/A/FMa/FMa037/Case4/ResultC/"
    path = "G:\\Toyota\\Data\\Compressible_Invicid\\Case4\\ResultC\\"
    p_n = len(path)
    if cp:
        dir = glob.glob(path + "CP*")
        for f_n in dir:
            if (f_n[p_n+11] == "_"):    # 4digit
                new_name = f_n[:p_n + 19] + ".dat"
                rename_on_win(f_n, new_name)
            else:   # 5digit
                new_name = f_n[:p_n + 20] + ".dat"
                rename_on_win(f_n, new_name)

    else:
        dir = glob.glob(path + "NACA*.dat")
        for f_n in dir:
            if (f_n[p_n+8] == "_"): # 4digit
                new_name = f_n[:p_n + 16] + ".dat"
                rename_on_win(f_n, new_name)
            else:
                new_name = f_n[:p_n + 17] + ".dat"
                rename_on_win(f_n, new_name)


if __name__ == '__main__':
    # main()
    # rc('text', usetex=True)
    # cp_plot_test()
    # cdcl_plot_test()
    # plot_residual_graph()
    make_restart_list(digit=4, deg_list=[0,3,6,9,12,15,18,21,24,27,30,33,36,39])
    # make_restart_list(digit = 5, deg_list = [0,3,6,9,12,15,18,21,24,27,30,33,36,39])
    # del_duplication(cp = True)