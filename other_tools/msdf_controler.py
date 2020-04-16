# coding: utf-8
from other_tools.setup_foam import load_patternList
from other_tools.generate_runOF import gen_casename
from training.read_training_data import load_df_without_shape
from naca_4digit_test import Naca_4_digit, Naca_5_digit
import numpy as np
import os
from pathlib import Path
import cv2
import glob
import pandas as pd

def test():
    msdf_reses = ["256", "512", "1024", "2048"]
    naca = Naca_4_digit("2611", attack_angle_deg=4, resolution=1000, quasi_equidistant=False)
    for msdf_res in msdf_reses:
        fname = str(msdf_res) + "_" + "c_i_" + str(4) + "_NACA" + str(2611) + "_aoa" + str(4).zfill(3) + ".svg"
        gen_cmd(naca, fname, msdf_res, msdf_res, debug=True)

def gen_cmd(wing, svg_name, size_x, size_y, dir="D:\\Toyota\\Downloads\\msdfgen-win64\\msdfgen\\NACA", debug=False):
    if not debug:
        dir += str(size_x) + "\\"
    else:
        dir += "test\\"
    png_name = svg_name[:-4]
    png_name += ".png"
    if not os.path.exists(dir + svg_name):
        wing.generate_svg(dir + svg_name)
    cmd = "G: && cd " + dir + "&& "
    cmd += 'msdfgen.exe -svg ' + svg_name + ' -o ' + png_name + ' -size ' + str(size_x) + ' ' + str(size_y) + ' -autoframe'
    os.system(cmd)

def naca_msdf_generator(number = "0012", aoa=0.0, fileHeader = "", msdf_res=16, dir=""):
    type = len(number)
    if type == 4:
        wing_maker = Naca_4_digit
        dir = "G:\\Toyota\\Data\\grid_vtk\\NACA4\\" + dir
    elif type == 5:
        wing_maker = Naca_5_digit
        dir = "G:\\Toyota\\Data\\grid_vtk\\NACA5\\" + dir
    naca = wing_maker(number, attack_angle_deg = aoa, resolution = 1000, quasi_equidistant = False)
    fname = str(msdf_res) + "_" + fileHeader + str(type) + "_NACA" + number + "_aoa" + str(aoa).zfill(3) + ".svg"
    gen_cmd(naca, fname, msdf_res, msdf_res, dir=dir)

def main_iter(msdf_res = 16, case = "CV_mixed"):
    naca4list = [True, False]
    for naca4 in naca4list:
        reynoldsList, machList, angleList, i1_list, i2_list, i3_list, i2_len = load_patternList(naca4, case = case)
        for i1 in i1_list:
            for i2 in i2_list:
                if ((i1 == 0) and (i2 == 0)) or (i1 != 0):
                    for i34 in i3_list:
                        int_4 = str(i1) + str(i2).zfill(i2_len) + str(i34)
                        print(int_4)
                        for angle in angleList:
                            naca_msdf_generator(number = int_4, aoa = angle, fileHeader = "c_v_", msdf_res = msdf_res, dir = case + "\\")


def incomp_invicsid_case(naca4=True, msdf_res = 16):
    fileHeader = "i_i_" #incomp inviscid
    if naca4:
        for i in range(5000):
            wing_number = str(2 * i + 1).zfill(4)
            print(wing_number)
            for angle in range(40):
                deg = - 40 + 2 * angle
                naca_msdf_generator(number = wing_number, aoa = deg, fileHeader = fileHeader, msdf_res = msdf_res)
    
    else:
        head_int3 = [210, 220, 230, 240, 250, 221, 231, 241, 251]
        for int3 in head_int3:
            for int2 in range(1, 100):
                wing_number = str(int3) + str(int2).zfill(2)
                print(wing_number)
                for angle in range(40):
                    deg = - 40 + 2 * angle
                    naca_msdf_generator(number = wing_number, aoa = deg, fileHeader = fileHeader, msdf_res = msdf_res)

def msdf_pack(naca4=True, msdf_res=16, case=""):

    def main(fname_list, msdf_res=msdf_res, type = 4, case=case):
        """

        :param fname_list:  画像ファイルのファイル名リスト
        :param img_path:  画像ファイルの入っているディレクトリフルパス
        :param csv_name: 保存するcsvのフルパス
        :return:
        """
    
        def get_number_and_aoa(fname):
            """

            :param fname:  ex:"16_i_i_4_NACA0001_aoa-02.png"
            :return:
            """
            param = fname.split("_")
            naca = int(param[4].replace("NACA", ""))
            aoa = float(param[5].replace("aoa", "").replace(".png", ""))
            return naca, aoa
    
        def read_img(fname):
            return cv2.imread(fname)
    
        def merge2lift_data(df_l, save_data_i, save_data_f, source, size = 512, case=""):
            if case == "":
                name_s = ["naca4"]
                name_s2 = ["angle"]
            elif case == "CV_mixed":
                name_s = ["NACA"]
                name_s2 = ["AoA"]
            # common
            name_s.extend([str(i).zfill(3) for i in range((size ** 2)*3)]) # x * r * RGB

            df_s1 = pd.DataFrame(save_data_i, columns = name_s, dtype = "int32")
            del save_data_i

            df_s2 = pd.DataFrame(save_data_f, columns = name_s2, dtype = "float32")
            del save_data_f

            df_s1 = df_s1.join(df_s2)
            del df_s2

            df_s1 = pd.merge(df_l, df_s1)
            del df_l

            if not df_s1.empty:
                if case == "":
                    y_train = df_s1["lift_coef"].values
                elif case == "CV_mixed":
                    y_train = df_s1[["Cl", "Cd"]].values
                
                data_num = len(df_s1)
                if case == "":
                    drop_col = ["naca4", "angle", "lift_coef"]
                elif case == "CV_mixed":
                    drop_col = ["NACA", "AoA", "Ma", "Re", "Cl", "Cd"]
                    x_train_label = df_s1[["Ma", "Re"]].values

                x_train_img = df_s1.drop(drop_col, axis = 1).values.reshape(data_num, size, size, 3)
                del df_s1
            
                # strstep = str(step).zfill(3)
            
                if type == 4:
                    dir = "NACA4\\"
                    footer = "4"
                # elif type == 5:
                else:
                    dir = "NACA5\\"
                    footer = "5"
                    
                if case == "":
                    code_name = "_i_i_"
                elif case == "CV_mixed":
                    code_name = "_c_v_"
                shape = str(size) + code_name
            
                if case == "":
                    np.savez_compressed(source + dir + shape + footer, x_train = x_train_img, y_train = y_train)
                elif case == "CV_mixed":
                    np.savez_compressed(source + dir + shape + footer, x_img = x_train_img, x_label = x_train_label, y_train = y_train)
    
        total_data = len(fname_list)

        data_length = read_img(str(fname_list[0])).reshape(-1).shape[0]
        if case == "":
            source = "G:\\Toyota\\Data\\Incompressible_Invicid\\training_data\\"
            if type == 4:
                fpath_lift = "NACA4\\s0000_e5000_a040_odd.csv"
            elif type == 5:
                fpath_lift = "NACA5\\s21001_e25199_a040_(angle_-40_38).csv"
        
            name = ("naca4", "angle", "lift_coef")
            data_type = {"naca4": int, "angle": float, "lift_coef": float}
    
            df_l = pd.read_csv(source + fpath_lift, header = None,
                               usecols = [1, 3, 4], names = name, dtype = data_type)
            
        elif case == "CV_mixed":
            source = "G:\\Toyota\\Data\\Compressible_Viscos\\rhoSimpleFoam\\training_data\\input_vectors\\"
            df_l = load_df_without_shape(concat = True, keyLists = ["Cl", "Cd"])

        else:
            raise ValueError
        
        save_data_i = np.zeros((total_data, 1 + data_length),
                               dtype = int)  # naca_number + imgRGB
        save_data_f = np.zeros((total_data), dtype = float)

        for i in range(total_data):
            fname = fname_list[i]
            
            naca4, aoa = get_number_and_aoa(str(fname.name))
            img = read_img(str(fname))
            shape = img.shape
            # NACA 4or5 digit, x_width, y_height, aspect_ratio, rgb_pixel(only "r" component)
            save_data_i[i, 0] = naca4
            save_data_i[i, 1:] = img.reshape(-1)
            save_data_f[i] = aoa
            
        merge2lift_data(df_l, save_data_i, save_data_f, source, size = msdf_res, case=case)
    
        # np.savetxt(csv_name, save_data, delimiter=",")

    if naca4:
        type = 4
    else:
        type = 5

    folder = Path("NACA" + str(type))
    dir = Path("G:\\Toyota\\Data\\grid_vtk")
    if case == "":
        dir = dir / folder / Path(str(msdf_res))
    elif case == "CV_mixed":
        dir = dir / folder / Path("CV_mixed\\" + str(msdf_res))
    else:
        raise ValueError
    
    flist = list(dir.glob(str(msdf_res) + "?????" + str(type) + "*.png"))

    main(flist, msdf_res, type = type, case = case)
    

if __name__ == '__main__':
    # main_iter(msdf_res=32)
    # main_iter(msdf_res = 64)
    # main_iter(msdf_res = 128)
    # exit()
    # for size in range(4):
    # test()
    # exit()
    size = 3
    msdf_res = 2 ** (size + 4)
    print(msdf_res)
    # incomp_invicsid_case(msdf_res = msdf_res)
    # incomp_invicsid_case(naca4 = False, msdf_res = msdf_res)
    # exit()
    msdf_pack(naca4 = True, msdf_res = msdf_res, case = "CV_mixed")
    