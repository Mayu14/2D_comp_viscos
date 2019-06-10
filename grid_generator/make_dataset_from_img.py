# -- coding: utf-8 --
import numpy as np
import cv2
import glob
import pandas as pd

def main(fname_list, img_path, size):
    """

    :param fname_list:  画像ファイルのファイル名リスト
    :param img_path:  画像ファイルの入っているディレクトリフルパス
    :param csv_name: 保存するcsvのフルパス
    :return:
    """
    def get_number_and_aspect(fname):
        """

        :param fname:  ex:"NACA62132__7.401631137049297__.png"
        :return:
        """
        param = fname.split("__")
        naca = int(param[0].replace("NACA", ""))
        aspect = float(param[1])
        return naca, aspect

    def read_img_and_drop_channel(fname):
        return cv2.imread(fname)[:,:,0]
    
    def merge2lift_data(df_l, save_data_i, save_data_f, source, size = 512):
        name_s = ["naca4", "shape0", "shape1"]
        name_s.extend([str(i).zfill(3) for i in range(size**2)])

        df_s1 = pd.DataFrame(save_data_i, columns = name_s, dtype = "int16")
        del save_data_i

        name_s2 = ["aspect"]
        df_s2 = pd.DataFrame(save_data_f, columns = name_s2, dtype = "float32")
        del save_data_f
        
        df_s1 = df_s1.join(df_s2)
        del df_s2
        # print(df_s1)
        # print(df_s2)
        df_s1 = pd.merge(df_l, df_s1)
        del df_l
        if not df_s1.empty:
            y_train = np.concatenate(
                (df_s1["lift_coef"].values.reshape(-1, 1).T, df_s1["drag_coef"].values.reshape(-1, 1).T)).T
    
            # df_s1 = df_s1.drop("lift_coef", axis = 1).drop("drag_coef", axis = 1).drop("naca4", axis=1).values
            x_train_param = np.concatenate(
                (df_s1["angle"].values.reshape(-1, 1).T, df_s1["aspect"].values.reshape(-1, 1).T)).T   # need to add Mach, Reynolds
            
            shape0 = df_s1["shape0"][0]
            shape1 = df_s1["shape1"][0]
            data_num = len(df_s1)
            drop_col = ["naca4", "angle", "lift_coef", "drag_coef", "shape0", "shape1", "aspect"]
            x_train_img = df_s1.drop(drop_col, axis=1).values.reshape(data_num, shape0, shape1)
            del df_s1
            
            # strstep = str(step).zfill(3)
            dir = "NACA4\\"
            shape = "_" + str(shape0).zfill(3) + "_" + str(shape1).zfill(3)
            """
            fnameXimg = "x_train_img"
            fnameXprm = "x_train_param"
            fnamey = "y_train"
            
            np.save(source + dir + fnameXimg + strstep, x_train_img)
            np.save(source + dir + fnameXprm + strstep, x_train_param)
            np.save(source + dir + fnamey + strstep, y_train)
            """
            np.savez_compressed(source + dir + "train" + shape, x_train_img=x_train_img, x_train_param=x_train_param, y_train=y_train)
            
        
    total_data = len(fname_list)
    
    data_length = read_img_and_drop_channel(fname_list[0]).reshape(-1).shape[0]

    source = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\"
    fpath_lift = "NACA4\\s1122_e9988_s4_a014.csv"
    
    name = ("naca4", "angle", "lift_coef", "drag_coef")
    data_type = {"naca4": int, "angle": float, "lift_coef": float, "drag_coef": float}

    df_l = pd.read_csv(source + fpath_lift, header = None,
                       usecols = [1, 3, 4, 5], names = name, dtype = data_type)

    
    save_data_i = np.zeros((total_data, data_length + 3), dtype=int) # naca_number + shape[0] + shape[1] + aspect = 4
    save_data_f = np.zeros((total_data, 1), dtype=float)  # only aspect = 1

    for i in range(total_data):
        fname = fname_list[i]
        
        naca4, aspect_ratio = get_number_and_aspect(fname.replace(img_path, ""))
        img = read_img_and_drop_channel(fname)
        shape = img.shape
        # NACA 4or5 digit, x_width, y_height, aspect_ratio, rgb_pixel(only "r" component)
        save_data_i[i, 0] = naca4
        save_data_i[i, 1] = shape[0]
        save_data_i[i, 2] = shape[1]
        save_data_i[i, 3:] = img.reshape(-1)
        save_data_f[i, 0] = aspect_ratio
        
    merge2lift_data(df_l, save_data_i, save_data_f, source, size=size)
        
    # np.savetxt(csv_name, save_data, delimiter=",")

def auto_job(type=4, size=512):
    save_path = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\"

    if type == 4:
        dir = "NACA4\\"
    elif type == 5:
        dir = "NACA5\\"

    #csv_name = save_path + dir + "shape_img_all.csv"
    img_path = save_path + dir + "cnn\\image\\" + str(size) + "\\"
    flist = glob.glob(img_path + "*.png")

    main(flist, img_path, size=size)


def test():
    fname_list = ["NACA62132__7.401631137049297__.png"]
    csv_name = "test.csv"
    main(fname_list, "", csv_name)

if __name__ == '__main__':
    auto_job(type=4, size=512)
