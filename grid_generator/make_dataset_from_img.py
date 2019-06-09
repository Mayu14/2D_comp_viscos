# -- coding: utf-8 --
import numpy as np
import cv2
import glob

def main(fname_list, img_path, csv_name):
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
        naca = float(param[0].replace("NACA", ""))
        aspect = float(param[1])
        return naca, aspect

    def read_img_and_drop_channel(fname):
        return cv2.imread(fname)[:,:,0]

    total_data = len(fname_list)
    data_length = read_img_and_drop_channel(fname_list[0]).reshape(-1).shape[0]
    save_data = np.zeros((total_data, data_length + 4), dtype=float) # naca_number + shape[0] + shape[1] + aspect = 4

    for i in range(total_data):
        fname = fname_list[i]
        naca4, aspect_ratio = get_number_and_aspect(fname.replace(img_path, ""))
        img = read_img_and_drop_channel(fname)
        shape = img.shape
        # NACA 4or5 digit, x_width, y_height, aspect_ratio, rgb_pixel(only "r" component)
        save_data[i, 0] = naca4
        save_data[i, 1] = shape[0]
        save_data[i, 2] = shape[1]
        save_data[i, 3] = aspect_ratio
        save_data[i, 4:] = img.reshape(-1)

    np.savetxt(csv_name, save_data, delimiter=",")

def auto_job(type=4, size=512):
    save_path = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\"

    if type == 4:
        dir = "NACA4\\"
    elif type == 5:
        dir = "NACA5\\"

    csv_name = save_path + dir + "shape_img_all.csv"
    img_path = save_path + dir + "cnn\\image\\" + str(size) + "\\"
    flist = glob.glob(img_path + "*.png")

    main(flist, img_path, csv_name=csv_name)


def test():
    fname_list = ["NACA62132__7.401631137049297__.png"]
    csv_name = "test.csv"
    main(fname_list, "", csv_name)

if __name__ == '__main__':
    auto_job(type=5)
