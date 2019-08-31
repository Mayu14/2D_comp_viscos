# coding: utf-8
import numpy as np
from other_tools.legacy_vtk_writer import MyStructuredVTKWriter
from grid_generator.naca_4digit_test import Naca_4_digit as N4d  # use only test
from grid_generator.naca_4digit_test import Naca_5_digit as N5d  # use only test
from PIL import Image
from grid_generator.mask_maker import VPmask as MM
import skfmm
import cv2

# 標準物体形状は[0, 1]×[0, 1]領域に描画されるものとする
# 即ち，広義の物体形状をmax(xmax - xmin, ymax - ymin)にて規格化したのち，原点が(0.5, 0.5)となるよう平行移動したものを標準物体形状と定める
# Shapeオブジェクトは，(上面関数，下面関数，親格子解像度，マージンセル数)によって初期化される

class Shape(object):
    def __init__(self, y_upper, y_lower, resolution, aspect_ratio=1.0, resize=True, name="Shape", auto_reshape=True,
                 path="", figformat="png"):
        """
        
        :param y_upper: (ndarray) y上側の形状
        :param y_lower: (ndarray) y下側の形状
        :param resolution: (integer) 画像の解像度
        :param aspect_ratio: (グリッドのアスペクト比) 格子生成する場合は解像度を上げたい方向に応じて変更，画像生成の場合は1.0で固定
        :param resize: (binary) 格子生成時に，物体長軸方向両端に隙間を開けるかどうか
        :param name: (character)  生成される画像の名前(header), resize=Trueで拡大の比率がfooterとして付加される
        :param auto_reshape: (binary) 画像生成時に，形状を画像の幅いっぱいに拡大するかどうか
        :param path: (character)    画像の保存先ディレクトリ
        :param figformat:  画像の保存形式
        """
        self.upper = y_upper  # ndarray
        self.lower = y_lower  # ndarray
        self.aspect = np.array([1.0, aspect_ratio])
        self.resolution = resolution
        self.grid_dx = 1.0 / resolution  # [0, 1]区間をresolution分割    # 厳密には(size + 1)であるが，形状誤差の許容値dx/sizeを与えて簡略化することで境界における問題を無くす
        self.grid_dy = self.grid_dx * aspect_ratio
        self.resize = resize
        self.name = name
        self.path = path
        self.figformat = figformat
        self.auto_reshape = auto_reshape
        self.get_mask()

    def change_aspect(self):
        eps = self.grid_dy
        y_range = np.max(self.upper) - np.min(self.lower)   # yの範囲を求めて
        y_center = (np.max(self.upper) + np.min(self.lower)) / 2.0  # yの中心を求めて
        ratio =(1.0 - 2.0 * eps) / y_range  # y方向の拡大率を求めて   # y方向上下端にdyの隙間を作成
        self.upper = (self.upper - y_center) * ratio + y_center  # 中心を0に合わせてから拡大，中心位置を元に戻す
        self.lower = (self.lower - y_center) * ratio + y_center  # 中心を0に合わせてから拡大，中心位置を元に戻す
        self.magnification = ratio

    def y_shift(self):
        center = (np.max(self.upper) + np.min(self.lower)) / 2.0 - 0.5
        self.lower -= center
        self.upper -= center

    def get_mask(self):
        size = self.resolution  # 格子の分割数
        # """
        if self.resize == True:
            y_u = self.upper
            y_l = self.lower
            self.upper = np.zeros(size)
            self.lower = np.zeros(size)
            # self.upper[int(size/4):int(size/4*3)] = 0.5 * (y_u[0::2] - 0.5) + 0.5 * self.aspect[1]
            # self.lower[int(size / 4):int(size / 4 * 3)] = 0.5 * (y_l[0::2] - 0.5) + 0.5 * self.aspect[1]
            self.upper[int(size / 8 * 3):int(size / 8 * 5)] = 0.25 * (y_u[0::4] - 0.5) + 0.5 * self.aspect[1]
            self.lower[int(size / 8 * 3):int(size / 8 * 5)] = 0.25 * (y_l[0::4] - 0.5) + 0.5 * self.aspect[1]
        # """
        if self.auto_reshape:
            self.change_aspect()
        self.y_shift()
        mask = np.full((size, size), 1.0)  # マスク関数1:流体，0:物体
        # print(mask)

        obj_resolution = size  # 物体形状の解像度
        obj_x = np.linspace(start=0, stop=1, num=obj_resolution)  # margin分を除去した数列
        obj_dx = 1.0 / (obj_resolution + 1)
        y_u_by_dy = np.floor(self.upper / self.grid_dy)  # y上側点は格子何個分なのか
        y_l_by_dy = np.floor(self.lower / self.grid_dy)  # y下側点は格子何個分なのか
        y_width = self.upper - self.lower  # y上側点とy下側点との間に格子は何個あるのか

        grid_x = np.linspace(start=0, stop=1, num=size)  # 元格子の格子列

        obj_cell = 0  # 物体形状の格子確認用
        epsilon = 0.01 * self.grid_dy  # 上下の点の間の幅がこの値以下なら格子の存在を認めない
        for grid_cell in range(size):
            if ((grid_cell >= 0.0) or ((grid_cell < (size - 1) - 0.0))):
                if (y_width[obj_cell] > epsilon):
                    mask[grid_cell, int(y_l_by_dy[obj_cell]):int(y_u_by_dy[obj_cell] + 1)] = 0
                obj_cell += 1

            if (obj_cell == y_width.shape[0]):
                break

        self.mask = np.flipud(mask)[-1::-1]

    def test_output(self):
        x = np.linspace(start=0, stop=1, num=self.mask.shape[0] + 1)
        z = [0]
        writer = MyStructuredVTKWriter(xcoords=x, ycoords=x, zcoords=z, data_dimension=2)
        writer.filename = "mask"
        writer.casename = "mask_naca" + self.name
        writer.append_data(name="object", SCALARS=True, data=self.mask)
        writer.output("mask")

    def mask2pict(self):
        def set_color(white=True):
            if white:
                return 255, 255, 255
            else:
                return 0, 0, 0

        size = self.mask.shape
        im = Image.new("RGB", size)
        for x in range(size[0]):
            for y in range(size[1]):
                r, g, b = set_color((self.mask[x, y] == 0))
                im.putpixel((x, y), (r, g, b))

        im_name = self.path + self.name
        if self.auto_reshape:
            im_name += "__" + str(self.magnification) + "__"
        im_name += "." + self.figformat
        im.save(im_name, format=self.figformat)


def y_upper_test(x, kind="real"):
    if kind == "real":
        return np.where(x < 0.5, 0.5, x)
    elif kind == "int":
        y = np.linspace(start=0, stop=1, num=100)
        return y[int(x)]


def y_upper(x, kind="real"):
    if kind == "real":
        return np.where(x < 0.5, 0.5, x)
    elif kind == "int":
        y = np.linspace(start=0, stop=1, num=100)
        return y[int(x)]


def y_lower(x):
    return np.where(x < 0.5, 0.5, 0.5)


def main():
    naca_4digit_name = "0003"
    attack_angle_deg = 0
    np.set_printoptions(threshold=np.inf)
    grid_resolution = 2 ** 9
    naca = N4d(naca_4digit_name, attack_angle_deg, grid_resolution)

    aspectratio = 1.0  # y/x
    resize = False  # int(grid_resolution / 8)
    s = Shape(naca.equidistant_y_u, naca.equidistant_y_l, grid_resolution, aspectratio, resize, naca_4digit_name)
    s.mask2pict()
    exit()
    s.test_output()

    """
    x = np.linspace(start = 0, stop = 1, num = 100)
    y_u = s.upper(x)
    y_l = s.lower(x)

    print(y_u)
    print(y_l)
    """

def make_dataset(type=4, grid_resolution = 2 ** 9, sdf=False):
    def concat_u_and_l(naca):
        x = np.concatenate([naca.x_u, naca.x_l[::-1]])
        y = np.concatenate([naca.y_u, naca.y_l[::-1]])
        return [x], [y]

    def main_process(naca4, path, grid_resolution, NACA=N4d, sdf=False):
        attack_angle_deg = 0
        naca = NACA(naca4, attack_angle_deg, grid_resolution)
        x, y = concat_u_and_l(naca)
        """
        s = Shape(naca.equidistant_y_u, naca.equidistant_y_l, grid_resolution, resize = False,
                  name = "NACA" + naca4, auto_reshape = True, path = path)
        s.mask2pict()
        """
        canvas = [grid_resolution, grid_resolution]
        fname = path + naca4

        if sdf:
            imglim = [[-0.1, 1.1], [-0.1, 1.1]]
            mask = MM(x, y, canvas, default_canvas = canvas, deform = "None", start_point = "rb", imglim = imglim)
            mask.save_sdf(fname, save_contour = True)
            
        else:
            mask = MM(x, y, canvas, deform="Fit", start_point="rt")
            mask.save_img(fname)
        
    
    path = "G:\\Toyota\\Data\\Compressible_Invicid\\training_data\\"
    if type == 4:
        path += "NACA4\\cnn\\image\\" + str(grid_resolution) + "\\"
        for i1 in range(0,10):
            for i2 in range(0,10):
                for i34 in range(1, 100):
                    naca4 = str(i1) + str(i2) + str(i34).zfill(2)
                    print(naca4)
                    main_process(naca4, path, grid_resolution, sdf = sdf)
                    
                    
                        
    elif type == 5:
        path += "NACA5\\cnn\\image\\" + str(grid_resolution) + "\\"
        mid_int2 = [10, 20, 30, 40, 50, 21, 31, 41, 51]
        for int1 in range(1, 6):
            for int23 in mid_int2:
                for i45 in range(1, 100):
                    naca5 = str(int1) + str(int23) + str(i45).zfill(2)
                    main_process(naca5, path, grid_resolution, NACA=N5d, sdf = sdf)
    
    
if __name__ == '__main__':
    # main()
    # for j in range(9, 10):
        # for i in range(4,6):
    sdf = True
    for i in range(5, 6):
        make_dataset(type = i, grid_resolution=2**10, sdf=sdf)
        