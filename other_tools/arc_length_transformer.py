# coding: utf-8
from scipy.interpolate import PchipInterpolator
from scipy.integrate import trapz
import numpy as np
from numba import jit
from matplotlib import pyplot as plt
from naca_4digit_test import Naca_4_digit, Naca_5_digit

class ArcLength():
    def __init__(self, x=None, y=None, normalize=1, min_dif =10e-6,
                 vector_dimension=200, base_resolution=5000, encoder_mode=True,
                 encoded_obj=None, overwrite_bn=True):
        # common
        self.machine_epsilon = min_dif
        self.normalize = normalize  # 曲線を[0,1]にする場合1, [0, pi]にする場合np.pi
        self.n_max = int(vector_dimension / 2)
        self.base_resolution = base_resolution
        self.overwrite_bn = overwrite_bn
        
        if encoder_mode:
            if (type(x) == type(None)) or (type(y) == type(None)):
                print("encoder mode requires x and y inputs")
                exit()

            self.check_array(x, y)

            self.get_arc_length_param()
            self.set_arc_vs_xy()
            self.fourier_sin_serires()
            self.encode_obj_vector()
        else:
            if (type(encoded_obj) == type(None)):
                print("The decoder mode requires input of encoded_obj")
            self.decode_obj_vector(encoded_obj)
            self.decrypt(base_resolution)
        #self.plot_s_xy_plane(origin=False)

    def check_array(self, x, y):
        size = x.shape[0]
        size_y = y.shape[0]
        if size != size_y:
            print("input arrays x and y have difference shape")
            exit()
        
        if min(abs(x[size - 1] - x[0]), abs(y[size_y - 1] - y[0])) > self.machine_epsilon:
            x = np.concatenate([x, [x[0]]])
            y = np.concatenate([y, [y[0]]])

        self.x = x
        self.y = y
        self.size = x.shape[0]
    
    def get_arc_length_param(self):
        start = np.zeros(1)
        dif_x = np.concatenate([start, np.diff(self.x)])    # 座標の差分
        dif_y = np.concatenate([start, np.diff(self.y)])
        dist = np.sqrt(dif_x**2 + dif_y**2) # 三平方の定理
        dist = dist.cumsum()    # 部分和
        self.length = dist[-1]  # 曲線の長さ
        self.arcLength = dist / self.length * self.normalize # 長さをnormalizeの値に調整
        
    def set_arc_vs_xy(self):
        self.fx = PchipInterpolator(self.arcLength, self.x)
        self.fy = PchipInterpolator(self.arcLength, self.y)

    #@jit
    def fourier_sin_serires(self):
        L = self.normalize
        n_max = self.n_max
        self.bn_x = np.zeros(n_max)
        self.bn_y = np.zeros(n_max)
        s = np.linspace(start=0, stop=1, num=self.base_resolution)
        x = self.fx(s)
        y = self.fy(s)
        for n in range(n_max):
            x_star = 2.0 * x * np.sin((n + 1) * np.pi * s / L)
            y_star = 2.0 * y * np.sin((n + 1) * np.pi * s / L)
            self.bn_x[n] = trapz(y=x_star, x=s, axis=-1)
            self.bn_y[n] = trapz(y=y_star, x=s, axis=-1)

    def decrypt(self, resolution=100):
        decryption_x = np.zeros(resolution)
        decryption_y = np.zeros(resolution)
        s = np.linspace(start=0, stop=1, num=resolution)
        L = self.normalize
        n_max = self.n_max

        for n in range(n_max):
            decryption_x += self.bn_x[n] * np.sin((n + 1) * np.pi * s / L)
            decryption_y += self.bn_y[n] * np.sin((n + 1) * np.pi * s / L)

        self.decryption_s = s
        self.decryption_x = decryption_x
        self.decryption_y = decryption_y
        self.decryption = np.concatenate([decryption_x, decryption_y]).reshape(2, -1)

    def plot_s_xy_plane(self, origin=True):
        s_rough = np.linspace(start=0, stop=1, num=20)
        s_fine = np.linspace(start=0, stop=1, num=2000)
        if origin:
            plt.plot(self.arcLength, self.x)
            plt.plot(self.arcLength, self.fx(self.arcLength))
            plt.plot(self.arcLength, self.y)
            plt.plot(self.arcLength, self.fy(self.arcLength))

            fig = plt.figure(figsize=[16, 3])
            ax = fig.add_subplot(141)
            ax.plot(self.x, self.y, label="original")
            ax.set_title("NACA6189")
            ax.set_xlabel("$x/c$")
            ax.set_ylabel("$y/c$")

            ax = fig.add_subplot(142)
            ax.plot(s_fine, self.fx(s_fine), label="$x=x(t)$")
            ax.set_title("$x=x(t)$")
            ax.set_xlabel("arc-length")
            ax.set_ylabel("$x/c$")

            ax = fig.add_subplot(143)
            ax.plot(s_fine, self.fy(s_fine), label="$y=y(t)$")
            ax.set_title("$y=y(t)$")
            ax.set_xlabel("arc-length")
            ax.set_ylabel("$y/c$")

            self.decrypt()
            ax = fig.add_subplot(144)
            ax.plot(self.decryption_x, self.decryption_y, label="decoded")
            ax.set_title("Decoded NACA6189")
            ax.set_xlabel("$x/c$")
            ax.set_ylabel("$y/c$")
            ax.legend()
            plt.tight_layout()
            plt.show()

            plt.plot()
            plt.plot(s_fine, self.fy(s_fine))
            plt.show()
        else:
            """
            plt.plot(self.fx(self.arcLength), self.fy(self.arcLength))
            plt.plot(self.x, self.y, "x")
            plt.show()
            """
            self.decrypt()
            plt.plot(self.decryption[0], self.decryption[1])
            plt.plot(self.x, self.y)
            plt.show()

    def encode_obj_vector(self):
        self.obj = np.concatenate([self.bn_x, self.bn_y])

    def decode_obj_vector(self, obj_vector):
        size = obj_vector.shape[0]
        self.n_max = int(size/2)
        self.decoded_bn_x = obj_vector[:self.n_max]
        self.decoded_bn_y = obj_vector[self.n_max:]
        if self.overwrite_bn:
            self.bn_x = self.decoded_bn_x
            self.bn_y = self.decoded_bn_y
            
    def fixed_plot(self, x, y, imgWidth=200, imgHeight=200, xLable="", yLabel="", dataLabel="", title="", legend=False, no_axis=True, return_array=False):
        
        if return_array:
            fig = plt.figure(figsize = (max(int(imgWidth/100), 1), max(int(imgHeight/100), 1)))
        else:
            fig = plt.figure(figsize = (4, 4))
            
        ax = fig.add_subplot(111)
        ax.plot(x, y, label=dataLabel)
        ax.set_title(title)
        ax.set_xlabel(xlabel = xLable)
        ax.set_ylabel(ylabel = yLabel)
        plt.axis("equal")
        if no_axis:
            plt.axis("off")
        if legend:
            ax.legend()
        # plt.tight_layout()
        if not return_array:
            plt.show()
        else:
            fig.canvas.draw()
            img = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
            img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))
            plt.close(fig)
            return img


def make_mod_fourier_from_naca(naca4, resolution, digit=4):
    if digit == 4:
        naca = Naca_4_digit(naca4, attack_angle_deg=0, resolution=resolution, quasi_equidistant=False, position="start_0")
    elif digit == 5:
        naca = Naca_5_digit(naca4, attack_angle_deg=0, resolution=resolution, quasi_equidistant=False, position="start_0")
    else:
        print("line130: make_mod_fourier_from_naca: please check digit")
    naca.one_stroke_ccw()
    arc = ArcLength(x=naca.ccw_x, y=naca.ccw_y)
    return arc.obj

def check_decode_input(bnxy, from_concat_bn_xy, bnx, bny, from_bn_x_and_bn_y):
    if (from_bn_x_and_bn_y) and (from_concat_bn_xy):
        raise ValueError  # choose 1 mode only

    if from_bn_x_and_bn_y:
        if (bnx is None) or (bny is None):
            raise ValueError  # set bnx and bny
        bnxy = np.concatenate([bnx, bny])

    elif from_concat_bn_xy:
        if bnxy is None:
            raise ValueError  # set bnxy
    else:
        raise ValueError  # set decode mode

    return bnxy

def rotate(t, decodedAL):
    xy = np.concatenate([decodedAL.decryption_x.reshape(-1, 1), decodedAL.decryption_y.reshape(-1, 1)], axis = 1).T
    gc = np.average(xy, axis = 1).reshape(-1, 1)
    rot = [[np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)]]
    xy = (np.dot(rot, xy - gc) + gc).T
    return xy[:, 0], xy[:, 1]

def get_obj_rotated(bnxy, angle, deg):
    decode = ArcLength(encoder_mode = False, encoded_obj = bnxy)
    if deg:
        t = angle * np.pi / 180.0
    else:
        t = angle
    decode.decryption_x, decode.decryption_y = rotate(t, decode)
    return decode

def get_obj_cplx(bnxy, angle, deg, epsilon = 0.001):
    decode = get_obj_rotated(bnxy, angle, deg)
    z = decode.decryption_x + 1j * decode.decryption_y
    dist1 = (decode.decryption_x[0] - decode.decryption_x[-1]) ** 2 + (decode.decryption_y[0] - decode.decryption_y[-1])**2
    dist2 = (decode.decryption_x[-2] - decode.decryption_x[-1]) ** 2 + (decode.decryption_y[-2] - decode.decryption_y[-1])**2

    if dist1 < epsilon * dist2:
        return z[1:]
    else:
        return z

def get_decode_img_from_bnxy(bnxy=None, angle=0.0, from_concat_bn_xy=False, bnx=None, bny=None, from_bn_x_and_bn_y=False, imgWidth=200, imgHeight=200, deg=False):
    bnxy = check_decode_input(bnxy, from_concat_bn_xy, bnx, bny, from_bn_x_and_bn_y)
    decode = get_obj_rotated(bnxy, angle, deg)
    return decode.fixed_plot(decode.decryption_x, decode.decryption_y, return_array = True, imgWidth=imgWidth, imgHeight=imgHeight)

def decode_into_cplx(bnxy=None, angle=0.0, from_concat_bn_xy=False, bnx=None, bny=None, from_bn_x_and_bn_y=False, deg=False):
    bnxy = check_decode_input(bnxy, from_concat_bn_xy, bnx, bny, from_bn_x_and_bn_y)
    z = get_obj_cplx(bnxy, angle, deg)
    return z

def main(x,y):
    # 始点 != 終点のとき，点を追加する

    # 長さの差分求める
    
    pass

if __name__ == '__main__':
    naca = Naca_4_digit("6189", attack_angle_deg = 0, resolution = 200, quasi_equidistant = False, position = "start_0")
    naca.one_stroke_ccw()
    x = naca.ccw_x
    y = naca.ccw_y

    arc = ArcLength(x,y, vector_dimension=200, base_resolution=5000)
    # arc.plot_s_xy_plane(origin=False)
    # arc.plot_s_xy_plane(origin=True)
    
    obj = arc.obj
    decode = ArcLength(encoder_mode = False, encoded_obj = obj)
    # plt.plot(decode.decryption_x, decode.decryption_y)
    # plt.show()
    # z = decode_into_cplx(bnxy=obj, from_concat_bn_xy=True, angle=90.0, deg=True)
    # plt.plot(np.real(z), np.imag(z))
    # plt.show()
    # exit()
    # decode.fixed_plot(decode.decryption_x, decode.decryption_y, return_array = False)
    img = decode.fixed_plot(decode.decryption_x, decode.decryption_y, return_array = True)
    
    img = get_decode_img_from_bnxy(bnxy=obj, from_concat_bn_xy = True, angle = 90.0)
    print(img.shape)
    import cv2
    from PIL import Image
    pilImg = Image.fromarray(np.uint8(img)[:, :, :])  # RBG -> RGB
    # pilImg.save("sample.png")
    # exit()
    cv2.imshow("img", img)
    cv2.waitKey()
    #plt.plot(x, y)
    #plt.show()
