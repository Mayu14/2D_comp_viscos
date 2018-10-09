# -- coding: utf-8 --
import numpy as np
import matplotlib.pyplot as plt
from naca_4digit_test import Naca_4_digit
def re_numbering(z):
    start = np.argmax(np.real(z))
    return np.concatenate([z[start:z.shape[0] - 1], z[:start + 1]])

def joukowski_tr(z, a):
    return z + a**2/z

def karman_trefftz_tr(z, b, n=1.94):
    return n * b * ((z + b)**n + (z - b)**n) / ((z + b)**n - (z - b)**n)

def circle_center_a(size, center_x, center_y):
    R = np.sqrt((1.0 - center_x) ** 2 + center_y ** 2)
    b = center_x + 1j * center_y
    t = np.linspace(start = 0, stop = 2.0*np.pi, num = size + 1)
    return R * np.exp(1j * t) + b

def validation(size, center_x, center_y):
    z1 = circle_center_a(size, center_x, center_y)
    z2 = circle_center_a_mk2(z1.shape[0], center_x, center_y)
    plt.plot(np.real(z1), np.imag(z1), "x")
    plt.plot(np.real(z2), np.imag(z2), "o")
    plt.show()
    

def circle_center_a_mk2(size, center_x, center_y):
    R = np.sqrt((1.0 - center_x) ** 2 + center_y ** 2)
    x_pre = np.linspace(start=-R+center_x, stop=R+center_x, num=size)

    y_u = np.sqrt(R**2 - (x_pre - center_x)**2) + center_y
    y_l = - np.sqrt(R**2 - (x_pre - center_x)**2) + center_y

    # x_pre = 1, y_l = 0の点を先頭に半時計まわりに並べ替えて，x, yを返したい
    # つまり，argmaxで上述の1を特定，その番号を起点に，スライスを作り、concatする
    # RE:number and concat
    y_u = y_u[::-1]
    epsilon = 1.9 * R / size

    # find index
    a_x = np.where(((1.0 - epsilon < x_pre) & (x_pre < 1.0 + epsilon)), 1, 0)
    if center_y >= 0:
        a_y = np.where(((- epsilon < y_l) & (y_l < epsilon)), 1, 0)
    else:
        a_y = np.where(((- epsilon < y_u) & (y_u < epsilon)), 1, 0)
    # index = np.argmax(a_x * a_y)
    index = np.argmax(a_x)

    # renumber index

    start = 0
    stop = size
    if (y_l[size - 1] == y_u[0]):
        start += 1

    if (y_u[size - 1] == y_l[0]):
        stop -= 1
        
    x_inv = x_pre[::-1]

    if center_y >= 0:
        y = np.concatenate([y_l[index:], y_u[start:stop], y_l[:index]])
        x = np.concatenate([x_pre[index:], x_inv[start:stop], x_pre[:index]])

    else:
        print("error, please check center_y >= 0")
        exit()
        x = np.concatenate([x_inv[index:], x_pre[start:stop], x_inv[:index]])
        y = np.concatenate([y_u[index:], y_l[start:stop], y_u[:index]])


    z = x + 1j * y
    return re_numbering(z)

def joukowski_wing(size, center_x, center_y):
    z = joukowski_wing_complex(size, center_x, center_y)
    return np.real(z), np.imag(z)

def joukowski_wing_complex(size, center_x, center_y):
    a = 1.0
    z = circle_center_a(size, center_x, center_y)
    z = joukowski_tr(z, a)
    return re_numbering(z)

def karman_trefftz_wing_complex(size, center_x, center_y):
    b = 1.0
    z = circle_center_a(size, center_x, center_y)
    z = karman_trefftz_tr(z, b)
    return re_numbering(z)

def naca_4_complex(size, int4):
    naca = Naca_4_digit(int_4=int4, attack_angle_deg=0, resolution=size)

    start = 0
    stop = size
    if (naca.equidistant_y_l[0] == naca.equidistant_y_u[size - 1]):
        start += 1

    if (naca.equidistant_y_u[0] == naca.equidistant_y_l[size - 1]):
        stop -= 1
    # これはうまくいかない(外れ値などの対応が必要なので後日書き直す)
    x = np.concatenate([naca.x_u[::-1], naca.x_l[start:stop]])
    y = np.concatenate([naca.equidistant_y_u[::-1], naca.equidistant_y_l[start:stop]])
    return x + 1j * y

def main():
    # ジューコフスキー翼のξ-η平面における円周上の点を格納する(引数は半径Rと，中心とaの角度β)
    # x-y平面における座標値z = w + a^2/wを格納
    # ReとImを分離してx, yにそれぞれ格納

    # 複素数配列を作るときは
    # 実部の配列x，虚部の配列yをそれぞれ作ったあと，z = x + 1j * yで代入する
    size = 50
    center_x = -0.08
    center_y = 0.08

    z = circle_center_a(size, center_x, center_y)
    # z = circle_center_a_mk2(size, center_x, center_y)
    z = joukowski_tr(z, 1.0)

    x = np.real(z)
    y = np.imag(z)
    plt.plot(x, y)
    plt.show()

    x_u = np.real(z[:size])
    y_u = np.imag(z[:size])
    x_l = np.real(z[size:])
    y_l = np.imag(z[size:])
    plt.plot(x_u, y_u, "o")
    plt.plot(x_l, y_l)
    plt.show()

    w = 3 + 4j
    a = 1.0
    z = w + a**2/w
    print(z)


if __name__ == '__main__':
    main()