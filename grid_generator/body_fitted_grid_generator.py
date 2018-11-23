# coding: utf-8
from math import sqrt
from scipy import interpolate
from scipy.spatial import Delaunay
import numpy as np
from numpy.linalg import norm
from naca_4digit_test import Naca_4_digit, Naca_5_digit
from joukowski_wing import joukowski_wing_complex, karman_trefftz_wing_complex
import matplotlib.pyplot as plt



# 物体表面の複素座標を取得する
def get_complex_coords(type, size, center_x=-0.08, center_y=0.08, naca4="0012"):
    def reshape_z(z):
        if z[0] != z[z.shape[0] - 1]:
            return np.concatenate([z, z[0].reshape(-1)]), z.shape[0] + 1
        else:
            return z, z.shape[0]

    # 極端に距離の近い制御点を除去する
    def adjust_length(z):
        len = np.zeros_like(z, dtype=float)
        len[:z.shape[0] - 1] = np.abs(z[1:] - z[:z.shape[0] - 1])
        len[z.shape[0] - 1] = np.abs(z[0] - z[z.shape[0] - 1])
        average_len = np.average(len)

        put_out = lambda x, count: np.hstack((x[:count], x[count + 1]))
        count = z.shape[0] - 2
        while count > 0:
            if len[count] < 0.1 * average_len:
                z = put_out(z, count)
            count -= 1
        return z

    if type == 0:
        t = np.linspace(start=0, stop=2.0 * np.pi, num=size + 1)
        z = np.exp(1j * t)[:size]
    elif type == 1:
        z = joukowski_wing_complex(size, center_x, center_y)
    elif type == 2:
        z = karman_trefftz_wing_complex(size, center_x, center_y)
    elif type == 3:
        naca = Naca_4_digit(int_4=naca4, attack_angle_deg=0.0, resolution=size, quasi_equidistant=False)
        z = naca.transform2complex()
    elif type == 4:
        naca = Naca_5_digit(int_5=naca4, attack_angle_deg=0.0, resolution=size, quasi_equidistant=False,
                            length_adjust=True)
        z = naca.transform2complex()
    else:
        print("type error")
        exit()

    z = adjust_length(z)
    if type < 3:
        return reshape_z(z)
    else:
        return z, z.shape[0]


# 物体の代表長さを得る(x,y方向どちらかが最大と仮定しており，最長部長さを求めているわけでないことに注意されたし)
def get_model_length(z, both=False):
    def get_length(x):
        return np.max(x) - np.min(x)

    x = np.real(z)
    y = np.imag(z)
    if both:
        return get_length(x), get_length(y)
    else:
        return max(get_length(x), get_length(y))

# 物体の中心位置を返す
def get_model_center(z):
    return np.average(np.real(z)), np.average(np.imag(z))

# 格子の外部境界(分割数は物体と同じの，物体形状の長軸長さのmagnification倍の円を返す)
def get_outer_boundary(z1, magnification=5, equidistant=False):

    model_length = get_model_length(z1)
    center_x, center_y = get_model_center(z1)
    zc = center_x + 1j * center_y
    radius = model_length * magnification
    if equidistant == False:
        # 法線の角度
        delta2 = get_delta2(z1)
        theta1 = np.angle(delta2 / (np.abs(delta2) * 1j))
        theta1 = np.where(theta1 > 0, theta1, theta1 + 2.0 * np.pi)

        # 物体を円に変換したときの換算角度
        theta2 = get_length_rate(z1) * 2.0 * np.pi
        average_theta = np.sort(0.5 * (theta1 + theta2))
        z3 = zc + radius * np.exp(1j * average_theta)
        return z3[::-1]  # Clock Wise and -1
    else:
        z3 = zc + radius * np.exp(1j * np.linspace(0, 2.0 * np.pi, z1.shape[0] + 1))
        return z3[1:]


# そこまでの累積長さが全体の長さに占める割合を返す
def get_length_rate(z1, output_total_length=False):
    size = z1.shape[0]
    delta1 = np.zeros_like(z1, dtype=complex)
    delta1[:size - 1] = z1[1:] - z1[:size - 1]
    delta1[size - 1] = z1[0] - z1[size - 1]
    len = np.abs(delta1)
    total_len = np.sum(len)
    len_rate = np.zeros_like(z1, dtype=float)
    accumulated_len = 0.0
    for i in range(size):
        len_rate[i] = accumulated_len / total_len
        accumulated_len += len[i]
    if output_total_length:
        return len_rate, total_len
    else:
        return len_rate


# 1点飛ばしでの座標の差分を返す
def get_delta2(z1):
    size = z1.shape[0]
    delta2 = np.zeros_like(z1, dtype=complex)
    delta2[0] = z1[1] - z1[size - 1]
    for i in range(1, size - 1):
        delta2[i] = z1[i + 1] - z1[i - 1]
    delta2[size - 1] = z1[0] - z1[size - 2]
    return delta2


# 物体と外周を結ぶ線分を返す
def get_connect_z1_to_z3(z1, z3, resolution=None, magnification=10):
    if resolution == None:
        resolution = z1.shape[0]
    else:
        resolution = resolution

    inner_end = z1[np.argmax(np.real(z1))]  # 内側のx方向最大位置
    outer_end = z3[np.argmax(np.real(z3))]  # 外側のx方向最大位置

    exp_end = np.log(magnification)  # 指定倍率:magnificationに達する指数関数上の位置:magnification = exp(exp_end) + 1

    delta_x = np.exp(np.linspace(0, exp_end, resolution - 2, dtype=complex))  # 指数関数の等間隔サンプリング
    raw_length = np.sum(delta_x)  # 等間隔サンプリングされた微小長さの総和
    delta_x = (outer_end - inner_end) / raw_length * delta_x  # 内から外への長さにスケーリング&方向を付ける

    z2 = np.zeros(resolution, dtype=complex)
    z2[0] = inner_end
    z2[resolution - 1] = outer_end
    for k in range(1, resolution - 1):
        z2[k] = z2[k - 1] + delta_x[k - 1]

    return z2


def deduplication(z, array_list=None):
    def put_out(x, count):
        return np.hstack((x[:count], x[count + 1]))

    def put_out_bound(x):
        return x[:x.shape[0] - 1]

    size = z.shape[0]

    if z[size - 1] == z[0]:
        size -= 1
        z = put_out_bound(z)
        if array_list != None:
            for i in range(len(array_list)):
                array_list[i] = put_out_bound(array_list[i])

    count = size - 2
    while count > 0:
        if z[count] == z[count + 1]:
            z = put_out(z, count)
            if array_list != None:
                for i in range(len(array_list)):
                    array_list[i] = put_out(array_list[i], count)
        count -= 1
    if array_list == None:
        return z
    else:
        return z, array_list

def Tri2vtk(path, fname, Tri_points, Tri_simplices):
    fname = path + fname + ".vtk"
    with open(fname, 'w') as f:
        point_number = Tri_points.shape[0]
        cell_number = Tri_simplices.shape[0]
        cell_vertex_number = 4 * Tri_simplices.shape[0]
        #cell_vertex_number = str(2 * 4 * (xi_max) * (eta_max - 1))


        f.write("# vtk DataFile Version 3.0\n")
        f.write("Unstructured Grid tri example\n")
        f.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
        f.write("POINTS " + str(point_number) + " double\n")
        # point coordinates
        for i in range(point_number):
            f.write(str(Tri_points[i, 0]) + " " + str(Tri_points[i, 1]) + " 0.0\n")

        # cell structure
        f.write("CELLS " + str(cell_number) + " " + str(cell_vertex_number) + "\n")
        for i in range(cell_number):
            f.write("3 " + str(Tri_simplices[i, 0]) + " " + str(Tri_simplices[i, 1]) + " " + str(Tri_simplices[i, 2]) + "\n")

        # cell types
        f.write("CELL_TYPES " + str(cell_number) + "\n")
        [f.write("5\n") for i in range(cell_number)]

# p1とp2, p3とp4が線分をなすとして
def line_intersect(p1, p2, p3, p4):
    flag = 0
    for i in range(2):
        if (p1[i] >= p2[i]):
            if((p1[i] < p3[i] and p1[i] < p4[i]) or (p2[i] > p3[i] and p2[i] > p4[i])):
                flag = 1
                break
        else:
            if((p2[i] < p3[i] and p2[i] < p4[i]) or (p1[i] > p3[i] and p1[i] > p4[i])):
                flag = 1
                break

    if flag == 1:
        return False

    c1 = (p3[0] - p4[0]) * (p1[1] - p3[1]) + (p3[1] - p4[1]) * (p3[0] - p1[0])
    c2 = (p3[0] - p4[0]) * (p2[1] - p3[1]) + (p3[1] - p4[1]) * (p3[0] - p2[0])
    c3 = (p1[0] - p2[0]) * (p3[1] - p1[1]) + (p1[1] - p2[1]) * (p1[0] - p3[0])
    c4 = (p1[0] - p2[0]) * (p4[1] - p1[1]) + (p1[1] - p2[1]) * (p1[0] - p4[0])
    return c3 * c4 < 0 and c1 * c2 < 0

def make_grid_seko(z1, path="", fname="sample", mg2=True, vtk=True, bdm=True, trianglation=True):
    xi_max = z1.shape[0]
    eta_max = z1.shape[0]  # int(0.5 * z1.shape[0])

    def make_point(z):
        return np.vstack([np.real(z), np.imag(z)]).T

    def convert_complex2real(comp):
        return np.vstack([np.real(comp), np.imag(comp)]).flatten()

    # pt_cd = [x, y] ... 1Dimension, pts_masktri = Tri.simplices[mask, :]
    def check_intersect(pt_cd, pts_masktri, tri_total):
        flag = 1
        for tri in range(tri_total):
            relative_pts = pts_masktri[tri] - pt_cd.flatten()
            sign1 = np.sign(np.cross(relative_pts[0], relative_pts[1]))
            sign2 = np.sign(np.cross(relative_pts[1], relative_pts[2]))
            sign3 = np.sign(np.cross(relative_pts[2], relative_pts[0]))
            if (sign1 == sign2) and (sign2 == sign3):
                flag = -1
                break

        if flag == 1:
            return False    # 交差なし
        else:
            return True # 交差あり

    grid_x = np.zeros((xi_max, eta_max))
    grid_y = np.zeros((xi_max, eta_max))

    def output_vtk_tri(fname, path):
        fname = path + fname + ".vtk"
        with open(fname, 'w') as f:
            point_number = str(xi_max * eta_max)
            cell_number = str(2 * (xi_max) * (eta_max - 1))
            cell_vertex_number = str(2 * 4 * (xi_max) * (eta_max - 1))
            pid = lambda i, j: i + xi_max * j

            # 四角形セルの対角線長さを計算し，左下→右上の対角線(diag1)が，左上→右下の対角線(diag2)より短いときTrueを返す
            def set_relatively_not_flat_tri(i, ip1, j):
                diag1 = (grid_x[ip1, j + 1] - grid_x[i, j]) ** 2 + (grid_y[ip1, j + 1] - grid_y[i, j]) ** 2
                diag2 = (grid_x[i, j + 1] - grid_x[ip1, j]) ** 2 + (grid_y[i, j + 1] - grid_y[ip1, j]) ** 2

                if diag1 < diag2:
                    return True
                else:
                    return False

            def cell_structure(i, j):
                # 短い対角線を新しい辺として四角形を三角形に分割する
                if i != xi_max - 1:
                    ip1 = i + 1
                else:
                    ip1 = 0
                flag = set_relatively_not_flat_tri(i, ip1, j)
                # flag == True:左下と右上を結ぶ対角線を新たな辺として，左上の三角形と右下の三角形に分割
                if flag:
                    tri1 = "3 " + str(pid(i, j)) + " " + str(pid(ip1, j + 1)) + " " + str(pid(i, j + 1))
                    tri2 = "3 " + str(pid(i, j)) + " " + str(pid(ip1, j)) + " " + str(pid(ip1, j + 1))
                else:
                    tri1 = "3 " + str(pid(i, j)) + " " + str(pid(ip1, j)) + " " + str(pid(i, j + 1))
                    tri2 = "3 " + str(pid(ip1, j)) + " " + str(pid(ip1, j + 1)) + " " + str(pid(i, j + 1))
                return tri1 + "\n" + tri2

            # header
            f.write("# vtk DataFile Version 3.0\n")
            f.write("Unstructured Grid tri example\n")
            f.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
            f.write("POINTS " + point_number + " double\n")
            # point coordinates
            for j in range(eta_max):
                for i in range(xi_max):
                    f.write(str(grid_x[i, j]) + " " + str(grid_y[i, j]) + " 0.0\n")

            # cell structure
            f.write("CELLS " + cell_number + " " + cell_vertex_number + "\n")
            for j in range(eta_max - 1):
                for i in range(xi_max):
                    f.write(cell_structure(i, j) + "\n")

            # cell types
            f.write("CELL_TYPES " + cell_number + "\n")
            for j in range(eta_max - 1):
                for i in range(xi_max):
                    f.write("5\n5\n")

    def plot_tmp(detail=True, zoom=False, overview=True):
        if detail:
            for i in range(xi_max):
                plt.plot(grid_x[i, :], grid_y[i, :])
            for j in range(eta_max):
                plt.plot(grid_x[:, j], grid_y[:, j])
            plt.xlim(-0.2, 1.2)
            plt.ylim(0.3, 0.7)
            plt.show()
        if zoom:
            for i in range(xi_max):
                plt.plot(grid_x[i, :], grid_y[i, :])
            for j in range(eta_max):
                plt.plot(grid_x[:, j], grid_y[:, j])
            plt.xlim(-2.3, -1.7)
            plt.ylim(-0.3, 0.3)
            plt.show()
        if overview:
            for i in range(xi_max):
                plt.plot(grid_x[i, :], grid_y[i, :])
            for j in range(eta_max):
                plt.plot(grid_x[:, j], grid_y[:, j])
            plt.show()

    # 物体表面を押し出す(η格子線の複素座標，オフセット方向(外側or内側)，最大のオフセット量，オフセット倍率(遠方領域で1以上の値を与えて計算領域を効率的に広げる))
    def offset_surface(z, outer=False, max_incremental=0.1, accel=1.0):
        size = z.shape[0]

        delta = np.zeros(size, dtype=complex)
        delta[0] = z[1] - z[size - 1]
        for i in range(1, size - 1):
            delta[i] = z[i + 1] - z[i - 1]
        delta[size - 1] = z[0] - z[size - 2]

        normal = -1j * delta / np.abs(delta)

        # 格子の裏返り防止(隣接する格子点から外向きξ方向へ伸びる2つの格子線がなす角度が90°以上開いている場合に，新しいη格子線が既存のη格子線と被らないように角度を修正する)
        def prevent_turn_over_cell(i, imp1, downwind=True):
            dot_normal = np.real(normal[i] * np.conj(normal[imp1]))  # 隣接する格子線同士の内積を取る
            flag = 0
            if dot_normal < 0:  # 内積が負のとき修正対象
                flag = 1
                angle1 = np.angle(normal[i])
                angle0 = np.angle(normal[imp1])
                if angle1 < angle0:
                    angle1 += 2.0 * np.pi
                if downwind:  # ξが小さい側を修正する場合
                    normal[i] -= (angle1 - angle0) / 100 * 5  # 角度差の5%だけ寄せる
                else:  # ξが大きい側を変更する場合(風下と風上とを2回セットとして修正する)
                    normal[i] += (angle1 - angle0) / 100 * 5  # 角度差の5%だけ寄せる
            return flag

        get_positive_angle = lambda angle: np.where(np.angle(angle) < 0.0, np.angle(angle) + 2.0 * np.pi,
                                                    np.angle(angle))
        flag = 0
        count = 0
        tmp_normal = normal.copy()
        while flag != 0:  # flag=0となるまで格子の裏返り防止処理を掛ける
            flag = 0
            count += 1
            flag += prevent_turn_over_cell(0, xi_max - 1)
            for i in range(1, size):
                flag += prevent_turn_over_cell(i, i - 1)
            for i in range(size - 1):
                flag += prevent_turn_over_cell(i, i + 1, downwind=False)
            flag += prevent_turn_over_cell(xi_max - 1, 0, downwind=False)
            if count == 1000:
                flag = True

        if flag:
            normal = tmp_normal

        if np.min(np.abs(delta)) > 0.1 * np.average(np.abs(delta)):
            incremental = accel * min(min(2.0 / np.pi * np.min(np.abs(delta)), np.average(np.abs(delta))),
                                      max_incremental)  # obj sizeとboundaryサイズを均等に分割したときの幅で置換すべき(0.1)
        else:
            incremental = accel * 0.5 * np.average(np.abs(delta))  # obj sizeとboundaryサイズを均等に分割したときの幅で置換すべき(0.1)

        if outer == True:
            return z - normal * incremental
        else:
            return z + normal * incremental

    # jacobi mk2(point wise iterate)
    def get_im2_im1_ip1_and_ip2(i):
        if i == 1:
            return xi_max - 1, 0, 2, 3
        if i == 0:
            return xi_max - 2, xi_max - 1, 1, 2
        elif i == xi_max - 1:
            return xi_max - 3, xi_max - 2, 0, 1
        elif i == xi_max - 2:
            return xi_max - 4, xi_max - 3, xi_max - 1, 0
        else:
            return i - 2, i - 1, i + 1, i + 2

    sum_i_x = lambda im1, ip1, j: grid_x[im1, j] + grid_x[ip1, j]
    sum_i_y = lambda im1, ip1, j: grid_y[im1, j] + grid_y[ip1, j]

    sum_j_x = lambda i, j: grid_x[i, j - 1] + grid_x[i, j + 1]
    sum_j_y = lambda i, j: grid_y[i, j - 1] + grid_y[i, j + 1]

    sum_i_x2 = lambda im2, im1, ip1, ip2, j: 0.5 * (sum_i_x(im1, ip1, j) + sum_i_x(im2, ip2, j))
    sum_i_y2 = lambda im2, im1, ip1, ip2, j: 0.5 * (sum_i_y(im1, ip1, j) + sum_i_y(im2, ip2, j))

    def sum_j_x2(i, j):
        if (j == 1) or (j == eta_max - 2):
            return sum_j_x(i, j)
        else:
            return 0.5 * (sum_j_x(i, j) + grid_x[i, j - 2] + grid_x[i, j + 2])

    def sum_j_y2(i, j):
        if (j == 1) or (j == eta_max - 2):
            return sum_j_y(i, j)
        else:
            return 0.5 * (sum_j_y(i, j) + grid_y[i, j - 2] + grid_y[i, j + 2])

    sum_i_j_x = lambda im1, ip1, j: grid_x[ip1, j + 1] - grid_x[ip1, j - 1] - grid_x[im1, j + 1] + grid_x[im1, j - 1]
    sum_i_j_y = lambda im1, ip1, j: grid_y[ip1, j + 1] - grid_y[ip1, j - 1] - grid_y[im1, j + 1] + grid_y[im1, j - 1]

    def get_equidistant_curve(z2):
        def func4(x):
            return np.where(x < 1 / 4, 4 * x ** 2,
                            np.where(x < 2 / 4, -4 * x ** 2 + 4 * x - 1 / 2,
                                     np.where(x < 3 / 4, 4 * x ** 2 - 4 * x + 3 / 2, -4 * x ** 2 + 8 * x - 3)))

        t, total_len = get_length_rate(z2, output_total_length=True)
        fx = interpolate.PchipInterpolator(np.hstack((t, np.array([1.0]))), np.real(np.hstack((z2, z2[0]))))
        fy = interpolate.PchipInterpolator(np.hstack((t, np.array([1.0]))), np.imag(np.hstack((z2, z2[0]))))
        equidistant_t = np.linspace(0, 1, z2.shape[0] + 1)[:z2.shape[0]]
        return fx(equidistant_t) + 1j * fy(equidistant_t)

    def equidistant_offset(z2, max_incremental, accel):
        z2 = offset_surface(z2, outer=True, max_incremental=max_incremental, accel=accel)
        return get_equidistant_curve(z2)

    get_object_center = lambda z2: np.average(np.real(z2)) + 1j * np.average(np.imag(z2))

    magnification = 5.0
    max_incremental = 10.0  # radius1 * (magnification - 1) / eta_max
    accel_parameter = 1.3  # 物体遠方領域でオフセット量を増やす際の割合

    # accel量の設定
    def set_accel(j, accel_parameter):
        if j < int(eta_max / 2):
            return 1.0
        else:
            return accel_parameter

    z1_eq = get_equidistant_curve(z1)
    grid_x[:, 0] = np.real(z1_eq)
    grid_y[:, 0] = np.imag(z1_eq)

    z2 = equidistant_offset(z1_eq, max_incremental, accel=1.0)
    z2_eq = z2

    # z2 = np.hstack((z2[1:], z2[0]))
    circular_grid = False
    if circular_grid:
        for j in range(1, eta_max - 1):
            grid_x[:, j] = np.real(z2)
            grid_y[:, j] = np.imag(z2)
            z2 = equidistant_offset(z2, max_incremental, accel=1.0)

        obj_center = get_object_center(z2)
        phi = np.angle(z2[0] - obj_center)

        radius = 0.5 * get_model_length(z1) * magnification
        z3 = obj_center + radius * np.exp(1j * np.linspace(phi, phi + 2.0 * np.pi, z1.shape[0] + 1))[1:][::-1]
    else:
        model_length = get_model_length(z1)
        for j in range(1, eta_max):
            grid_x[:, j] = np.real(z2)
            grid_y[:, j] = np.imag(z2)

            flag = 0
            mix_rate = 1.0 - max(np.exp(-0.18 * float(j - 1)),0.7)  # j番目のη線と直交するように与えたz2_orthogonalと，Δηが均等になるように与えたz2_equidistantの配合比率
            while flag >= 0:
                flag += 1
                z2_equidistant = equidistant_offset(z2, max_incremental, set_accel(j, accel_parameter))
                z2_orthogonal = offset_surface(z2, outer=True, max_incremental=max_incremental,
                                               accel=set_accel(j, accel_parameter))
                fix_z2 = (1.0 - mix_rate) * z2_orthogonal + mix_rate * z2_equidistant
                delta_j__ = np.hstack((z2[1:] - z2[:xi_max - 1], z2[0] - z2[xi_max - 1]))
                delta_jp1 = np.hstack((fix_z2[1:] - fix_z2[:xi_max - 1], fix_z2[0] - fix_z2[xi_max - 1]))


                if np.any(np.real(delta_j__ * np.conj(delta_jp1)) < 0):
                    mix_rate = 1.0 - max(-0.18 * np.exp(-float(j - 2)), 0.7 - 0.1 * flag)
                    print(mix_rate)
                    if mix_rate > 0.9:
                        z2 = fix_z2
                        break
                else:
                    flag = -1
                    z2 = fix_z2

            if (np.max(np.real(z2)) - np.min(np.real(z2)) > 40.0 * model_length):
                eta_max = j + 2
                break

        z3 = z2

    grid_x[:, eta_max - 1] = np.real(z3)
    grid_y[:, eta_max - 1] = np.imag(z3)
    # grid_x, grid_y = jacobi_pre(grid_x, grid_y)
    grid_x, grid_y = grid_x[:, :eta_max], grid_y[:, :eta_max]

    # ここまででおおよその点の追加が終了
    # このまま三角形を形成すると段階では物体を貫通していたり，都合のよくない線が存在している可能性がある
    # そこで、物体を三角形領域に分割し，その三角形群と格子線の接触判定を行う
    def set_long_axis_direction(z1):
        if (np.max(np.real(z1)) - np.min(np.real(z1))) > np.max(np.imag(z1)) - np.min(np.imag(z1)):
            return True
        else:
            return False

    def split_surface(z):
        direction = set_long_axis_direction(z)
        if direction:
            start = np.argmin(np.real(z))
            end = np.argmax(np.real(z))
        else:
            start = np.argmin(np.imag(z))
            end = np.argmax(np.imag(z))

        if start < end:
            z_upper = z[start:end + 1]
            z_lower = np.concatenate([z[end:], z[:start + 1]])
        else:
            z_upper = z[end:start + 1]
            z_lower = np.concatenate([z[start:], z[:end + 1]])
        return z_upper, z_lower

    def obj_trianglization(z1):
        z_u, z_l = split_surface(z1)
        new_z_eq = np.concatenate([z_u, z_l[1:z_l.shape[0]-1]])
        size = new_z_eq.shape[0]

        def plot_complex(z):
            plt.plot(np.real(z), np.imag(z), "x")

        def set_r(i):
            return i

        def set_l(i, size):
            return  size - i    # -1ではない

        cell_total = new_z_eq.shape[0] - 2
        tri_point = np.vstack([np.real(new_z_eq), np.imag(new_z_eq)]).T
        tri_block = np.zeros((cell_total, 3), dtype=int)
        process = 0
        for i in range(int(cell_total/2)):
            tri_block[2 * i, :] = np.array([set_r(i), set_r(i+1), set_l(i + 1, size)])
            tri_block[2 * i + 1, :] = np.array([set_r(i + 1), set_l(i + 1, size), set_l(i + 2, size)])
            process += 2

        if process != cell_total:
            i = int(cell_total/2)
            tri_block[cell_total - 1, :] = np.array([set_r(i + 1), set_l(i + 1, size), set_l(i + 2, size)])

        return tri_point, tri_block

    # 物体内部の三角形化
    obj_tri_pts, obj_tri_spx = obj_trianglization(z1)

    # 物体外部+内部の三角形化
    pointj1 = np.vstack([grid_x[:, 1].flatten(), grid_y[:, 1].flatten()]).T # η=1線は独立させておく
    points2D = np.vstack([grid_x[:, 2:].flatten(), grid_y[:, 2:].flatten()]).T
    points2D = np.concatenate([pointj1, points2D])

    Tri2 = Delaunay(points2D)

    # 品質の悪いセルを除去する
    calc_Area2 = lambda pts3: np.cross(pts3[1, :] - pts3[0, :], pts3[2, :] - pts3[0, :])  # 実際は1/2にする必要がある
    minimum_edge = lambda pts3: np.argmin(np.array([norm(pts3[1, :] - pts3[0, :]), norm(pts3[2, :] - pts3[1, :]), norm(pts3[0, :] - pts3[2, :])]))
    minimum_area2_half = 0.5 * np.average(np.abs(z1_eq[1:] - z1_eq[:z1_eq.shape[0] - 1])) ** 2

    total_add_cell = Tri2.simplices.shape[0]
    total_add_point = points2D.shape[0]
    remove_point = np.zeros(total_add_point)
    for iCell in range(total_add_cell):
        if all(Tri2.simplices[iCell] > xi_max) and (calc_Area2(Tri2.points[Tri2.simplices[iCell]]) < minimum_area2_half):    # η=1線の格子点を含まず，物体近傍のセル面積の半分以下の面積しか持たない要素について
            target = minimum_edge(Tri2.points[Tri2.simplices[iCell]])
            if Tri2.simplices[iCell][target] > Tri2.simplices[iCell][(target + 1) % 3]:
                remove_point[Tri2.simplices[iCell][target]] = 1
            else:
                remove_point[Tri2.simplices[iCell][(target + 1) % 3]] = 1

    remove_point = np.where(remove_point == 1, False, True) # Falseを削除，Trueを残す
    points2D = points2D[remove_point]
    Tri2 = Delaunay(points2D)
    total_add_point = points2D.shape[0]

    # 追加する点の総数が分かった時点で，物体表面→η=1線の格子を先に切る
    # 物体表面のη=0線と物体表面から少し外側のη=1格子線を三角形で繋ぐ
    num = z1_eq.shape[0]
    edge_mask = np.ones((num, num), dtype=int)
    length = np.zeros((num, num), dtype=float)
    for i in range(num):
        for j in range(num):
            length[i, j] = np.abs(z1_eq[i] - z2_eq[j])

    z1_xy = make_point(z1_eq)
    z2_xy = make_point(z2_eq)

    # 1.そもそも物体表面と交差してるのを除外
    for i in range(num):
        for j in range(num):
            for k in range(num):
                kp1 = k + 1
                if kp1 == num:
                    kp1 = 0
                if (i != k and i != kp1):
                    if line_intersect(z1_xy[i], z2_xy[j], z1_xy[k], z1_xy[kp1]):
                        edge_mask[i, j] = 0


    # 2.線分同士で交差してるのを除外(bluteforth)
    for i in range(num):    # p1
        for j in range(num):    # p2
            if edge_mask[i, j] != 0:    # まだ除外されてない線分p1-p2のうち
                for k in range(num):    # p3
                    if k != i:  # p1 = p3は交差してないことにする
                        for l in range(num):    # p4
                            if l != j:  # p2 = p4も交差してないことにする
                                if(edge_mask[k, l] != 0):   # p3-p4がまだ除外されていなければ
                                    if line_intersect(z1_xy[i], z2_xy[j], z1_xy[k], z2_xy[l]):  # 交差してたら長い方を消す
                                        if length[i, j] > length[k, l]:
                                            edge_mask[i, j] = 0
                                        else:
                                            edge_mask[k, l] = 0

    # 3.残った辺から三角形を構築
    simplices = []  # 最終的な格子に登録する用
    for i in range(num):    # p1
        for j in range(num):    # p3
            if edge_mask[i, j] == 1:    # 辺p1-p3が生き残っていたときに
                jp1 = j + 1
                if jp1 == num:
                    jp1 = 0

                if edge_mask[i, jp1] == 1:
                    simplices.append([i + total_add_point, j, jp1])

    for i in range(num):    # p2
        for j in range(num):    # p1
            if edge_mask[j, i] == 1:    # 辺p1-p2が生き残っており
                jp1 = j + 1
                if jp1 == num:
                    jp1 = 0
                if edge_mask[jp1, i] == 1:    # 次に生き残っている辺p2-p3について
                    simplices.append([i, j + total_add_point, jp1 + total_add_point])

    simplices = (np.array(simplices))
    new_edge = []
    new_edge_num = 0
    for i in range(num):
        for j in range(num):
            if edge_mask[i, j] == 1:
                new_edge.append([i, j])
                new_edge_num += 1
    new_edge = np.array(new_edge)
    # ここまで物体側

    # 物体外部+内部の格子について物体近傍の三角形のみを取り出す(物体の大きさの1.1倍程度のboudanry box)
    len_x, len_y = get_model_length(z1, both=True)
    center = get_object_center(z1)
    param = 0.55
    min_x = np.real(center) - param * len_x
    max_x = np.real(center) + param * len_x
    min_y = np.imag(center) - param * len_y
    max_y = np.imag(center) + param * len_y

    total_add_cell = Tri2.simplices.shape[0]
    checK_target = np.zeros(total_add_cell)

    for iCell in range(total_add_cell):
        pts = Tri2.points[Tri2.simplices[iCell, :]]
        for iPoint in range(3):
            if (min_x < pts[iPoint, 0] < max_x) and (min_y < pts[iPoint, 1] < max_y):
                checK_target[iCell] = 1

    mask = np.where(checK_target == 1, True, False)


    def tri_intersect(pts3, obj_tri_pts, obj_tri_spx):
        pattern = np.array([[0,1], [1,2], [2,0]])
        cross = 0
        for tri in obj_tri_spx:
            obj_pts3 = obj_tri_pts[tri]
            for out_pat in pattern:
                for obj_pat in pattern:
                    chk = line_intersect(pts3[out_pat[0]], pts3[out_pat[1]], obj_pts3[obj_pat[0]], obj_pts3[obj_pat[1]])
                    if chk:
                        cross = 1
                        break
                if cross == 1:
                    break
            if cross == 1:
                break

        if cross == 0:
            for iEdge in range(new_edge_num):
                i = new_edge[iEdge, 0]
                j = new_edge[iEdge, 1]
                for out_pat in pattern:
                    chk = line_intersect(pts3[out_pat[0]], pts3[out_pat[1]], z1_xy[i], z2_xy[j])
                    if chk:
                        cross = 1
                        break
                if cross == 1:
                    break

        return cross == 1


    # 物体を通過する線分を保有する三角形要素 (&新規追加する辺と接触する三角形要素)の除去(これを行うためのマスク作成)
    for iCell in range(total_add_cell):
        if mask[iCell]:
            pts3 = Tri2.points[Tri2.simplices[iCell], :]
            if tri_intersect(pts3, obj_tri_pts, obj_tri_spx) == False:
                mask[iCell] = False


    # 4.最後方にη=0の格子点を追加
    pointj0 = np.vstack([grid_x[:, 0].flatten(), grid_y[:, 0].flatten()]).T
    grid_pts = np.concatenate([points2D, pointj0])
    grid_simplices = np.concatenate([Tri2.simplices[~mask, :], simplices])

    Tri2vtk(path="", fname="test", Tri_points=grid_pts, Tri_simplices=grid_simplices)
    exit()

    grid_x = np.vstack((grid_x[:, :], grid_x[0, :]))
    grid_y = np.vstack((grid_y[:, :], grid_y[0, :]))

    return grid_x, grid_y




def make_grid(fname, type, size=100, naca4="0012", center_x=0.08, center_y=0.08, mayugrid2=False, vtk=False, bdm=False,
              trianglation=True, path=""):
    z1, size = get_complex_coords(type=type, center_x=center_x, center_y=center_y, naca4=naca4, size=size)
    z1 = deduplication(z1)[::-1]
    make_grid_seko(z1, path, fname, mayugrid2, vtk, bdm, trianglation)

def main():
    z1, size = get_complex_coords(type=3, naca4="4101", size=50)
    # z1, size = get_complex_coords(type=0, center_x=0.08, center_y=0.3, naca4="4912", size=100)
    z1 = deduplication(z1)[::-1]
    make_grid_seko(z1)
    # plt.plot(np.real(z1), np.imag(z1))
    # plt.show()


def makeGridLoop():
    header = "NACA"
    for i in range(1, 9999, 2):
        naca4 = str(i).zfill(4)
        fname = header + naca4
        make_grid(fname, type=3, naca4=naca4)


if __name__ == '__main__':
    main()
    # makeGridLoop()