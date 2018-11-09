# coding: utf-8
from math import sqrt
from scipy import interpolate
import numpy as np
from naca_4digit_test import Naca_4_digit, Naca_5_digit
from joukowski_wing import joukowski_wing_complex, karman_trefftz_wing_complex
import matplotlib.pyplot as plt
import Lib.bisect as bisect

# 物体表面の複素座標を取得する
def get_complex_coords(type, size, center_x = -0.08, center_y = 0.08, naca4 = "0012"):
    def reshape_z(z):
        if z[0] != z[z.shape[0] - 1]:
            return np.concatenate([z, z[0].reshape(-1)]), z.shape[0] + 1
        else:
            return z, z.shape[0]
    
    # 極端に距離の近い制御点を除去する
    def adjust_length(z):
        len = np.zeros_like(z, dtype = float)
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
        t = np.linspace(start = 0, stop = 2.0 * np.pi, num = size + 1)
        z = np.exp(1j * t)[:size]
    elif type == 1:
        z = joukowski_wing_complex(size, center_x, center_y)
    elif type == 2:
        z = karman_trefftz_wing_complex(size, center_x, center_y)
    elif type == 3:
        naca = Naca_4_digit(int_4 = naca4, attack_angle_deg = 0.0, resolution = size, quasi_equidistant = False)
        z = naca.transform2complex()
    elif type == 4:
        naca = Naca_5_digit(int_5 = naca4, attack_angle_deg = 0.0, resolution = size, quasi_equidistant = False,
                            length_adjust = True)
        z = naca.transform2complex()
    else:
        print("type error")
        exit()
    
    z = adjust_length(z)
    if type < 3:
        return reshape_z(z)
    else:
        return z, z.shape[0]


# 格子の外部境界(分割数は物体と同じの，物体形状の長軸長さのmagnification倍の円を返す)
def get_outer_boundary(z1, magnification=5, equidistant=False):
    # 物体の代表長さを得る(x,y方向どちらかが最大と仮定しており，最長部長さを求めているわけでないことに注意されたし)
    def get_model_length(z):
        def get_length(x):
            return np.max(x) - np.min(x)

        x = np.real(z)
        y = np.real(z)
        return max(get_length(x), get_length(y))

    # 物体の中心位置を返す
    def get_model_center(z):
        return np.average(np.real(z)), np.average(np.imag(z))

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
        return z3[::-1] # Clock Wise and -1
    else:
        z3 = zc + radius * np.exp(1j * np.linspace(0, 2.0 * np.pi, z1.shape[0] + 1))
        return z3[1:]


# そこまでの累積長さが全体の長さに占める割合を返す
def get_length_rate(z1, output_total_length = False):
    size = z1.shape[0]
    delta1 = np.zeros_like(z1, dtype = complex)
    delta1[:size - 1] = z1[1:] - z1[:size - 1]
    delta1[size - 1] = z1[0] - z1[size - 1]
    len = np.abs(delta1)
    total_len = np.sum(len)
    len_rate = np.zeros_like(z1, dtype = float)
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
    delta2 = np.zeros_like(z1, dtype = complex)
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

    exp_end = np.log(magnification) # 指定倍率:magnificationに達する指数関数上の位置:magnification = exp(exp_end) + 1

    delta_x = np.exp(np.linspace(0, exp_end, resolution - 2, dtype=complex))    # 指数関数の等間隔サンプリング
    raw_length = np.sum(delta_x)    # 等間隔サンプリングされた微小長さの総和
    delta_x = (outer_end - inner_end) / raw_length * delta_x # 内から外への長さにスケーリング&方向を付ける

    z2 = np.zeros(resolution, dtype=complex)
    z2[0] = inner_end
    z2[resolution - 1] = outer_end
    for k in range(1, resolution-1):
        z2[k] = z2[k-1] + delta_x[k - 1]

    return z2
    
def deduplication(z, array_list=None):
    def put_out(x, count):
        return np.hstack((x[:count], x[count+1]))

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

def make_grid_seko(z1, path = "", fname="sample", mg2=True, vtk=True, bdm=True, trianglation=True):
    xi_max = z1.shape[0]
    eta_max = z1.shape[0]   # int(0.5 * z1.shape[0])

    grid_x = np.zeros((xi_max, eta_max))
    grid_y = np.zeros((xi_max, eta_max))

    def x_xi(i, j):
        if i == 0:
            return 0.5 * (grid_x[1, j] - grid_x[xi_max - 1, j]) # loop boundary
        elif i == xi_max - 1:
            return 0.5 * (grid_x[0, j] - grid_x[xi_max - 2, j])
        else:
            return 0.5 * (grid_x[i+1, j] - grid_x[i-1, j])

    def y_xi(i, j):
        if i == 0:
            return 0.5 * (grid_y[1, j] - grid_y[xi_max - 1, j])  # loop boundary
        elif i == xi_max - 1:
            return 0.5 * (grid_y[0, j] - grid_y[xi_max - 2, j])
        else:
            return 0.5 * (grid_y[i+1, j] - grid_y[i-1, j])

    def x_eta(i, j):
        if j == 0:
            return 0.5 * (-grid_x[i, 2] + 4.0 * grid_x[i, 1] - 3.0 * grid_x[i, 0])
        elif j == eta_max - 1:
            return - 0.5 * (-grid_x[i, eta_max - 3] + 4.0 * grid_x[i, eta_max - 2] - 3.0 * grid_x[i, eta_max - 1])
        else:
            return 0.5 * (grid_x[i, j+1] - grid_x[i, j-1])

    def y_eta(i, j):
        if j == 0:
            return 0.5 * (-grid_y[i, 2] + 4.0 * grid_y[i, 1] - 3.0 * grid_y[i, 0])
        elif j == eta_max - 1:
            return - 0.5 * (-grid_y[i, eta_max - 3] + 4.0 * grid_y[i, eta_max - 2] - 3.0 * grid_y[i, eta_max - 1])
        else:
            return 0.5 * (grid_y[i, j+1] - grid_y[i, j-1])

    def x_xixi(i, j):
        if i == 0:
            return grid_x[1, j] - 2.0 * grid_x[0, j] + grid_x[xi_max-1, j]
        elif i == xi_max - 1:
            return grid_x[0, j] - 2.0 * grid_x[xi_max-1, j] + grid_x[xi_max-2, j]
        else:
            return grid_x[i+1, j] - 2.0 * grid_x[i, j] + grid_x[i-1, j]

    def y_xixi(i, j):
        if i == 0:
            return grid_y[1, j] - 2.0 * grid_y[0, j] + grid_y[xi_max-1, j]
        elif i == xi_max - 1:
            return grid_y[0, j] - 2.0 * grid_y[xi_max-1, j] + grid_y[xi_max-2, j]
        else:
            return grid_y[i+1, j] - 2.0 * grid_y[i, j] + grid_y[i-1, j]

    def x_etaeta(i, j):
        if j == 0:
            return 2.0 * grid_x[i, 0] - 5.0 * grid_x[i, 1] + 4.0 * grid_x[i, 2] - grid_x[i, 3]
        elif j == eta_max - 1:
            return -(2.0 * grid_x[i, eta_max-1] - 5.0 * grid_x[i, eta_max-2] + 4.0 * grid_x[i, eta_max-3] - grid_x[i, eta_max-4])
        else:
            return grid_x[i, j+1] - 2.0 * grid_x[i, j] + grid_x[i, j-1]

    def y_etaeta(i, j):
        if j == 0:
            return 2.0 * grid_y[i, 0] - 5.0 * grid_y[i, 1] + 4.0 * grid_y[i, 2] - grid_y[i, 3]
        elif j == eta_max - 1:
            return -(2.0 * grid_y[i, eta_max-1] - 5.0 * grid_y[i, eta_max-2] + 4.0 * grid_y[i, eta_max-3] - grid_y[i, eta_max-4])
        else:
            return grid_y[i, j+1] - 2.0 * grid_y[i, j] + grid_y[i, j-1]

    def x_xieta(i, j):
        if i == 0:
            i_l = xi_max - 1
            i_r = 1
        elif i == xi_max - 1:
            i_l = xi_max - 2
            i_r = 0
        else:
            i_l = i - 1
            i_r = i + 1

        if j == 0:
            return 0.25 * ((-grid_x[i_r, 2] + 4.0*grid_x[i_r, 1] - 3.0 * grid_x[i_r, 0])
                           - (-grid_x[i_l, 2] + 4.0*grid_x[i_l, 1] - 3.0 * grid_x[i_l, 0]))
        elif j == eta_max - 1:
            return -0.25 * ((-grid_x[i_r, eta_max-3] + 4.0*grid_x[i_r, eta_max-2] - 3.0 * grid_x[i_r, eta_max-1])
                           - (-grid_x[i_l, eta_max-3] + 4.0*grid_x[i_l, eta_max-2] - 3.0 * grid_x[i_l, eta_max-1]))
        else:
            return 0.5 * (
                grid_x[i_r, j + 1] - grid_x[i_r, j] - grid_x[i, j + 1] + 2.0 * grid_x[i, j] - grid_x[i, j - 1] -
                grid_x[i_l, j] + grid_x[i_l, j - 1])

    def y_xieta(i, j):
        if i == 0:
            i_l = xi_max - 1
            i_r = 1
        elif i == xi_max - 1:
            i_l = xi_max - 2
            i_r = 0
        else:
            i_l = i - 1
            i_r = i + 1

        if j == 0:
            return 0.25 * ((-grid_y[i_r, 2] + 4.0 * grid_y[i_r, 1] - 3.0 * grid_y[i_r, 0])
                           - (-grid_y[i_l, 2] + 4.0 * grid_y[i_l, 1] - 3.0 * grid_y[i_l, 0]))
        elif j == eta_max - 1:
            return -0.25 * (
                        (-grid_y[i_r, eta_max - 3] + 4.0 * grid_y[i_r, eta_max - 2] - 3.0 * grid_y[i_r, eta_max - 1])
                        - (-grid_y[i_l, eta_max - 3] + 4.0 * grid_y[i_l, eta_max - 2] - 3.0 * grid_y[i_l, eta_max - 1]))
        else:
            return 0.5 * (
                    grid_y[i_r, j + 1] - grid_y[i_r, j] - grid_y[i, j + 1] + 2.0 * grid_y[i, j] - grid_y[i, j - 1] -
                    grid_y[i_l, j] + grid_y[i_l, j - 1])

    g11 = lambda i, j: x_xi(i, j)**2 + y_xi(i, j)**2
    g12 = lambda i, j: x_xi(i, j) * x_eta(i, j) + y_xi(i, j) * y_eta(i, j)
    g22 = lambda i, j: x_eta(i, j)**2 + y_eta(i, j)**2

    # xi線をxi0線へ近づける制御関数P
    def control_function_of_xi(xi, eta, xi_line, eta_line, a=10.0, b=0.0, c=1.0, d=0.0,  not_move_xi=False):
        if (eta != 1) and (eta != eta_max - 2):
            if not_move_xi:
                return 0
            else:
                p = 0
                for xi0 in xi_line:
                    p += - a * np.sign(xi - xi0) * np.exp(-c * np.abs(xi - xi0))
                return p
        else:
            i = xi
            j = eta
            return -(-(x_xi(i, j) * x_xixi(i, j) + y_xi(i, j) * y_xixi(i, j)) / g11(i, j)
                    -(x_xi(i, j) * x_etaeta(i, j) + y_xi(i, j) * y_etaeta(i, j)) / g22(i, j))

    # eta線をeta0線へ近づける制御関数Q
    def control_function_of_eta(xi, eta, xi_line, eta_line, a=10.0, b=0.0, c=1.0, d=0.0, not_move_eta=False):
        if not_move_eta:
            return 0
        if (eta != 1) and (eta != eta_max - 2):
            q = 0
            for eta0 in eta_line:
                q += - a * np.sign(eta - eta0) * np.exp(-c * np.abs(eta - eta0))
            return q
        else:
            i = xi
            j = eta
            return -(-(x_eta(i, j) * x_etaeta(i, j) + y_eta(i, j) * y_etaeta(i, j)) / g22(i, j)
                    -(x_eta(i, j) * x_xixi(i, j) + y_eta(i, j) * y_xixi(i, j)) / g11(i, j))

    def update_control_function(xi_line = [0], eta_line = [0], not_move_xi = False, not_move_eta = False):
        for i in range(xi_max):
            for j in range(eta_max):
                control_P[i, j] = control_function_of_xi(i, j, xi_line, eta_line, not_move_xi=not_move_xi)
                control_Q[i, j] = control_function_of_eta(i, j, xi_line, eta_line, not_move_eta=not_move_eta)
        return control_P, control_Q

    control_P = np.zeros((xi_max, eta_max))
    control_Q = np.zeros((xi_max, eta_max))

    def rhs_x(i, j):
        if (j != 1) and (j != eta_max-2):
            return ((g22(i, j) * (x_xixi(i, j) + control_P[i, j] * x_xi(i, j)))
                    + (g11(i, j) * (x_etaeta(i, j) + control_Q[i, j] * x_eta(i, j)))
                    - 2.0 * g12(i, j) * x_xieta(i, j))
        else:
            return (g22(i, j) * x_xixi(i, j) + control_P[i, j] * x_xi(i, j)
                    + g11(i, j) * x_etaeta(i, j) + control_Q[i, j] * x_eta(i, j))

    def rhs_y(i, j):
        if (j != 1) and (j != eta_max - 2):
            return ((g22(i, j) * (y_xixi(i, j) + control_P[i, j] * y_xi(i, j)))
                    + (g11(i, j) * (y_etaeta(i, j) + control_Q[i, j] * y_eta(i, j)))
                    - 2.0 * g12(i, j) * y_xieta(i, j))
        else:
            return (g22(i, j) * y_xixi(i, j) + control_P[i, j] * y_xi(i, j)
                    + g11(i, j) * y_etaeta(i, j) + control_Q[i, j] * y_eta(i, j))

    def output_vtk(fname, path):
        fname = path + fname + ".vtk"
        with open(fname, 'w') as f:
            point_number = str(xi_max * eta_max)
            cell_number = str((xi_max) * (eta_max - 1))
            cell_vertex_number = str(5 * (xi_max) * (eta_max - 1))
            pid = lambda i, j: i + xi_max * j

            def cell_structure(i, j):
                if i != xi_max - 1:
                    ip1 = i + 1
                else:
                    ip1 = 0
                return "4 " + str(pid(i, j)) + " " + str(pid(ip1, j)) + " " + str(pid(ip1, j + 1)) + " " + str(pid(i, j + 1))

            # header
            f.write("# vtk DataFile Version 3.0\n")
            f.write("Unstructured Grid example\n")
            f.write("ASCII\nDATASET UNSTRUCTURED_GRID\n")
            f.write("POINTS " + point_number + " float\n")
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
                    f.write("9\n")

    def output_vtk_tri(fname, path):
        fname = path + "tri_" + fname + ".vtk"
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
                    tri1 = "3 " + str(pid(i, j)) + " "  + str(pid(ip1, j + 1)) + " " + str(pid(i, j + 1))
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

    def output_bdm(fname, path):
        fname = path + fname + ".mayugrid2"
        def header_gen(area_name, point_number):
            head1 = "begin " + area_name + "\n"
            head2 = point_number + "\n"
            return head1 + head2

        with open(fname, 'w') as f:
            f.write("type ""2D""\n")
            f.write(header_gen(str(fname), str(xi_max * eta_max)))
            for j in range(1, eta_max):
                for i in range(xi_max):
                    f.write(str(grid_x[i, j]) + " " + str(grid_y[i, j]) + "\n")
            f.write("end " + str(fname))
            f.write(header_gen("object_surface", str(xi_max)))
            for i in range(xi_max):
                f.write(str(grid_x[i, 0]) + " " + str(grid_y[i, 0]) + "\n")
            f.write("end object_surface")

    def plot_tmp(detail = True, zoom = False, overview = True):
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
    def offset_surface(z, outer = False, max_incremental=0.1, accel=1.0):
        size = z.shape[0]

        delta = np.zeros(size, dtype = complex)
        delta[0] = z[1] - z[size - 1]
        for i in range(1, size - 1):
            delta[i] = z[i + 1] - z[i - 1]
        delta[size - 1] = z[0] - z[size - 2]

        normal = -1j * delta / np.abs(delta)
        # 格子の裏返り防止(隣接する格子点から外向きξ方向へ伸びる2つの格子線がなす角度が90°以上開いている場合に，新しいη格子線が既存のη格子線と被らないように角度を修正する)
        def prevent_turn_over_cell(i, imp1, downwind=True):
            dot_normal = np.real(normal[i] * np.conj(normal[imp1])) # 隣接する格子線同士の内積を取る
            flag = 0
            if dot_normal < 0:  # 内積が負のとき修正対象
                flag = 1
                angle1 = np.angle(normal[i])
                angle0 = np.angle(normal[imp1])
                if angle1 < angle0:
                    angle1 += 2.0 * np.pi
                if downwind:    # ξが小さい側を修正する場合
                    normal[i] -= (angle1 - angle0) / 100 * 5  # 角度差の5%だけ寄せる
                else:   # ξが大きい側を変更する場合(風下と風上とを2回セットとして修正する)
                    normal[i] += (angle1 - angle0) / 100 * 5  # 角度差の5%だけ寄せる
            return flag
        
        flag = 1
        while flag != 0:    # flag=0となるまで格子の裏返り防止処理を掛ける
            flag = 0
            flag += prevent_turn_over_cell(0, xi_max - 1)
            for i in range(1, size):
                flag += prevent_turn_over_cell(i, i - 1)
            for i in range(size-1):
                flag += prevent_turn_over_cell(i, i + 1, downwind=False)
            flag += prevent_turn_over_cell(xi_max - 1, 0, downwind=False)
        
        incremental = accel * min(min(2.0 / np.pi * np.min(np.abs(delta)), np.average(np.abs(delta))), max_incremental)  # obj sizeとboundaryサイズを均等に分割したときの幅で置換すべき(0.1)

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

    def jacobi_pre(grid_x, grid_y):
        tmp_grid_x = grid_x.copy()
        tmp_grid_y = grid_y.copy()
        xi_line = []    # [0, int(eta_max/2)]
        eta_line = [0]
        flag = 0
        while flag == 0:
            flag = 1
            control_P, control_Q = update_control_function(xi_line, eta_line)
            for j in range(eta_max-2, 0, -1):
                for i in range(xi_max-1, -1, -1):
                    g11ij = g11(i, j)
                    g22ij = g22(i, j)
                    g12ij = g12(i, j)
                    im2, im1, ip1, ip2 = get_im2_im1_ip1_and_ip2(i)
                    scale = 0.5 / (g11ij + g22ij)
                    term1_x = g22ij * (sum_i_x2(im2, im1, ip1, ip2, j))
                    term1_y = g22ij * (sum_i_y2(im2, im1, ip1, ip2, j))

                    term2_x = g11ij * (sum_j_x2(i, j))
                    term2_y = g11ij * (sum_j_y2(i, j))

                    term3_x = -0.5 * g12ij * (sum_i_j_x(im1, ip1, j))
                    term3_y = -0.5 * g12ij * (sum_i_j_y(im1, ip1, j))

                    term4_x = 0.5 * g22ij * control_P[i, j] * (grid_x[ip1, j] - grid_x[im1, j])
                    term4_y = 0.5 * g22ij * control_P[i, j] * (grid_y[ip1, j] - grid_y[im1, j])

                    term5_x = 0.5 * g11ij * control_Q[i, j] * (grid_x[i, j + 1] - grid_x[i, j - 1])
                    term5_y = 0.5 * g11ij * control_Q[i, j] * (grid_y[i, j + 1] - grid_y[i, j - 1])

                    tmp_grid_x[i, j] = scale * (term1_x + term2_x + term3_x + term4_x + term5_x)
                    tmp_grid_y[i, j] = scale * (term1_y + term2_y + term3_y + term4_y + term5_y)

            res = np.sum((grid_x - tmp_grid_x)**2 + (grid_y - tmp_grid_y)**2)
            if res > 0.0000001:
                flag = 0
            # print(res)
            grid_x = tmp_grid_x.copy()
            grid_y = tmp_grid_y.copy()
        return grid_x, grid_y

    def get_equidistant_curve(z2):
        t, total_len = get_length_rate(z2, output_total_length=True)
        fx = interpolate.PchipInterpolator(np.hstack((t, np.array([1.0]))), np.real(np.hstack((z2, z2[0]))))
        fy = interpolate.PchipInterpolator(np.hstack((t, np.array([1.0]))), np.imag(np.hstack((z2, z2[0]))))
        equidistant_t = np.linspace(0, 1, z2.shape[0] + 1)[:z2.shape[0]]
        return fx(equidistant_t) + 1j * fy(equidistant_t)

    def equidistant_offset(z2, max_incremental, accel):
        z2 = offset_surface(z2, outer=True, max_incremental=max_incremental, accel = accel)
        return get_equidistant_curve(z2)

    get_object_center = lambda z2: np.average(np.real(z2)) + 1j * np.average(np.imag(z2))
    get_length = lambda x: np.max(x) - np.min(x)
    get_model_length = lambda z: max(get_length(np.real(z)), get_length(np.imag(z)))

    magnification = 5.0
    max_incremental = 10.0   # radius1 * (magnification - 1) / eta_max
    accel_parameter = 1.3   # 物体遠方領域でオフセット量を増やす際の割合
    # accel量の設定
    def set_accel(j, accel_parameter):
        if j < int(eta_max / 2):
            return 1.0
        else:
            return accel_parameter
        
    z1_eq = get_equidistant_curve(z1)

    grid_x[:, 0] = np.real(z1_eq)
    grid_y[:, 0] = np.imag(z1_eq)
    z2 = equidistant_offset(z1_eq, max_incremental, accel = 1.0)

    # z2 = np.hstack((z2[1:], z2[0]))
    circular_grid = False
    if circular_grid:
        for j in range(1, eta_max-1):
            grid_x[:, j] = np.real(z2)
            grid_y[:, j] = np.imag(z2)
            z2 = equidistant_offset(z2, max_incremental, accel = 1.0)

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
            mix_rate = 1.0 - max(np.exp(-0.18 * float(j - 1)), 0.7)   # j番目のη線と直交するように与えたz2_orthogonalと，Δηが均等になるように与えたz2_equidistantの配合比率
            while flag >= 0:
                flag += 1
                z2_equidistant = equidistant_offset(z2, max_incremental, set_accel(j, accel_parameter))
                z2_orthogonal = offset_surface(z2, outer=True, max_incremental=max_incremental, accel = set_accel(j, accel_parameter))
                fix_z2 = (1.0 - mix_rate) * z2_orthogonal + mix_rate * z2_equidistant
                delta_j__ = np.hstack((z2[1:] - z2[:xi_max - 1], z2[0] - z2[xi_max - 1]))
                delta_jp1 = np.hstack((fix_z2[1:] - fix_z2[:xi_max - 1], fix_z2[0] - fix_z2[xi_max - 1]))
                if np.any(np.real(delta_j__*np.conj(delta_jp1)) < 0):
                    mix_rate = 1.0 - max(-0.18 * np.exp(-float(j - 2)), 0.7 - 0.1 * flag)
                    print(mix_rate)
                    if mix_rate > 0.9:
                        print("not converge")
                        exit()
                else:
                    flag = -1
                    z2 = fix_z2

            if (np.max(np.real(z2)) - np.min(np.real(z2)) > 40.0 * model_length):
                eta_max = j + 2
                break

        z3 = z2

    grid_x[:, eta_max-1] = np.real(z3)
    grid_y[:, eta_max - 1] = np.imag(z3)
    # grid_x, grid_y = jacobi_pre(grid_x, grid_y)
    grid_x, grid_y = grid_x[:, :eta_max], grid_y[:, :eta_max]
    if vtk:
        output_vtk(fname, path)
    if trianglation:
        output_vtk_tri(fname, path)
    if bdm == True:
        output_bdm(fname, path)
    grid_x = np.vstack((grid_x[:, :], grid_x[0, :]))
    grid_y = np.vstack((grid_y[:, :], grid_y[0, :]))

    # for RANS simulation   # get distance from wall
    def get_distance_from_object(grid_x, grid_y):
        wall2point = [[] for i in range(xi_max)]  # [i][n]    i:wall number, n:point number (in this list)
        point2wall = np.zeros((xi_max, eta_max, 2))  # (i, j, 0): wall number, (i, j, 1): distance between wall to point
        # edge_i is made by point_i and point_i+1
        make_e_angle_from_p_angle = lambda p_angle: 0.5 * np.hstack(
            ((p_angle[1:] + p_angle[:size - 1]), (p_angle[0] - p_angle[size - 1])))
        set_0to2pi = lambda e_angle: np.where(e_angle <= -np.pi, e_angle + 2.0 * np.pi,
                                              e_angle)  # [-pi, pi] �� [0, 2pi]

        eta0 = grid_x[:, 0] + 1j * grid_y[:, 0]  # object_surface (clock-wise)
        center = get_object_center(eta0)
        p_angle = np.angle(eta0 - center)  # grid_point's angle
        size = p_angle.shape[0]

        e_angle = make_e_angle_from_p_angle(p_angle)
        e_angle = set_0to2pi(e_angle)
        e_number = np.arange(size)  # clock wise
        # if (np.any(e_angle[1:] - e_angle[:size - 1] < 0.0) or (e_angle[0] - e_angle[size - 1] < 0.0)):  # if object center is outside of object
        e_number = e_number[np.argsort(e_angle)]
        e_angle = np.sort(e_angle).tolist()

        def get_candidate_point(nearest_edge):
            if (nearest_edge == 0) or (nearest_edge == size - 1):  # edge made by point0 and point1
                return [size - 2, 0, 1, 2]  # candidate_point 0, 1, 2, 3
            elif nearest_edge == size - 2:
                return [size - 3, size - 2, 0, 1]
            elif nearest_edge == size - 3:
                return [size - 4, size - 3, size - 2, 0]
            else:
                return [nearest_edge - 1, nearest_edge, nearest_edge + 1, nearest_edge + 2]

        get_SQ_distance_from_point = lambda i, j, cp: (grid_x[i, j] - grid_x[cp, 0]) ** 2 + (
                    grid_y[i, j] - grid_y[cp, 0]) ** 2

        for j in range(1, eta_max - 1):
            eta_jth_line = grid_x[:, j] + 1j * grid_y[:, j]  # jth eta line (complex form)
            # jth_p_angle = np.angle(eta_jth_line - center)
            jth_e_angle = set_0to2pi(
                make_e_angle_from_p_angle(np.angle(eta_jth_line - center)))  # jth edge angle (i=0, 1, ..., xi_max-1)

            for i in range(xi_max):
                tmp_nearest_edge_candidate = bisect.bisect_left(e_angle, jth_e_angle[i])  # nearest_edge from jth edge (edge_number)
                if tmp_nearest_edge_candidate == size - 1:
                    tmp_nearest_edge_candidate = 0
                nearest_edge_candidate = e_number[tmp_nearest_edge_candidate]
                candidate_point = get_candidate_point(nearest_edge_candidate)  # adjacent point from nearest edge (point number)
                distance = np.zeros(4)
                for k in range(4):
                    distance[k] = get_SQ_distance_from_point(i, j, candidate_point[k])  # get distance ** 2
                min_k = np.argmin(distance)
                point2wall[i, j, 0] = candidate_point[min_k]  # nearest point number = nearest_edge_number  # edge = wall
                point2wall[i, j, 1] = sqrt(distance[min_k])  # nearest_distance
                wall2point[candidate_point[min_k]].append([i, j])

        return point2wall, wall2point
    
    point2wall, wall2point = get_distance_from_object(grid_x, grid_y)
    if mg2:
        write_out_mayugrid2(fname, path, grid_x, grid_y, point2wall, wall2point)
    return grid_x, grid_y, point2wall, wall2point


def write_out_mayugrid2(fname, path, grid_x, grid_y, point2wall, wall2point):
    fname = path + fname + ".mayugrid2"
    xi_max = grid_x.shape[0]
    eta_max = grid_y.shape[1]
    with open(fname, 'w') as f:
        f.write("# mayugrid2 Version 0.9\n")
        f.write("xi_max " + str(xi_max) + "\n")
        f.write("eta_max " + str(eta_max) + "\n")
        f.write("xi = i, eta = j, N = i + xi_max * j, and include loopy elements.\n")
        f.write("x - y coordinate " + str(xi_max * eta_max) + " items inline\n")
        # point coordinates
        for j in range(eta_max):
            for i in range(xi_max):
                f.write(str(grid_x[i, j]) + " " + str(grid_y[i, j]) + "\n")

        f.write("The edge number of the object surface nearest to each grid point and the distance\n".replace(" ", "_"))
        for j in range(eta_max):
            for i in range(xi_max - 1):
                f.write(str(point2wall[i, j, 0]) + " " + str(point2wall[i, j, 1]) + "\n")

        f.write("Grid point numbers belonging to each edge of the object surface\n")
        for edge in range(xi_max - 1):
            p_size = len(wall2point[edge])
            f.write("Total_number_of_points_belonging_to_" + str(edge) + "th_edge: " + str(p_size) + "\n")
            for p in range(p_size):
                f.write(str(wall2point[edge][p]).replace("[", "").replace(",", "").replace("]", "") + "\n")


def make_grid(fname, type, size=100, naca4="0012", center_x=0.08, center_y=0.08, mayugrid2=True, vtk=True, bdm=True, trianglation=True, path = ""):
    z1, size = get_complex_coords(type = type, center_x = center_x, center_y=center_y, naca4 = naca4, size = size)
    z1 = deduplication(z1)[::-1]
    make_grid_seko(z1, path, fname, mayugrid2, vtk, bdm, trianglation)
    
def main():
    z1, size = get_complex_coords(type = 3, naca4 = "4912b", size = 100)
    # z1, size = get_complex_coords(type=1, center_x=0.08, center_y=0.3, naca4="4912", size=100)
    # plt.plot(np.real(z1), np.imag(z1))
    z1 = deduplication(z1)[::-1]
    make_grid_seko(z1)
    # plt.plot(np.real(z1), np.imag(z1))
    # plt.show()

if __name__ == '__main__':
    main()