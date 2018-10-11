# coding: utf-8
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import bicg, spsolve
from math import floor
import numpy as np
from naca_4digit_test import Naca_4_digit, Naca_5_digit
from joukowski_wing import joukowski_wing_complex, karman_trefftz_wing_complex
import matplotlib.pyplot as plt

# 物体表面の複素座標を取得する
def get_complex_coords(type, size, center_x = -0.08, center_y = 0.08, naca4 = "0012"):
    def reshape_z(z):
        if z[0] != z[z.shape[0] - 1]:
            return np.concatenate([z, z[0].reshape(-1)]), z.shape[0] + 1
        else:
            return z, z.shape[0]
    
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
    
    if type < 3:
        return reshape_z(z)
    else:
        return z, z.shape[0]



# 格子の外部境界(分割数は物体と同じの，物体形状の長軸長さのmagnification倍の円を返す)
def get_outer_boundary(z1, magnification=5):
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
    resolution = z1.shape[0]
    center = center_x + 1j * center_y
    radius = model_length * magnification
    t = np.linspace(start=0, stop=2.0 * np.pi, num=resolution)
    z3 = center + radius * np.exp(1j * t)   # Counter Clock Wise
    return z3[::-1] # Clock Wise and -1

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
# 内積の総和が小さければ直交性高い
def check_orthogonal(grid_x, grid_y):
    pass

def make_grid_seko(z1, z2, z3, z4):
    xi_max = z1.shape[0]  # xi方向の格子点数
    eta_max = z2.shape[0]  # eta方向の格子点数
    center_x = 0.25 * (np.average(np.real(z1)) + np.average(np.real(z2)) + np.average(np.real(z3)) + np.average(np.real(z4)))
    center_y = 0.25 * (np.average(np.imag(z1)) + np.average(np.imag(z2)) + np.average(np.imag(z3)) + np.average(np.imag(z4)))
    grid_x = np.zeros((xi_max, eta_max))  # x座標格納
    grid_y = np.zeros((xi_max, eta_max))  # y座標格納
    for i in range(eta_max):
        grid_x[i, :] = np.linspace(np.real(z1[i]), np.real(z3[eta_max-1-i]), eta_max)
        grid_y[i, :] = np.linspace(np.imag(z1[i]), np.imag(z3[eta_max-1-i]), eta_max)
    for i in range(xi_max):
        plt.plot(grid_x[i, :], grid_y[i, :])
        plt.plot(grid_x[:, i], grid_y[:, i])

    plt.show()
    # 境界条件適用
    grid_x[0, :] = np.real(z1)  # 底辺
    grid_x[:, 0] = np.real(z2)  # 左辺
    grid_x[xi_max - 1, :] = np.real(z3)  # 上辺
    grid_x[:, eta_max - 1] = np.real(z4)

    grid_y[0, :] = np.imag(z1)  # 底辺
    grid_y[:, 0] = np.imag(z2)[::-1]  # 左辺
    grid_y[xi_max - 1, :] = np.imag(z3) # 上辺
    grid_y[:, eta_max - 1] = np.imag(z4)    # 右辺

    # 以降係数は無視し，足し合わせる際に帳尻を合わせる
    x_xi = lambda i, j: (grid_x[i + 1, j] - grid_x[i - 1, j])   # xのξ微分(の2倍)
    x_eta = lambda i, j: (grid_x[i, j + 1] - grid_x[i, j - 1])  # xのη微分(の2倍)
    y_xi = lambda i, j: (grid_y[i + 1, j] - grid_y[i - 1, j])   # yのξ微分(の2倍)
    y_eta = lambda i, j:(grid_y[i, j + 1] - grid_y[i, j - 1])   # yのη微分(の2倍)
    
    Aij = lambda i, j: (x_eta(i, j) ** 2 + y_eta(i, j) ** 2)    # ξ2階微分の係数(の4倍)
    Bij = lambda i, j: (x_xi(i, j) * x_eta(i, j) + y_xi(i, j) * y_eta(i, j))    # ξη交差微分の係数(の4倍)
    Cij = lambda i, j: (x_xi(i, j) ** 2 + y_xi(i, j) ** 2)    # η2階微分の係数(の4倍)
    
    x_xixi = lambda i, j: grid_x[i + 1, j] - 2.0 * grid_x[i, j] + grid_x[i - 1, j]  # xのξ2階微分
    y_xixi = lambda i, j: grid_y[i + 1, j] - 2.0 * grid_y[i, j] + grid_y[i - 1, j]  # yのξ2階微分
    
    x_etaeta = lambda i, j: grid_x[i, j + 1] - 2.0 * grid_x[i, j] + grid_x[i, j - 1]    # xのη2階微分
    y_etaeta = lambda i, j: grid_y[i, j + 1] - 2.0 * grid_y[i, j] + grid_y[i, j - 1]    # yのη2階微分
    
    x_xieta = lambda i, j: (grid_x[i + 1, j + 1] - grid_x[i + 1, j - 1] - grid_x[i - 1, j + 1] + grid_x[i - 1, j - 1])  # xのξη交差微分(の4倍)
    y_xieta = lambda i, j: (grid_y[i + 1, j + 1] - grid_y[i + 1, j - 1] - grid_y[i - 1, j + 1] + grid_y[i - 1, j - 1])  # yのξη交差微分(の4倍)

    
    delta_t = 0.0001  # delta-formの時間発展用
    point_number = (xi_max - 2) * (eta_max - 2) # 未知数の総数
    Nj = eta_max - 2    # 係数行列の列方向ブロックサイズ
    
    # [I,I]からdeltaだけずれた位置の成分が何個目のブロック行列に属するか返す関数
    block_id = lambda delta: floor((I + delta) / Nj)
    # xについての式 + yについての式
    for iter in range(1000):
        # matrix & rhsの準備
        # rhs = np.zeros(2 * point_number)    # rhsベクトルの初期化
        sol = np.zeros(2 * point_number)    # 解ベクトルの初期化
        # matrix = lil_matrix((2 * point_number, 2 * point_number))   # 係数行列matrixの初期化
        dij = 0
        for I in range(point_number):
            block = block_id(0) # 現在のブロック
            i = I - block * Nj + 1  # 未知数の0 = 1番目の格子点(0番目の格子点は更新する必要がないため放置) # matrixの0番はgrid_xyの1番に相当
            j = block + 1  # aij対策で+1
            
            aij = Aij(i, j) # I行内では常に同じ添え字ijを用いるため最初に計算しておく
            bij = Bij(i, j)
            cij = Cij(i, j)
            dij = max(max(aij, bij), max(cij, dij))
            """
            matrix[I, I] = 1.0 - delta_t * 2.0 * (aij + cij)    # xの係数
            matrix[I + point_number, I + point_number] = 1.0 - delta_t * 2.0 * (aij + cij)  # yの係数
            
            # AΔt
            if (I > 0) and (block_id(-1) == block):  # 対角成分と同一ブロック内にないときは0
                matrix[I, I - 1] = delta_t * aij
                matrix[I + point_number, I + point_number - 1] = delta_t * aij
            if (I < point_number - 1) and (block_id(1) == block):
                matrix[I, I + 1] = delta_t * aij
                matrix[I + point_number, I + point_number + 1] = delta_t * aij
            
            # CΔt
            if I > Nj - 1:
                matrix[I, I - Nj] = delta_t * cij
                matrix[I + point_number, I + point_number - Nj] = delta_t * cij
            if I < point_number - Nj:
                matrix[I, I + Nj] = delta_t * cij
                matrix[I + point_number, I + point_number + Nj] = delta_t * cij
            
            # -0.5*bij
            if ((I > (Nj - 1)) and (block_id(-(Nj - 1)) == block - 1)):
                matrix[I, I - (Nj - 1)] = - delta_t * 0.5 * bij
                matrix[I + point_number, I + point_number - (Nj - 1)] = - delta_t * 0.5 * bij
            if ((I < point_number - (Nj - 1)) and (block_id(Nj - 1) == block + 1)):
                matrix[I, I + (Nj - 1)] = - delta_t * 0.5 * bij
                matrix[I + point_number, I + point_number + (Nj - 1)] = - delta_t * 0.5 * bij
            
            # +0.5*bij
            if ((I > (Nj + 1) - 1 and (block_id(-(Nj + 1)) == block - 1))):
                matrix[I, I - (Nj + 1)] = delta_t * 0.5 * bij
                matrix[I + point_number, I + point_number - (Nj + 1)] = delta_t * 0.5 * bij
            if ((I < point_number - (Nj + 1) - 1) and (block_id(Nj + 1) == block + 1)):
                matrix[I, I + (Nj + 1)] = delta_t * 0.5 * bij
                matrix[I + point_number, I + point_number + (Nj + 1)] = delta_t * 0.5 * bij
            #"""
            # rhs[I] = - delta_t * (aij * x_xixi(i, j) - 0.5 * bij * x_xieta(i, j) + cij * x_etaeta(i, j))
            # rhs[I + point_number] = -delta_t * (aij * y_xixi(i, j) - 0.5 * bij * y_xieta(i, j) + cij * y_etaeta(i, j))
            sol[I] = (aij * x_xixi(i, j) - 0.5 * bij * x_xieta(i, j) + cij * x_etaeta(i, j))
            sol[I + point_number] = (aij * y_xixi(i, j) - 0.5 * bij * y_xieta(i, j) + cij * y_etaeta(i, j))
        
        # delta_xy = bicg(matrix, rhs)[0]
        # matrix = matrix.tocsr()
        # delta_xy = spsolve(matrix, rhs)

        delta_t = 1.0 / (3.0 * dij)
        # grid_x[1:xi_max - 1, 1:eta_max - 1] += delta_xy[:point_number].reshape(xi_max - 2, -1).T
        # grid_y[1:xi_max - 1, 1:eta_max - 1] += delta_xy[point_number:].reshape(xi_max - 2, -1).T
        grid_x[1:xi_max - 1, 1:eta_max - 1] = delta_t * sol[:point_number].reshape(xi_max - 2, -1).T + grid_x[1:xi_max - 1, 1:eta_max - 1]
        grid_y[1:xi_max - 1, 1:eta_max - 1] = delta_t * sol[point_number:].reshape(xi_max - 2, -1).T + grid_y[1:xi_max - 1, 1:eta_max - 1]
        # maximum = max(np.max(np.abs(grid_x)), np.max(np.abs(grid_y)))
        # print(maximum)
        if iter % 50 == 0:
            for i in range(xi_max):
                plt.plot(grid_x[i, :], grid_y[i, :])
                plt.plot(grid_x[:, i], grid_y[:, i])

            plt.show()


# return grid_x, grid_y


# a, b, cからなる三重対角行列(n×n)および，n次元定数ベクトルdの連立方程式の解xを返す
# aは2行目からb行目まで，cは1行目からn-1行まで
def tridiagonal_matrix_algorithm(a, b, c, d, n):
    c_prime = np.zeros(n)
    d_prime = np.zeros(n)
    x = np.zeros(n)
    c_prime[0] = c[0] / b[0]
    d_prime[0] = d[0] / b[0]

    for i in range(1, n - 1):
        c_prime[i] = c[i] / (b[i] - a[i] * c_prime[i - 1])
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / (b[i] - a[i] * c_prime[i - 1])

    d_prime[n-1] = (d[n-1] - a[n-1] * d_prime[n-2]) / (b[n-1] - a[n-1] * c_prime[n-2])

    x[n-1] = d_prime[n-1]
    for i in range(n-2, -1, -1):
        x[i] = d_prime[i] - c_prime[i] * x[i + 1]

    return x

def check_TDMA():
    N = 5
    A = np.zeros((N, N))
    A[0, 0] = 3
    A[0, 1] = 4
    for i in range(1, N-1):
        A[i, i-1] = 2+i
        A[i, i] = 3
        A[i, i+1] = 4
    A[N-1, N-2] = 2+N-1
    A[N-1, N-1] = 3
    print(A)

    b = np.arange(5)
    print(linalg.det(A))
    print(b)
    print(linalg.solve(A, b))
    left = np.arange(5) + 3 - 1
    diag = np.ones(5) * 3
    right = np.ones(5) * 4

    print(tridiagonal_matrix_algorithm(left, diag, right, b, 5))
    exit()

def main():
    z1, size = get_complex_coords(type = 3, naca4 = "4912", size = 5)
    z3 = get_outer_boundary(z1, magnification=10)
    z2 = get_connect_z1_to_z3(z1, z3)
    z4 = z2
    """
    box = 10
    z1 = np.zeros(box) + 1j * np.arange(box)    #
    z2 = np.arange(box) + 1j * np.zeros(box)
    z3 = np.ones(box) * (box-1) + 1j * np.arange(box)
    z4 = np.arange(box) + 1j * np.ones(box) * (box-1)
    #"""
    plt.plot(np.real(z1), np.imag(z1))
    plt.plot(np.real(z2), np.imag(z2), "o")
    plt.plot(np.real(z3), np.imag(z3), "x")
    plt.plot(np.real(z4), np.imag(z4))
    plt.show()
    
    make_grid_seko(z1, z2, z3, z4)
    plt.plot(np.real(z1), np.imag(z1))
    plt.plot(np.real(z2), np.imag(z2), "o")
    plt.plot(np.real(z3), np.imag(z3), "x")
    plt.show()

    exit()
    """
    grid[0] = line_up_z2xy(z1[:size-1], two_rows = True)
    
    xy0 = line_up_z2xy(offset_surface(z1)[:size-1])
    
    
    # print(xy0)
    phi0 = np.zeros(size-1)
    phi1 = np.zeros(size-1)
    V0 = np.zeros(size-1) + 0.0001
    V1 = np.zeros(size-1)
    cM, cV = made_coefficient_matrix(xy0[0::2], grid[0, :, 0], xy0[1::2], grid[0, :, 1], phi0, phi1, V0, V1)

    # print(cM)
    xy1 = linalg.solve(cM, cV)
    xy1_ = spla.bicg(cM, cV)[0]
    print(xy1_)
    plt.plot(np.real(z1), np.imag(z1))
    plt.plot(xy1_[0::2], xy1_[1::2])
    plt.show()
    """
if __name__ == '__main__':
    main()

"""
def make_grid_seko(z1, z2, z3, z4):
    xi_max = z1.shape[0]  # xi方向の格子点数
    eta_max = z2.shape[0]  # eta方向の格子点数
    
    grid_x = np.zeros((xi_max, eta_max))  # x座標格納
    grid_y = np.zeros((xi_max, eta_max))  # y座標格納
    # 境界条件適用
    grid_x[0, :] = np.real(z1)  # 底辺
    grid_x[:, 0] = np.real(z2)  # 左辺
    grid_x[xi_max - 1, :] = np.real(z3)  # 上辺
    grid_x[:, eta_max - 1] = np.real(z4)

    grid_y[0, :] = np.imag(z1)  # 底辺
    grid_y[:, 0] = np.imag(z2)[::-1]  # 左辺
    grid_y[xi_max - 1, :] = np.imag(z3) # 上辺
    grid_y[:, eta_max - 1] = np.imag(z4)    # 右辺

    # 以降係数は無視し，足し合わせる際に帳尻を合わせる
    x_xi = lambda i, j: (grid_x[i + 1, j] - grid_x[i - 1, j])   # xのξ微分
    x_eta = lambda i, j: (grid_x[i, j + 1] - grid_x[i, j - 1])  # xのη微分
    y_xi = lambda i, j: (grid_y[i + 1, j] - grid_y[i - 1, j])   # yのξ微分
    y_eta = lambda i, j:(grid_y[i, j + 1] - grid_y[i, j - 1])   # yのη微分
    
    Aij = lambda i, j: (x_eta(i, j) ** 2 + y_eta(i, j) ** 2)    # ξ2階微分の係数
    Bij = lambda i, j: (x_xi(i, j) * x_eta(i, j) + y_xi(i, j) * y_eta(i, j))    # ξη交差微分の係数
    Cij = lambda i, j: (x_xi(i, j) ** 2 + y_xi(i, j) ** 2)    # η2階微分の係数
    
    x_xixi = lambda i, j: grid_x[i + 1, j] - 2.0 * grid_x[i, j] + grid_x[i - 1, j]  # xのξ2階微分
    y_xixi = lambda i, j: grid_y[i + 1, j] - 2.0 * grid_y[i, j] + grid_y[i - 1, j]  # yのξ2階微分
    
    x_etaeta = lambda i, j: grid_x[i, j + 1] - 2.0 * grid_x[i, j] + grid_x[i, j - 1]    # xのη2階微分
    y_etaeta = lambda i, j: grid_y[i, j + 1] - 2.0 * grid_y[i, j] + grid_y[i, j - 1]    # yのη2階微分
    
    x_xieta = lambda i, j: (grid_x[i + 1, j + 1] - grid_x[i + 1, j - 1] - grid_x[i - 1, j + 1] + grid_x[i - 1, j - 1])  # xのξη交差微分
    y_xieta = lambda i, j: (grid_y[i + 1, j + 1] - grid_y[i + 1, j - 1] - grid_y[i - 1, j + 1] + grid_y[i - 1, j - 1])  # yのξη交差微分

    
    delta_t = 0.1  # delta-formの時間発展用
    point_number = (xi_max - 2) * (eta_max - 2) # 未知数の総数
    Nj = eta_max - 2    # 係数行列の列方向ブロックサイズ
    
    # [I,I]からdeltaだけずれた位置の成分が何個目のブロック行列に属するか返す関数
    block_id = lambda delta: floor((I + delta) / Nj)
    # xについての式 + yについての式
    for iter in range(1000):
        # matrix & rhsの準備
        rhs = np.zeros(2 * point_number)    # rhsベクトルの初期化
        matrix = lil_matrix((2 * point_number, 2 * point_number))   # 係数行列matrixの初期化
        for I in range(point_number):
            block = block_id(0) # 現在のブロック
            i = I - block * Nj + 1  # 未知数の0 = 1番目の格子点(0番目の格子点は更新する必要がないため放置) # matrixの0番はgrid_xyの1番に相当
            j = block + 1  # aij対策で+1
            
            aij = Aij(i, j) # I行内では常に同じ添え字ijを用いるため最初に計算しておく
            bij = Bij(i, j)
            cij = Cij(i, j)
            
            matrix[I, I] = 1.0 - delta_t * 2.0 * (aij + cij)    # xの係数
            matrix[I + point_number, I + point_number] = 1.0 - delta_t * 2.0 * (aij + cij)  # yの係数
            
            # AΔt
            if (I > 0) and (block_id(-1) == block):  # 対角成分と同一ブロック内にないときは0
                matrix[I, I - 1] = delta_t * aij
                matrix[I + point_number, I + point_number - 1] = delta_t * aij
            if (I < point_number - 1) and (block_id(1) == block):
                matrix[I, I + 1] = delta_t * aij
                matrix[I + point_number, I + point_number + 1] = delta_t * aij
            
            # CΔt
            if I > Nj - 1:
                matrix[I, I - Nj] = delta_t * cij
                matrix[I + point_number, I + point_number - Nj] = delta_t * cij
            if I < point_number - Nj:
                matrix[I, I + Nj] = delta_t * cij
                matrix[I + point_number, I + point_number + Nj] = delta_t * cij
            
            # -0.5*bij
            if ((I > (Nj - 1)) and (block_id(-(Nj - 1)) == block - 1)):
                matrix[I, I - (Nj - 1)] = - delta_t * 0.5 * bij
                matrix[I + point_number, I + point_number - (Nj - 1)] = - delta_t * 0.5 * bij
            if ((I < point_number - (Nj - 1)) and (block_id(Nj - 1) == block + 1)):
                matrix[I, I + (Nj - 1)] = - delta_t * 0.5 * bij
                matrix[I + point_number, I + point_number + (Nj - 1)] = - delta_t * 0.5 * bij
            
            # +0.5*bij
            if ((I > (Nj + 1) - 1 and (block_id(-(Nj + 1)) == block - 1))):
                matrix[I, I - (Nj + 1)] = delta_t * 0.5 * bij
                matrix[I + point_number, I + point_number - (Nj + 1)] = delta_t * 0.5 * bij
            if ((I < point_number - (Nj + 1) - 1) and (block_id(Nj + 1) == block + 1)):
                matrix[I, I + (Nj + 1)] = delta_t * 0.5 * bij
                matrix[I + point_number, I + point_number + (Nj + 1)] = delta_t * 0.5 * bij

            rhs[I] = - delta_t * (aij * x_xixi(i, j) - 0.5 * bij * x_xieta(i, j) + cij * x_etaeta(i, j))
            rhs[I + point_number] = -delta_t * (aij * y_xixi(i, j) - 0.5 * bij * y_xieta(i, j) + cij * y_etaeta(i, j))
        
        delta_xy = bicg(matrix, rhs)[0]
        
        grid_x[1:xi_max - 1, 1:eta_max - 1] += delta_xy[:point_number].reshape(xi_max - 2, -1).T
        grid_y[1:xi_max - 1, 1:eta_max - 1] += delta_xy[point_number:].reshape(xi_max - 2, -1).T
        maximum = max(np.max(np.abs(grid_x)), np.max(np.abs(grid_y)))
        print(maximum)
        if maximum < 8:
            for i in range(xi_max):
                plt.plot(grid_x[i, :], grid_y[i, :])
            plt.show()
            
            for j in range(eta_max):
                plt.plot(grid_x[:, j], grid_y[:, j])
            plt.show()
"""

"""
def make_grid_seko(z1, z2, z3):
    # 極端に近い点を除去(内側を優先的に除去)
    def delete_near_point(z, eps=0.0001):
        flag = 1
        while flag == 1:
            flag = 0
            len = np.abs(z[1:] - z[:z.shape[0]-1])
            num = np.argmin(len)
            if len[num] < eps:
                z = np.hstack([z[:num], z[num+1:]])
                flag = 1
        return z
    
    # 終端の重複している座標点を除去
    def deduplication(z, eps=0.0001):
        z = delete_near_point(z)
        flag = 1
        while flag == 1:
            flag = 0
            if np.abs(z[0] - z[z.shape[0] - 1]) <= eps:
                z = z[:z.shape[0] - 1]
                flag = 1
        return z
    
    # 終端に重複する座標点を追加
    def add_duplicates(z, eps=0.0001):
        if np.abs(z[0] - z[z.shape[0] - 1]) > eps:
            z = np.hstack([z, z[0]])
        return z
    
    # 中心からの角度基準で座標点を並べ替える
    def sort_by_angle(z):
        center = np.average(np.real(z)) + 1j * np.average(np.imag(z))
        return z[np.argsort(np.angle(z - center))]
    
    def sort_for_cross_delete(z0, z1):
        def find_intersection_line_by_line(p1, p2, p3, p4):
            def det123(z1, z2, z3):
                # a = np.array([[1, 1, 1], [np.real(z1), np.real(z2), np.real(z3)], [np.imag(z1), np.imag(z2), np.imag(z3)]])
                return np.linalg.det(np.array(
                    [[1, 1, 1], [np.real(z1), np.real(z2), np.real(z3)], [np.imag(z1), np.imag(z2), np.imag(z3)]]))
            
            if det123(p1, p2, p3) * det123(p1, p2, p4) < 0:
                if det123(p3, p4, p1) * det123(p3, p4, p2) < 0:
                    return True
            
            return False
        
        def swap_a_b(a, b):
            return b, a
        
        def half_length(dz, z):
            return 0.5 * (dz + z)
        
        size = z0.shape[0]
        flag = 1
        while flag == 1:
            flag = 0
            for i in range(size-1):
                if find_intersection_line_by_line(z0[i], z0[i+1], z1[i], z1[i+1]):
                    # z1[i], z1[i + 1] = swap_a_b(z1[i], z1[i+1])
                    z1[i] = half_length(z1[i], z0[i])
                    z1[i+1] = half_length(z1[i+1], z0[i+1])
                    flag = 1
            if find_intersection_line_by_line(z0[size-1], z0[0], z1[size-1], z1[0]):
                # z1[size-1], z1[0] = swap_a_b(z1[size-1], z1[0])
                z1[size-1] = half_length(z1[size-1], z0[0])
                z1[0] = half_length(z1[0], z0[size-1])
                flag = 1
                
        return z1
        
    # 物体表面を押し出す
    def offset_surface(z):
        def get_angle(z2, z1, z0):
            # if (np.imag(z2 - z1) >= 0) and (np.imag(z1 - z0)):
            az21 = np.angle(z2 - z1)
            az01 = np.angle(z0 - z1)
            print(az21, az01)
            #if az21 < 0:
            #    return 0.5 * (az21 + az01) + np.pi
            return 0.5 * (az21 + az01)

        size = z.shape[0]
        
        delta = np.zeros(size, dtype = complex)
        delta[0] = z[1] - z[size -1]
        for i in range(1, size-1):
            delta[i] = z[i + 1] - z[i - 1]
        delta[size - 1] = z[0] - z[size - 2]

        normal = -1j * delta / np.abs(delta)
        incremental = np.min(np.abs(delta))
        dz = z + normal * incremental
        center = np.average(np.real(z)) + 1j * np.average(np.imag(z))
        # np.argsort(dz - center)
        dz = sort_for_cross_delete(z, dz)
        return z + normal * incremental, incremental
    
    z1 = deduplication(z1)
    # z1 = sort_by_angle(z1)
    
    size = z1.shape[0]
    z = np.zeros((size, size), dtype = complex)
    z[:, 0] = z1
    distance = 0
    plt.plot(np.real(z[:, 0]), np.imag(z[:, 0]))
    for j in range(1, size):
        z[:, j], incremental = offset_surface(z[:, j - 1])
        distance += incremental
        # plt.plot(np.real(z[:, j]), np.imag(z[:, j]))
        
    # for i in range(size):
        # plt.plot(np.real(z[i, :]), np.imag(z[i, :]))
    
    
    print(distance)
    # plt.show()

    for i in range(int(size / 7)):
        plt.plot(np.real(z[i, :]), np.imag(z[i, :]))
        plt.plot(np.real(z[:, i]), np.imag(z[:, i]))
    xbot = 0.8
    xtop = 1.0
    ybot = 0.3
    ytop = 0.5
    plt.xlim(xbot, xtop)
    plt.ylim(ybot, ytop)
    plt.show()
    exit()
"""
"""
    xi_max = z1.shape[0]    # xi方向の格子点数
    eta_max = z2.shape[0]   # eta方向の格子点数

    grid_x = np.zeros((xi_max, eta_max))    # x座標格納
    grid_y = np.zeros((xi_max, eta_max))    # y座標格納

    center = np.average(np.real(z1))
    eps = 0.01
    min_rad = np.max(np.abs(z1 + center)) + eps
    max_rad = np.min(np.abs(z3 + center)) - eps
    radius = np.linspace(min_rad, max_rad, xi_max)

    for i in range(1, xi_max-1):
        circle = center + radius[i-1] * np.exp(1j*np.linspace(0, 2.0*np.pi, eta_max))
        grid_x[i, :] = np.real(circle)
        grid_y[i, :] = np.imag(circle)
    z3 = z3[::-1]

    # 境界条件適用
    grid_x[:, 0] = np.real(z2)  # 左辺
    grid_x[0, :] = np.real(z1)  # 底辺
    grid_x[xi_max-1, :] = np.real(z3)    # 上辺
    grid_x[:, eta_max - 1] = np.real(z2)

    grid_y[0, :] = np.imag(z1)  # 底辺
    grid_y[:, 0] = np.imag(z2) + eps    # 左辺
    grid_y[xi_max - 1, :] = np.imag(z3)
    grid_y[:, eta_max - 1] = np.imag(z2) - eps


    for i in range(xi_max):
        plt.plot(grid_x[i, :], grid_y[i, :])
    plt.show()
    for j in range(eta_max):
        plt.plot(grid_x[:, j], grid_y[:, j])
    plt.show()

    diff_x_xi = lambda i, j: 0.5 * (grid_x[i+1, j] - grid_x[i-1, j])
    diff_y_xi = lambda i, j: 0.5 * (grid_y[i+1, j] - grid_y[i-1, j])
    diff_x_eta = lambda i, j: 0.5 * (grid_x[i, j+1] - grid_x[i, j-1])
    diff_y_eta = lambda i, j: 0.5 * (grid_y[i, j+1] - grid_y[i, j-1])

    Aij = lambda i, j: (diff_x_eta(i, j)**2 + diff_y_eta(i, j)**2)
    Bij = lambda i, j: (diff_x_xi(i, j) * diff_x_eta(i, j) + diff_y_xi(i, j) * diff_y_eta(i, j))
    Cij = lambda i, j: (diff_x_xi(i, j)**2 + diff_y_xi(i, j)**2)

    diff_x_xixi = lambda i, j: grid_x[i+1, j] - 2.0 * grid_x[i, j] + grid_x[i-1, j]
    diff_y_xixi = lambda i, j: grid_y[i+1, j] - 2.0 * grid_y[i, j] + grid_y[i-1, j]

    diff_x_etaeta = lambda i, j: grid_x[i, j+1] - 2.0 * grid_x[i, j] + grid_x[i, j-1]
    diff_y_etaeta = lambda i, j: grid_y[i, j+1] - 2.0 * grid_y[i, j] + grid_y[i, j-1]

    diff_x_xieta = lambda i, j: grid_x[i+1, j+1] - grid_x[i+1, j-1] - grid_x[i-1, j+1] + grid_x[i-1, j-1]
    diff_y_xieta = lambda i, j: grid_y[i+1, j+1] - grid_y[i+1, j-1] - grid_y[i-1, j+1] + grid_y[i-1, j-1]

    delta_t = 0.01
    theta = 0.5
    delta_t *= theta
    point_number = (xi_max - 2) * (eta_max - 2)
    Ni = xi_max - 2
    Nj = eta_max - 2

    # [I,I]からdeltaだけずれた位置の成分が何列目のブロック行列に属するか返す
    block_id = lambda delta: floor((I + delta) / Nj)

    for iter in range(1000):
        # matrix & rhsの準備
        rhs = np.zeros(2 * point_number)
        matrix = lil_matrix((2 * point_number, 2 * point_number))
        for I in range(point_number):
            block = block_id(0)
            i = I - block * Nj + 1    # aij対策で+1  # matrixの0番はgrid_xyの1番に相当
            j = block + 1 # aij対策で+1

            aij = Aij(i, j)
            bij = Bij(i,j)
            cij = Cij(i, j)

            matrix[I, I] = 1.0 - delta_t * 2.0 * (aij + cij)
            matrix[I + point_number, I + point_number] = 1.0 - delta_t * 2.0 * (aij + cij)
            # x-direction
            if (I > 0) and (block_id(-1) == block): # 対角成分と同一ブロック内にないときは0
                matrix[I, I-1] = delta_t * aij
                matrix[I + point_number, I + point_number - 1] = delta_t * aij
            if (I < point_number-1) and (block_id(1) == block):
                matrix[I, I+1] = delta_t * aij
                matrix[I + point_number, I + point_number + 1] = delta_t * aij

            if I > Nj - 1:
                matrix[I, I-Nj] = delta_t * cij
                matrix[I + point_number, I + point_number - Nj] = delta_t * cij
            if I < point_number - Nj:
                matrix[I, I + Nj] = delta_t * cij
                matrix[I + point_number, I + point_number + Nj] = delta_t * cij

            if ((I > (Nj - 1)) and (block_id(-(Nj - 1)) == block-1)):
                matrix[I, I-(Nj-1)] = delta_t * 2.0 * bij
                matrix[I + point_number, I + point_number - (Nj - 1)] = delta_t * 2.0 * bij
            if ((I < point_number - (Nj - 1)) and (block_id(Nj - 1) == block + 1)):
                matrix[I, I+(Nj-1)] = delta_t * 2.0 * bij
                matrix[I + point_number, I + point_number + (Nj - 1)] = delta_t * 2.0 * bij


            if ((I > (Nj + 1) - 1 and (block_id(-(Nj+1)) == block - 1))):
                matrix[I, I-(Nj+1)] = - delta_t * 2.0 * bij
                matrix[I + point_number, I + point_number - (Nj + 1)] = -delta_t * 2.0 * bij
            if ((I < point_number - (Nj + 1) - 1) and (block_id(Nj + 1) == block + 1)):
                matrix[I, I+(Nj+1)] = -delta_t * 2.0 * bij
                matrix[I + point_number, I + point_number + (Nj + 1)] = - delta_t * 2.0 * bij
            if I==10:
                print(aij)
            rhs[I] = -delta_t * (aij * diff_x_xixi(i, j) - 2.0 * bij * diff_x_xieta(i, j) + cij * diff_x_etaeta(i, j))
            rhs[I + point_number] = -delta_t * (aij * diff_y_xixi(i, j) - 2.0 * bij * diff_y_xieta(i, j) + cij * diff_y_etaeta(i, j))

        rhs /= theta
        delta_xy = bicg(matrix, rhs)[0]

        grid_x[1:xi_max-1, 1:eta_max-1] += delta_xy[:point_number].reshape(xi_max-2, -1)
        grid_y[1:xi_max-1, 1:eta_max-1] += delta_xy[point_number:].reshape(xi_max-2, -1)

        for i in range(xi_max):
            plt.plot(grid_x[i, :], grid_y[i, :])
        plt.show()

        for j in range(eta_max):
            plt.plot(grid_x[:, j], grid_y[:, j])
        plt.show()

    # return grid_x, grid_y
"""
"""
def make_grid(z1, z2, z3):
    xi_max = z1.shape[0]    # xi方向の格子点数
    eta_max = z2.shape[0]   # eta方向の格子点数

    grid_x = np.zeros((xi_max, eta_max))    # x座標格納
    grid_y = np.zeros((xi_max, eta_max))    # y座標格納

    center = np.average(np.real(z1))
    eps = 0.01
    min_rad = np.max(np.abs(z1 + center)) + eps
    max_rad = np.min(np.abs(z3 + center)) - eps
    radius = np.linspace(min_rad, max_rad, xi_max)

    for i in range(1, xi_max-1):
        circle = center + radius[i-1] * np.exp(1j*np.linspace(0, 2.0*np.pi, eta_max))
        grid_x[i, :] = np.real(circle)
        grid_y[i, :] = np.imag(circle)
    z3 = z3[::-1]

    # 境界条件適用
    grid_x[:, 0] = np.real(z2)  # 左辺
    grid_x[0, :] = np.real(z1)  # 底辺
    grid_x[xi_max-1, :] = np.real(z3)    # 上辺
    grid_x[:, eta_max - 1] = np.real(z2)

    grid_y[0, :] = np.imag(z1)  # 底辺
    grid_y[:, 0] = np.imag(z2) + eps    # 左辺
    grid_y[xi_max - 1, :] = np.imag(z3)
    grid_y[:, eta_max - 1] = np.imag(z2) - eps


    for i in range(xi_max):
        plt.plot(grid_x[i, :], grid_y[i, :])
    plt.show()
    for j in range(eta_max):
        plt.plot(grid_x[:, j], grid_y[:, j])
    plt.show()

    diff_x_xi = lambda i, j: 0.5 * (grid_x[i+1, j] - grid_x[i-1, j])
    diff_y_xi = lambda i, j: 0.5 * (grid_y[i+1, j] - grid_y[i-1, j])
    diff_x_eta = lambda i, j: 0.5 * (grid_x[i, j+1] - grid_x[i, j-1])
    diff_y_eta = lambda i, j: 0.5 * (grid_y[i, j+1] - grid_y[i, j-1])

    Aij = lambda i, j: (diff_x_eta(i, j)**2 + diff_y_eta(i, j)**2)
    Bij = lambda i, j: (diff_x_xi(i, j) * diff_x_eta(i, j) + diff_y_xi(i, j) * diff_y_eta(i, j))
    Cij = lambda i, j: (diff_x_xi(i, j)**2 + diff_y_xi(i, j)**2)

    diff_x_xixi = lambda i, j: grid_x[i+1, j] - 2.0 * grid_x[i, j] + grid_x[i-1, j]
    diff_y_xixi = lambda i, j: grid_y[i+1, j] - 2.0 * grid_y[i, j] + grid_y[i-1, j]

    diff_x_etaeta = lambda i, j: grid_x[i, j+1] - 2.0 * grid_x[i, j] + grid_x[i, j-1]
    diff_y_etaeta = lambda i, j: grid_y[i, j+1] - 2.0 * grid_y[i, j] + grid_y[i, j-1]

    diff_x_xieta = lambda i, j: grid_x[i+1, j+1] - grid_x[i+1, j-1] - grid_x[i-1, j+1] + grid_x[i-1, j-1]
    diff_y_xieta = lambda i, j: grid_y[i+1, j+1] - grid_y[i+1, j-1] - grid_y[i-1, j+1] + grid_y[i-1, j-1]

    delta_t = 0.01
    theta = 0.5
    delta_t *= theta
    point_number = (xi_max - 2) * (eta_max - 2)
    Ni = xi_max - 2
    Nj = eta_max - 2

    # [I,I]からdeltaだけずれた位置の成分が何列目のブロック行列に属するか返す
    block_id = lambda delta: floor((I + delta) / Nj)

    for iter in range(1000):
        # matrix & rhsの準備
        rhs = np.zeros(2 * point_number)
        matrix = lil_matrix((2 * point_number, 2 * point_number))
        for I in range(point_number):
            block = block_id(0)
            i = I - block * Nj + 1    # aij対策で+1  # matrixの0番はgrid_xyの1番に相当
            j = block + 1 # aij対策で+1

            aij = Aij(i, j)
            bij = Bij(i,j)
            cij = Cij(i, j)

            matrix[I, I] = 1.0 - delta_t * 2.0 * (aij + cij)
            matrix[I + point_number, I + point_number] = 1.0 - delta_t * 2.0 * (aij + cij)
            # x-direction
            if (I > 0) and (block_id(-1) == block): # 対角成分と同一ブロック内にないときは0
                matrix[I, I-1] = delta_t * aij
                matrix[I + point_number, I + point_number - 1] = delta_t * aij
            if (I < point_number-1) and (block_id(1) == block):
                matrix[I, I+1] = delta_t * aij
                matrix[I + point_number, I + point_number + 1] = delta_t * aij

            if I > Nj - 1:
                matrix[I, I-Nj] = delta_t * cij
                matrix[I + point_number, I + point_number - Nj] = delta_t * cij
            if I < point_number - Nj:
                matrix[I, I + Nj] = delta_t * cij
                matrix[I + point_number, I + point_number + Nj] = delta_t * cij

            if ((I > (Nj - 1)) and (block_id(-(Nj - 1)) == block-1)):
                matrix[I, I-(Nj-1)] = delta_t * 2.0 * bij
                matrix[I + point_number, I + point_number - (Nj - 1)] = delta_t * 2.0 * bij
            if ((I < point_number - (Nj - 1)) and (block_id(Nj - 1) == block + 1)):
                matrix[I, I+(Nj-1)] = delta_t * 2.0 * bij
                matrix[I + point_number, I + point_number + (Nj - 1)] = delta_t * 2.0 * bij


            if ((I > (Nj + 1) - 1 and (block_id(-(Nj+1)) == block - 1))):
                matrix[I, I-(Nj+1)] = - delta_t * 2.0 * bij
                matrix[I + point_number, I + point_number - (Nj + 1)] = -delta_t * 2.0 * bij
            if ((I < point_number - (Nj + 1) - 1) and (block_id(Nj + 1) == block + 1)):
                matrix[I, I+(Nj+1)] = -delta_t * 2.0 * bij
                matrix[I + point_number, I + point_number + (Nj + 1)] = - delta_t * 2.0 * bij
            if I==10:
                print(aij)
            rhs[I] = -delta_t * (aij * diff_x_xixi(i, j) - 2.0 * bij * diff_x_xieta(i, j) + cij * diff_x_etaeta(i, j))
            rhs[I + point_number] = -delta_t * (aij * diff_y_xixi(i, j) - 2.0 * bij * diff_y_xieta(i, j) + cij * diff_y_etaeta(i, j))

        rhs /= theta
        delta_xy = bicg(matrix, rhs)[0]

        grid_x[1:xi_max-1, 1:eta_max-1] += delta_xy[:point_number].reshape(xi_max-2, -1)
        grid_y[1:xi_max-1, 1:eta_max-1] += delta_xy[point_number:].reshape(xi_max-2, -1)

        for i in range(xi_max):
            plt.plot(grid_x[i, :], grid_y[i, :])
        plt.show()

        for j in range(eta_max):
            plt.plot(grid_x[:, j], grid_y[:, j])
        plt.show()
"""
"""
def make_grid(z1, z2, z3):
    xi_max = z1.shape[0]    # xi方向の格子点数
    eta_max = z2.shape[0]   # eta方向の格子点数

    grid_x = np.zeros((xi_max, eta_max))    # x座標格納
    grid_y = np.zeros((xi_max, eta_max))    # y座標格納
    # 境界条件適用
    grid_x[:, 0] = np.real(z1)  # 底辺
    grid_x[0, :] = np.real(z2)  # 左辺
    grid_x[:, eta_max - 1] = np.real(z3)    # 上辺
    grid_x[xi_max - 1, :] = np.real(z2)

    grid_y[:, 0] = np.imag(z1)
    grid_y[0, :] = np.imag(z2)
    grid_y[:, eta_max - 1] = np.imag(z3)
    grid_y[xi_max - 1, :] = np.imag(z2)

    diff_x_xi = lambda i, j: 0.5 * (grid_x[i+1, j] - grid_x[i-1, j])
    diff_y_xi = lambda i, j: 0.5 * (grid_y[i+1, j] - grid_y[i-1, j])
    diff_x_eta = lambda i, j: 0.5 * (grid_x[i, j+1] - grid_x[i, j-1])
    diff_y_eta = lambda i, j: 0.5 * (grid_y[i, j+1] - grid_y[i, j-1])

    Aij = lambda i, j: (diff_x_eta(i, j)**2 + diff_y_eta(i, j)**2)
    Bij = lambda i, j: (diff_x_xi(i, j) * diff_x_eta(i, j) + diff_y_xi(i, j) * diff_y_eta(i, j))
    Cij = lambda i, j: (diff_x_xi(i, j)**2 + diff_y_xi(i, j)**2)

    point_number = (xi_max - 2) * (eta_max - 2)
    Ni = xi_max - 2
    Nj = eta_max - 2
    matrix = lil_matrix((2*point_number, 2*point_number))

    for iter in range(100):
        # matrix & rhsの準備
        rhs = np.zeros(2 * point_number)
        for J in range(point_number):
            row = floor(J / Nj)
            i = J - row * Nj    # aij対策で+1  # matrixの0番はgrid_xyの1番に相当
            j = row + 1 # aij対策で+1

            aij = Aij(i, j)
            bij = Bij(i,j)
            cij = Cij(i, j)

            matrix[J, J] = -2.0 * (aij + cij)
            matrix[J + point_number, J + point_number] = -2.0 * (aij + cij)
            # x-direction
            if J > 0:
                matrix[J, J-1] = aij
                matrix[J + point_number, J + point_number - 1] = aij
            if J < point_number-1:
                matrix[J, J+1] = aij
                matrix[J + point_number, J + point_number + 1] = aij

            if J > Nj:
                matrix[J, J-Nj] = cij
                matrix[J + point_number, J + point_number - Nj] = cij
            if J < point_number - 1 - Nj:
                matrix[J, J+Nj] = cij
                matrix[J + point_number, J + point_number + Nj] = cij

            if J > (Nj - 1):
                matrix[J, J-(Nj-1)] = 2.0 * bij
                matrix[J + point_number, J + point_number - (Nj - 1)] = 2.0 * bij
            if J < point_number - 1 - (Nj - 1):
                matrix[J, J+(Nj-1)] = 2.0 * bij
                matrix[J + point_number, J + point_number + (Nj - 1)] = 2.0 * bij

            if J > (Nj + 1):
                matrix[J, J-(Nj+1)] = -2.0 * bij
                matrix[J + point_number, J + point_number - (Nj + 1)] = -2.0 * bij
            if J < point_number - 1 - (Nj + 1) - 1:
                matrix[J, J+(Nj+1)] = -2.0 * bij
                matrix[J + point_number, J + point_number + (Nj + 1)] = -2.0 * bij

            boundary = 0
            if i == 1:
                boundary += 1
            if i == Ni:
                boundary += 2
            if j == 1:
                boundary += 3
            if j == Nj:
                boundary += 7

            if boundary == 1:   # x0境界
                rhs[J] = 2.0 * bij * grid_x[0, j-1] - aij * grid_x[0, j] - 2.0 * bij * grid_x[0, j + 1]
                rhs[J + point_number] = 2.0 * bij * grid_y[0, j - 1] - aij * grid_y[0, j] - 2.0 * bij * grid_y[0, j + 1]

            if boundary == 3:   # y0境界
                rhs[J] = 2.0 * bij * (grid_x[i-1, 0] - grid_x[i+1, 0]) - cij * grid_x[i, 0]
                rhs[J + point_number] = 2.0 * bij * (grid_y[i - 1, 0] - grid_y[i + 1, 0]) - cij * grid_y[i, 0]

            if boundary == 2:   # xN境界
                rhs[J] = -2.0 * bij * grid_x[xi_max - 1, j - 1] - aij * grid_x[xi_max - 1, j] + 2.0 * bij * grid_x[xi_max-1, j + 1]
                rhs[J + point_number] = -2.0 * bij * grid_y[xi_max - 1, j - 1] - aij * grid_y[xi_max - 1, j] + 2.0 * bij * grid_y[xi_max - 1, j + 1]

            if boundary == 7:   # yN境界
                rhs[J] = -2.0 * bij * (grid_x[i - 1, eta_max-1] - grid_x[i + 1, eta_max-1]) - cij * grid_x[i, eta_max-1]
                rhs[J + point_number] = -2.0 * bij * (grid_y[i - 1, eta_max-1] - grid_y[i + 1, eta_max-1]) - cij * grid_y[i, eta_max-1]

            if boundary == 4:   # x0 and y0
                rhs[J] = 2.0 * bij * (grid_x[0, 0] - grid_x[2, 0] - grid_x[0, 2]) - aij * grid_x[0, 1] - cij * grid_x[1, 0]
                rhs[J + point_number] = 2.0 * bij * (grid_y[0, 0] - grid_y[2, 0] - grid_y[0, 2]) - aij * grid_y[0, 1] - cij * grid_y[1, 0]

            if boundary == 8:  # x0 and yN
                rhs[J] = 2.0 * bij * (grid_x[0, eta_max-3] - grid_x[0, eta_max-1] + grid_x[2, eta_max-1]) - aij * grid_x[0, eta_max-2] -cij * grid_x[1, eta_max-1]
                rhs[J + point_number] = 2.0 * bij * (grid_y[0, eta_max - 3] - grid_y[0, eta_max - 1] + grid_y[2, eta_max - 1]) - aij * grid_y[0, eta_max - 2] - cij * grid_y[1, eta_max - 1]

            if boundary == 5:   # y0 and xN
                rhs[J] = 2.0 * bij * (grid_x[xi_max-3, 0] - grid_x[xi_max-1, 0] + grid_x[xi_max-1, 2]) - aij * grid_x[xi_max-1, 1] - cij * grid_x[xi_max-2, 0]
                rhs[J + point_number] = 2.0 * bij * (grid_y[xi_max - 3, 0] - grid_y[xi_max - 1, 0] + grid_y[xi_max - 1, 2]) - aij * grid_y[xi_max - 1, 1] - cij * grid_y[xi_max - 2, 0]

            if boundary == 9:   # xN and yN
                rhs[J] = 2.0 * bij * (-grid_x[xi_max-1, eta_max-2] - grid_x[xi_max-3, eta_max-1] + grid_x[xi_max-1, eta_max-1]) - aij * grid_x[xi_max-1, eta_max-2] - cij * grid_x[xi_max-2, eta_max-1]
                rhs[J + point_number] = 2.0 * bij * (-grid_y[xi_max - 1, eta_max - 2] - grid_y[xi_max - 3, eta_max - 1] + grid_y[xi_max - 1, eta_max - 1]) - aij * grid_y[xi_max - 1, eta_max - 2] - cij * grid_y[xi_max - 2, eta_max - 1]

        xy = bicg(matrix, rhs)[0]

        grid_x[1:xi_max-1, 1:eta_max-1] = xy[:point_number].reshape(xi_max-2, -1).T
        grid_y[1:xi_max-1, 1:eta_max-1] = xy[point_number:].reshape(xi_max-2, -1).T

        plt.plot(grid_x.reshape(-1), grid_y.reshape(-1), "x")
        plt.show()


    # return grid_x, grid_y
"""

"""
def make_grid(z1, z2, z3):
    xi_max = z1.shape[0]    # xi方向の格子点数
    eta_max = z2.shape[0]   # eta方向の格子点数

    grid_x = np.zeros((xi_max, eta_max))    # x座標格納
    grid_y = np.zeros((xi_max, eta_max))    # y座標格納
    # 境界条件適用
    grid_x[:, 0] = np.real(z1)  # 底辺
    grid_x[0, :] = np.real(z2)  # 左辺
    grid_x[:, eta_max - 1] = np.real(z3)    # 上辺
    grid_x[xi_max - 1, :] = np.real(z2)

    grid_y[:, 0] = np.imag(z1)
    grid_y[0, :] = np.imag(z2)
    grid_y[:, eta_max - 1] = np.imag(z3)
    grid_y[xi_max - 1, :] = np.imag(z2)
    # 計算準備
    rhs_x = np.zeros((xi_max-2, eta_max-2))   # δxに関する方程式の右辺格納
    rhs_y = np.zeros((xi_max-2, eta_max-2))   # δyに関する方程式の右辺格納
    mid_x = np.zeros((xi_max-2, eta_max-2))   # δxに関する方程式の中間解格納
    mid_y = np.zeros((xi_max-2, eta_max-2))   # δyに関する方程式の中間解格納
    del_x = np.zeros((xi_max - 2, eta_max - 2))  # δxに関する方程式の解格納
    del_y = np.zeros((xi_max - 2, eta_max - 2))  # δyに関する方程式の解格納

    mXIdiag = np.zeros(eta_max-2)  # xi方向に対応する3重対角行列の対角成分格納用
    mXIside = np.zeros(eta_max-2) # xi方向に対応する3重対角行列の左右成分格納用

    mETAdiag = np.zeros(xi_max-2) # eta方向に対応する3重対角行列の対角成分格納用
    mETAside = np.zeros(xi_max-2)  # eta方向に対応する3重対角行列の左右成分格納用

    diff_x_xi = lambda i, j: 0.5 * (grid_x[i+1, j] - grid_x[i-1, j])
    diff_y_xi = lambda i, j: 0.5 * (grid_y[i+1, j] - grid_y[i-1, j])
    diff_x_eta = lambda i, j: 0.5 * (grid_x[i, j+1] - grid_x[i, j-1])
    diff_y_eta = lambda i, j: 0.5 * (grid_y[i, j+1] - grid_y[i, j-1])

    Aij = lambda i, j: (diff_x_eta(i, j)**2 + diff_y_eta(i, j)**2)
    Bij = lambda i, j: (diff_x_xi(i, j) * diff_x_eta(i, j) + diff_y_xi(i, j) * diff_y_eta(i, j))
    Cij = lambda i, j: (diff_x_xi(i, j)**2 + diff_y_xi(i, j)**2)

    diff_x_xixi = lambda i, j: grid_x[i+1, j] - 2.0 * grid_x[i, j] + grid_x[i-1, j]
    diff_y_xixi = lambda i, j: grid_y[i+1, j] - 2.0 * grid_y[i, j] + grid_y[i-1, j]

    diff_x_etaeta = lambda i, j: grid_x[i, j+1] - 2.0 * grid_x[i, j] + grid_x[i, j-1]
    diff_y_etaeta = lambda i, j: grid_y[i, j+1] - 2.0 * grid_y[i, j] + grid_y[i, j-1]

    diff_x_xieta = lambda i, j: grid_x[i+1, j+1] - grid_x[i+1, j-1] - grid_x[i-1, j+1] + grid_x[i-1, j-1]
    diff_y_xieta = lambda i, j: grid_y[i+1, j+1] - grid_y[i+1, j-1] - grid_y[i-1, j+1] + grid_y[i-1, j-1]

    delta_t = 0.1

    for iter in range(100):
        #rhs計算
        for i in range(1, xi_max-1):
            for j in range(1, eta_max-1):
                aij = Aij(i, j)
                bij = Bij(i,j)
                cij = Cij(i, j)

                rhs_x[i-1, j-1] = - delta_t * (aij * diff_x_xixi(i, j)
                                               - 2.0 * bij * diff_x_xieta(i, j)
                                               + cij * diff_x_etaeta(i, j))

                rhs_y[i-1, j-1] = - delta_t * (aij * diff_y_xixi(i, j)
                                               - 2.0 * bij * diff_y_xieta(i, j)
                                               + cij * diff_y_etaeta(i, j))

        # ADI 1st step
        for i in range(1, xi_max-1):
            for j in range(1, eta_max-1):
                aij = Aij(i, j)
                mXIdiag[j-1] = 1.0 - 2.0*delta_t*aij
                mXIside[j-1] = delta_t*aij

            mid_x[i-1, :] = tridiagonal_matrix_algorithm(mXIside, mXIdiag, mXIside, rhs_x[i-1, :], eta_max-2)
            # print(mid_x[i-1, :])
            matrix = np.zeros((eta_max-2, eta_max-2))
            matrix[0, 0] = 1.0 - 2.0*delta_t*Aij(1, 1)
            matrix[0, 1] = delta_t*Aij(1, 1)
            for k in range(1,eta_max-3):
                aij = Aij(k+1, k+1)
                matrix[k, k-1] = delta_t*aij
                matrix[k, k] = 1.0 - 2.0*delta_t*aij
                matrix[k, k+1] = delta_t*aij
            matrix[eta_max-3, eta_max-3] = 1.0 - 2.0*delta_t*Aij(1, eta_max-2)
            matrix[eta_max - 3, eta_max - 4] = delta_t * Aij(1, eta_max-2)
            mid_y[i-1, :] = tridiagonal_matrix_algorithm(mXIside, mXIdiag, mXIside, rhs_y[i-1, :], eta_max-2)

        # ADI 2nd step
        for j in range(1, eta_max-1):
            for i in range(1, xi_max-1):
                cij = Cij(i, j)
                mETAdiag[i-1] = 1.0 - 2.0 * delta_t * cij
                mETAside[i-1] = delta_t * cij

            del_x[:, j-1] = tridiagonal_matrix_algorithm(mETAside, mETAdiag, mETAside, mid_x[:, j-1], xi_max-2)
            del_y[:, j-1] = tridiagonal_matrix_algorithm(mETAside, mETAdiag, mETAside, mid_y[:, j-1], xi_max-2)


        # delta-form
        grid_x[1:xi_max-1, 1:eta_max-1] += del_x
        grid_y[1:xi_max-1, 1:eta_max-1] += del_y

        res = (np.sum(np.abs(del_x)) + np.sum(np.abs(del_y)))
        print(np.argmax(np.abs(del_x)), np.argmax(np.argmin(del_y)))

        print(res)
        plt.plot(grid_x.reshape(-1), grid_y.reshape(-1), "x")
        plt.show()
        if res < 0.001:
            break
"""


"""
# 複素数列z1, z2,...を, 2倍の長さの実数列x1, y1, x2, y2, ...に変換
def line_up_z2xy(z, two_rows=False):
    
x = np.real(z).reshape(1, -1)
y = np.imag(z).reshape(1, -1)
xy = np.vstack((x, y)).T.reshape(-1)
return xy

if two_rows:
    return np.vstack((np.real(z).reshape(1, -1), np.imag(z).reshape(1, -1))).T
else:
    return np.vstack((np.real(z).reshape(1, -1), np.imag(z).reshape(1, -1))).T.reshape(-1)


# 物体表面座標点列を内側に縮小し，0番目の仮データを用意する
def offset_surface(z):
size = z.shape[0] - 1
delta = z[1:] - z[:size]
theta = np.zeros(size)
theta[0] = 0.5 * np.angle(-delta[size-1]/delta[0]) + np.angle(delta[0])
theta[1:] = 0.5 * np.angle(-delta[:size-1]/delta[1:]) + np.angle(delta[1:])
offset_quantity = np.average(np.abs(delta))
z0 = np.zeros_like(z)
z0[:size] = z[:size] + offset_quantity * np.exp(1j * theta)
z0[size] = z0[0]
return z0

# 係数行列を計算
def made_coefficient_matrix(x0, x1, y0, y1, phi0, phi1, V0, V1):
def make_invA(i):
    detA = 4.0 / ((x1[i] - x0[i])**2 + (y1[i] - y0[i])**2)
    invA[0, 0] = detA * (x1[i] - x0[i])
    invA[0, 1] = detA * (y1[i] - y0[i])
    invA[1, 0] = detA * (y1[i] - y0[i])
    invA[1, 1] = - detA * (x1[i] - x0[i])
    return invA

def make_r(i):
    r[0] = -phi1[i]*V1[i]**2 - phi0[i]*V0[i]**2
    r[1] = V1[i] + V0[i]
    return r

size = x0.shape[0]
coefMat = np.zeros((2*size, 2*size))
const_vector = np.zeros(2*size)
invA = np.zeros((2, 2))
B = np.zeros((2, 2))
r = np.zeros(2)

# i = 0
invA = make_invA(0)
B[0, 0] = 0.5 * (x1[1] - x1[size-1])
B[0, 1] = 0.5 * (y1[1] - y1[size-1])
B[1, 0] = -0.5 * (y1[1] - y1[size-1])
B[1, 1] = 0.5 * (x1[1] - x1[size-1])
B *= 2
coefMat[0:2, 2*size - 2:2*size] = -1.0
coefMat[0:2, 0:2] = np.dot(invA, B)
coefMat[0:2, 2:4] = 1.0

r = make_r(0)
stabilizer = np.array([x1[1] - 2.0 * x1[0] + x1[size - 1], y1[1] - 2.0 * y1[0] + y1[size - 1]])
const_vector[0:2] = np.dot(invA, r + np.dot(B, np.array([x1[0], y1[0]]))) + stabilizer

for i in range(1, size-1):
    k = 2 * i # - 1
    invA = make_invA(i)
    B[0, 0] = 0.5 * (x1[i+1] - x1[i-1])
    B[0, 1] = 0.5 * (y1[i+1] - y1[i-1])
    B[1, 0] = -0.5 * (y1[i+1] - y1[i-1])
    B[1, 1] = 0.5 * (x1[i+1] - x1[i-1])
    B *= 2
    coefMat[k:k+2, k-2:k] = -1.0
    coefMat[k:k+2, k:k+2] = np.dot(invA, B)
    coefMat[k:k+2, k+2:k+4] = 1.0
    r = make_r(i)
    stabilizer = np.array([x1[i+1] - 2.0 * x1[i] + x1[i-1], y1[i+1] - 2.0 * y1[i] + y1[i-1]])
    const_vector[k:k+2] = np.dot(invA, r + np.dot(B, np.array([x1[i], y1[i]]))) + stabilizer

# i = size - 1
invA = make_invA(size-1)
B[0, 0] = 0.5 * (x1[0] - x1[size-2])
B[0, 1] = 0.5 * (y1[0] - y1[size-2])
B[1, 0] = -0.5 * (y1[0] - y1[size-2])
B[1, 1] = 0.5 * (x1[0] - x1[size-2])
B *= 2
coefMat[2*size-2:2*size, 2*size-4:2*size-2] = -1.0
coefMat[2*size-2:2*size, 2*size-2:2*size] = np.dot(invA, B)
coefMat[2*size-2:2*size, 0:2] = 1.0
r = make_r(size-1)
stabilizer = np.array([x1[0] - 2.0 * x1[size-1] + x1[size-2], y1[0] - 2.0 * y1[size-1] + y1[size-2]])
const_vector[2*size-2:2*size] = np.dot(invA, r + np.dot(B, np.array([x1[size-1], y1[size-1]]))) + stabilizer

return coefMat, const_vector

# 計算済みの閉曲線から一回り大きな閉曲線を求める
def get_next_closed_curve(z0, z1):
a = 1
"""
