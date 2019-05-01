# -- coding: utf-8 --
# 入力されたデータ列a,bから散布図を作成する
import numpy as np
import matplotlib.pyplot as plt

def make_scatter_plot(data_a, data_b, label_a, label_b, resolution=200, path="", fname="scat", log_scale=True, dash_line=True, fix_scale=False):
    size = data_a.shape[0]
    size_b = data_b.shape[0]
    if size != size_b:
        return "error:invalid data size!"

    max_data = max(np.max(data_a), np.max(data_b))
    min_data = min(np.min(data_a), np.min(data_b))

    # データの調整
    data_a = data_a - min_data
    data_b = data_b - min_data
    # データラベル
    label = np.linspace(start=min_data, stop=max_data, num=resolution + 1)
    # データプロット用の増分
    delta = (max_data - min_data) / resolution
    # 散布図作図用配列
    scatter = np.zeros((resolution + 1, resolution + 1))
    # 1次元境界箱による接触判定
    boundary_box_1d = lambda data, delta: np.rint(data / delta)

    data_scat_a = boundary_box_1d(data_a, delta)
    data_scat_b = boundary_box_1d(data_b, delta)

    for i in range(size):
        scatter[int(data_scat_a[i]), int(data_scat_b[i])] += 1
    
    if log_scale:
        scatter += 1
        scatter = np.log10(scatter)

    imshow = False
    plt.figure()
    if imshow:
        plt.imshow(scatter)
        plt.colorbar()
    else:
        plt.pcolormesh(label, label, scatter, cmap='jet')
        plt.clim(0,2.0)
        pp = plt.colorbar(orientation="vertical")
        pp.set_label("Number of times (log)", fontname = "Arial", fontsize = 12)
    plt.xlabel(label_a)
    plt.ylabel(label_b)
    if fix_scale:
        plt.xlim(-4, 4)
        plt.ylim(-4, 4)
        fname += "fixed"
    plt.savefig(path + fname + "_original.png")
    plt.close()
    
    plt.figure()
    if imshow:
        plt.imshow(scatter)
        plt.colorbar()
    else:
        plt.pcolormesh(label, label, scatter, cmap='jet')
        plt.clim(0, 2.0)
        pp = plt.colorbar(orientation = "vertical")
        pp.set_label("Number of times (log)", fontname = "Arial", fontsize = 12)
    plt.xlabel(label_a)
    plt.ylabel(label_b)
    if fix_scale:
        plt.xlim(-4, 4)
        plt.ylim(-4, 4)
        fname += "fixed"
    plt.plot(label, label, color="white", linestyle="dashed")
    plt.savefig(path + fname + "_with_center_line.png")
    plt.close()
    
def test():
    def test_func(x, r):
        return x ** 2 + 3 * x + r * np.random.standard_normal(x.shape[0])

    size = 100000
    x = np.random.normal(loc=0.0, scale=1.0, size=size)
    a = test_func(x, 0.0)
    b = test_func(x, 1.0)
    # plt.plot(x, a, "x")
    # plt.plot(x, b, "x")
    # plt.show()
    return a, b


if __name__ == '__main__':
    a, b = test()
    label_a = "a"
    label_b = "b"
    make_scatter_plot(a, b, label_a, label_b, resolution=200)