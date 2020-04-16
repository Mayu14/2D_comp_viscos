# -- coding: utf-8 --
# 入力されたデータ列a,bから散布図を作成する
import numpy as np
import matplotlib.pyplot as plt


def make_scatter_plot(data_a, data_b, label_a, label_b, resolution=200, path="", fname="scat", log_scale=True,
                      dash_line="both", prioritize_a=False, max_value=None, min_value=None, title=None):
    if (dash_line != "both") and (dash_line != "yes") and (dash_line != "no"):
        raise ValueError

    size = data_a.shape[0]
    size_b = data_b.shape[0]
    if size != size_b:
        return "error:invalid data size!"

    if not max_value is None:
        max_data = max_value
    else:
        max_data = max(np.max(data_a), np.max(data_b))
        if prioritize_a:
            max_data = np.max(data_a)

    if not min_value is None:
        min_data = min_value
    else:
        min_data = min(np.min(data_a), np.min(data_b))
        if prioritize_a:
            min_data = np.min(data_a)

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
        if (not int(data_scat_b[i]) < 0) and (not int(data_scat_b[i]) > resolution - 1):
            if (not int(data_scat_a[i]) < 0) and (not int(data_scat_a[i]) > resolution - 1):
                scatter[int(data_scat_a[i]), int(data_scat_b[i])] += 1

    pp_label = "Number of times"
    if log_scale:
        scatter += 1
        scatter = np.log10(scatter)
        pp_label += " (log)"
    imshow = False
    if (dash_line == "both") or (dash_line == "no"):
        plt.figure(figsize = (10, 8))
        if not title is None:
            plt.title(title)
        if imshow:
            plt.imshow(scatter)
            plt.colorbar()
        else:
            plt.pcolormesh(label, label, scatter, cmap='jet')
            plt.clim(0, 2.0)
            pp = plt.colorbar(orientation="vertical")
            pp.set_label(pp_label)#, fontname="Arial", fontsize=12)
        plt.xlabel(label_a)
        plt.ylabel(label_b)
        footer = ""
        if (dash_line == "both") or (not "png" in fname):
            footer += "_original.png"
        plt.savefig(path + fname + footer)
        plt.close()

    if (dash_line == "both") or (dash_line == "yes"):
        plt.figure()
        if imshow:
            plt.imshow(scatter)
            plt.colorbar()
        else:
            plt.pcolormesh(label, label, scatter, cmap='jet')
            plt.clim(0, 2.0)
            pp = plt.colorbar(orientation="vertical")
            pp.set_label(pp_label, fontname="Arial", fontsize=12)
        plt.xlabel(label_a)
        plt.ylabel(label_b)

        plt.plot(label, label, color="white", linestyle="dashed")
        footer = ""
        if (dash_line == "both") or (not "png" in fname):
            fname += "_with_center_line.png"
        plt.savefig(path + fname + footer)
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