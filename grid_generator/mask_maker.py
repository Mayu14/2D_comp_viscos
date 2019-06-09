# -- coding: utf-8 --
import cv2
import numpy as np
import matplotlib.pyplot as plt
import sys

class VPmask(object):
    def __init__(self, x, y, canvas_size, start_point="lb", imglim=None, deform="None", obj_margin=4,
                 default_canvas=[1024,1024], outerframe_margin=3, line_width=0.5):
        """
        物体表面のx,yの座標の配列を入力し，rgb配列に変換する(x,yは一筆書きできるように整列済みであるとする)
        :param x:   (list as [ndarray x-coords obj1, ndarray x-coords obj2, ... ])    # object surface points coordinates of x
        :param y:   (list as [ndarray y-coords obj1, ndarray y-coords obj2, ... ])    # object surface points coordinates of y
        :param canvas_size: (list as [x_resolution, y_resolution])  # image canvas size (pixel)
        :param start_point:   (character) # A start point of floodFill (start point must be outside of the object) "lb" > "leftbottom", "rt" > "righttop", "lt", "rb"
        :param imglim:    (list as [[x_min, x_max], [y_min, y_max]])    # upper and lower limit of drawing range (x & y-direction)  # 指定がないときは，長軸方向に合わせた正方形領域を選ぶ
        :param deform: (character or float) "None", "Fit", or x/y (aspect_ratio)
        :param obj_margin: (int) set margin from object to boundary (margin = obj_margin * number of pixels)
        :param outerframe_margin: (int) outer frame size of matplotlib's graph
        :param line_width: (float) thickness of curve representing object   (大きいと形状表現精度低下・小さすぎるとマスク生成失敗率上昇)
        :return:
        """
        self.x = x
        self.y = y
        self.object_number = len(x)
        self.__get_min_and_max__()


        self.canvas_size = default_canvas
        self.__set_canvas__(canvas_size)

        self.obj_margin = obj_margin
        self.outerframe_margin = outerframe_margin

        # 閉曲線を変形する
        self.deform_curve(deform)

        if type(imglim) == type(None):
            imglim = self.__get_imglim__()

        self.xlim = imglim[0]
        self.ylim = imglim[1]

        self.__set_start_points__(start_point)

        self.lw = line_width
        # 閉曲線x, yが一筆書きになるよう並べ替え処理を入れたい

        # 閉曲線から画像を生成する
        self.img = self.get_curve_img()
        self.mask = self.get_mask(self.img)
        self.resize()

    def __get_min_and_max__(self):
        min_x = np.inf  # 初期化
        min_y = np.inf
        max_x = -np.inf
        max_y = -np.inf
        for i in range(self.object_number): # すべての物体について座標の最大値・最小値を求める
            min_x = min(min_x, np.min(self.x[i]))
            min_y = min(min_y, np.min(self.y[i]))
            max_x = max(max_x, np.max(self.x[i]))
            max_y = max(max_y, np.max(self.y[i]))

        self.min = [min_x, min_y]
        self.max = [max_x, max_y]

    def __set_canvas__(self, canvas_size):
        for i in range(2):
            if self.canvas_size[i] < canvas_size[i]:
                self.canvas_size[i] = int(1.2 * canvas_size[i])
        self.final_canvas = canvas_size


    def deform_curve(self, deform):
        self.deform = False
        if deform != "None":
            self.deform = True
            if deform == "Fit":
                deform = (self.max[0] - self.min[0]) / (self.max[1] - self.min[1])  # x/y
            center_y = 0.5 * (self.max[1] + self.min[1])

            for i in range(self.object_number):
                self.y[i] = deform * (self.y[i] - center_y) + center_y
            self.aspect_ratio = deform
            self.__get_min_and_max__()

    def __get_imglim__(self):
        min_x = self.min[0]
        min_y = self.min[1]
        max_x = self.max[0]
        max_y = self.max[1]
        dx = max_x - min_x    # x方向の幅
        dy = max_y - min_y    # y方向の幅

        if dx > dy:
            eps = dx / float(self.canvas_size[0])
        else:
            eps = dy / float(self.canvas_size[1])

        eps *= float(self.obj_margin)   # marginで指定されたpixel数だけmarginを確保する
        if self.deform:
            return [[min_x - eps, max_x + eps], [min_y - eps, max_y + eps]]  # 長軸方向に余白を与える
        else:
            lim = [(min(min_x, min_y) - eps), (max(max_x, max_y) + eps)]  # 長軸方向に余白を与える
            return [lim, lim]    # 長軸方向基準の正方形領域とする

    def __set_internal_points__(self, internal_point):
        int_pts = []
        if type(internal_point) == type(None):
            for i in range(self.object_number):
                mid_x = (np.min(self.x[i]) + np.max(self.x[i])) / 2.0   # 物体の中心座標を求めて
                mid_y = (np.min(self.y[i]) + np.max(self.y[i])) / 2.0
                pix_xy = self.coords2pixel([mid_x, mid_y])  # 対応する画素に変換
                int_pts.append(pix_xy)

        elif internal_point[0][0].dtype == np.zeros(1, dtype=float).dtype:
            for i in range(self.object_number):
                pts_x = internal_point[i][0]    # 登録されている座標を読み込んで
                pts_y = internal_point[i][1]
                pix_xy = self.coords2pixel([pts_x, pts_y])  # 対応する画素に変換
                int_pts.append(pix_xy)

        return int_pts

    def __set_start_points__(self, start_points):
        def left_or_bottom():
            return 5

        def right_or_top(right=True):
            if right:
                direction = 0
            else:
                direction = 1
            return self.canvas_size[direction] - 5

        if start_points == "lb":
            xs = left_or_bottom()
            ys = left_or_bottom()
        elif start_points == "lt":
            xs = left_or_bottom()
            ys = right_or_top(right=False)
        elif start_points == "rb":
            xs = right_or_top()
            ys = left_or_bottom()
        else:
            xs = right_or_top()
            ys = right_or_top(right=False)

        self.start_points = [xs, ys]


    def coords2pixel(self, pts_xy):
        """
        座標を入力すると，それがおよそどのあたりの画素に対応しているかを返す
        :param pts_xy: (list as [x, y]) # 座標
        :return: 画素[i, j]
        """
        percent_x = (pts_xy[0] - self.min[0]) / (self.max[0] - self.min[0])  # その座標が空間全体において何%の位置にいるか求めて
        percent_y = (pts_xy[1] - self.min[1]) / (self.max[1] - self.min[1])
        pix_x = int(self.canvas_size[0] * percent_x)  # 画素数に位置の%を掛けて，対応する画素を登録
        pix_y = int(self.canvas_size[1] * percent_y)
        return [pix_x, pix_y]

    def get_curve_img(self):
        out_margin = 2*self.outerframe_margin
        # *** プロット処理 ***
        fig = plt.figure(figsize=(self.canvas_size[0] + out_margin, self.canvas_size[1] + out_margin), dpi=1)  # figureオブジェクトの初期化
        fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0, wspace=0.0, hspace=0.0)
        ax = fig.add_subplot(1,1,1) # subplotの追加
        for i in range(self.object_number):
            ax.plot(self.x[i], self.y[i], lw=self.lw)   # x,yをプロットする
        ax.set_xlim(self.xlim[0], self.xlim[1])
        ax.set_ylim(self.ylim[0], self.ylim[1])
        ax.tick_params(axis="both", which="both",
                       bottom=False, top=False, left=False, right=False,
                       labelbottom=False, labeltop=False, labelleft=False, labelright=False)

        # matplotlibの出力結果をnumpy配列に変換
        fig.canvas.draw()
        img = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
        img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        plt.close(fig)
        return img

    def get_mask(self, img):
        # *** エッジ検出画像をマスク画像として用意 ***
        #  mask画像の画素値が0の場所がfloodFillの処理対象となる
        # 従ってエッジ(物体境界)検出画像をmaskにすることでエッジの内部, 外部を自動で分離することができる
        mask = cv2.Canny(img, 100, 200) # 1st:input_img, 2nd:hysteresis min criteria 3rd:hysteresis max criteria
        # *** マスク画像の拡張。（画像の上下左右に1ラインずつ拡張）
        mask = cv2.copyMakeBorder(mask, 1, 1, 1, 1, cv2.BORDER_REPLICATE)   #cv2.BORDER_REPLICATEにより，端部は手前のピクセルと同じ値を取る

        # *** 領域塗りつぶしと結果の表示
        # cv2.floodFill( img, mask, (120,120), (0,0,0))   # 1st:img, 2nd:mask 3rd:seed point 4th:new value

        # cv2.floodFill(self.img, mask, (int(sz[0]/2), int(sz[1]/2)), (0,0,0))
        # cv2.floodFill(self.img, mask, (self.internal_point[0][0], self.internal_point[0][1]), (0,0,0))
        cv2.floodFill(img, mask, (self.start_points[0], self.start_points[1]), (0, 0, 0))    # 物体じゃない側を塗りつぶす
        # cv2.imshow("canny",mask)
        # self.img = 255 - self.img   # reverse
        xy_min = self.outerframe_margin
        x_max = (img.shape[0]) - xy_min
        y_max = (img.shape[1]) - xy_min
        return img[xy_min:x_max, xy_min:y_max, 0]

    def resize(self, interpolate = cv2.INTER_LANCZOS4):
        width = self.final_canvas[0]
        height = self.final_canvas[1]
        self.mask = cv2.resize(self.mask, (width, height), interpolation=interpolate)
        self.img = cv2.resize(self.img, (width, height), interpolation=interpolate)

    def show_mask(self):
        cv2.imshow("mask", self.mask)
        cv2.waitKey()

    def show_img(self):
        cv2.imshow("img", self.img)
        cv2.waitKey()

    def save_img(self, fname_without_extension, format="png"):
        if self.deform:
            fname_without_extension += "__" + str(self.aspect_ratio) + "__"
        fname = fname_without_extension + "." + format
        cv2.imwrite(fname, self.img)

def test_independent():
    x=[]
    y=[]
    for itht in range(361):
        tht = itht*np.pi/180.0
        num_polygon = 7
        r = 2+1*np.sin(tht*num_polygon)
        x.append(r*np.cos(tht))
        y.append(r*np.sin(tht))


def test_naca():
    from grid_generator.naca_4digit_test import Naca_4_digit
    attack_angle_deg = 0
    grid_resolution = 200
    naca4 = "12099"
    naca = Naca_4_digit(naca4, attack_angle_deg, grid_resolution, quasi_equidistant=True)
    naca.plot()
    x = np.concatenate([naca.x_u, naca.x_l[::-1]])
    y = np.concatenate([naca.y_u, naca.y_l[::-1]])

    mid_res = int(grid_resolution/2)

    int_pts = [naca.equidistant_x[mid_res], 0.5 * (naca.equidistant_y_u[mid_res] + naca.equidistant_y_l[mid_res])]

    return [x],[y]

def demo1_xy(radius = 0.2, center=[-0.4,-0.2], aspect_ratio=0.7):
    theta = np.linspace(start=0, stop=2.0*np.pi, num=360)
    x = radius * np.cos(theta) + center[0]
    y = radius * aspect_ratio * np.sin(theta) + center[1]
    return [x],[y]

def demo2_xy(center=0.3 + 1j * 0.5, outer_radius = 0.4):
    num_vertex = 5
    radius_ratio = 0.35
    inner_radius = radius_ratio * outer_radius
    outer_vertex = outer_radius * np.exp(1j * np.linspace(start=0, stop=2.0 * np.pi, num=num_vertex + 1)).reshape(1,-1) + center
    phase_differnce = 2.0 * np.pi / (2.0 * num_vertex)

    inner_vertex = inner_radius * np.exp(1j * np.linspace(start=0 - phase_differnce, stop=2.0 * np.pi - phase_differnce,
                         num=num_vertex + 1)).reshape(1, -1) + center
    vertex = np.concatenate((inner_vertex, outer_vertex)).T.reshape(-1)
    return [np.real(vertex)], [np.imag(vertex)]


if __name__ == '__main__':
    args = sys.argv

    if len(args) == 1:
        x,y = test_naca()
    elif args[1] == 0:
        x = []
        y = []
        for i in range(2, len(args)):
            fname = args[i]
            x.append(np.loadtxt(fname + "x.csv", deliminater=","))
            y.append(np.loadtxt(fname + "y.csv", deliminater=","))

    elif args[1] == 1:
        print("demo mode 1")
        x, y = demo1_xy()
    elif args[1] == 2:
        print("demo mode 2")
        x, y = demo2_xy()
    elif args[1] == 3:
        print("demo mode 2")
        x_1, y_1 = demo1_xy()
        x_2, y_2 = demo2_xy()
        x = [x_1[0], x_2[0]]
        y = [y_1[0], y_2[0]]
    else:
        exit()

    size = 512
    canvas = [size,size]
    #imglim = [[-1.1,1.1],[-1.1,1.1]]
    mask = VPmask(x, y, canvas, deform="Fit", obj_margin=4, start_point="rt")
    mask.show_img()
    # main()

