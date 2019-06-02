# coding: utf-8
import numpy as np
from other_tools.legacy_vtk_writer import MyStructuredVTKWriter
from naca_4digit_test import Naca_4_digit as N4d    # use only test
# 標準物体形状は[0, 1]×[0, 1]領域に描画されるものとする
# 即ち，広義の物体形状をmax(xmax - xmin, ymax - ymin)にて規格化したのち，原点が(0.5, 0.5)となるよう平行移動したものを標準物体形状と定める
# Shapeオブジェクトは，(上面関数，下面関数，親格子解像度，マージンセル数)によって初期化される

class Shape(object):
	def __init__(self, y_upper, y_lower, resolution, aspect_ratio=1.0, resize=True, name="Shape", auto_reshape=True):
		self.upper = y_upper    # ndarray
		self.lower = y_lower    # ndarray
		self.aspect = np.array([1.0, aspect_ratio])
		self.resolution = resolution
		self.grid_dx = 1.0 / resolution # [0, 1]区間をresolution分割    # 厳密には(size + 1)であるが，形状誤差の許容値dx/sizeを与えて簡略化することで境界における問題を無くす
		self.grid_dy = self.grid_dx * aspect_ratio
		self.resize = resize
		self.name = name
		self.auto_reshape = auto_reshape
		self.get_mask()

	def change_aspect(self):
		y_range = np.max(self.upper) - np.min(self.lower)	# yの範囲を求めて
		y_center = (np.max(self.upper) + np.min(self.lower)) / 2.0	# yの中心を求めて
		ratio = 1.0 / y_range	# y方向の拡大率を求めて
		self.upper = (self.upper - y_center) * ratio + y_center	# 中心を0に合わせてから拡大，中心位置を元に戻す
		self.lower = (self.lower - y_center) * ratio + y_center	# 中心を0に合わせてから拡大，中心位置を元に戻す

	def y_shift(self):
		center = (np.max(self.upper) + np.min(self.lower)) / 2.0 - 0.5
		self.lower -= center
		self.upper -= center

	def get_mask(self):
		size = self.resolution  # 格子の分割数
		#"""
		if self.resize == True:
			y_u = self.upper
			y_l = self.lower
			self.upper = np.zeros(size)
			self.lower = np.zeros(size)
			# self.upper[int(size/4):int(size/4*3)] = 0.5 * (y_u[0::2] - 0.5) + 0.5 * self.aspect[1]
			# self.lower[int(size / 4):int(size / 4 * 3)] = 0.5 * (y_l[0::2] - 0.5) + 0.5 * self.aspect[1]
			self.upper[int(size / 8 * 3):int(size / 8 * 5)] = 0.25 * (y_u[0::4] - 0.5) + 0.5 * self.aspect[1]
			self.lower[int(size / 8 * 3):int(size / 8 * 5)] = 0.25 * (y_l[0::4] - 0.5) + 0.5 * self.aspect[1]
		#"""
		if self.auto_reshape:
			self.change_aspect()
		self.y_shift()
		mask = np.full((size, size), 1.0) # マスク関数1:流体，0:物体
		# print(mask)
		
		obj_resolution = size  # 物体形状の解像度
		obj_x = np.linspace(start = 0, stop = 1, num = obj_resolution)    # margin分を除去した数列
		obj_dx = 1.0 / (obj_resolution + 1)
		y_u_by_dy = np.floor(self.upper / self.grid_dy)  # y上側点は格子何個分なのか
		y_l_by_dy = np.floor(self.lower / self.grid_dy)  # y下側点は格子何個分なのか
		y_width = self.upper - self.lower   # y上側点とy下側点との間に格子は何個あるのか
		
		grid_x = np.linspace(start = 0, stop = 1, num = size)   # 元格子の格子列

		obj_cell = 0    # 物体形状の格子確認用
		epsilon = 0.01 * self.grid_dy    # 上下の点の間の幅がこの値以下なら格子の存在を認めない
		for grid_cell in range(size):
			if ((grid_cell >= 0.0) or ((grid_cell < (size - 1) - 0.0))):
				if (y_width[obj_cell] > epsilon):
					mask[grid_cell, int(y_l_by_dy[obj_cell]):int(y_u_by_dy[obj_cell] + 1)] = 0
				obj_cell += 1
			
			if (obj_cell == y_width.shape[0]):
				break
				
		self.mask = np.flipud(mask)[-1::-1]
	
	def test_output(self):
		x = np.linspace(start=0, stop=1, num = self.mask.shape[0] + 1)
		z = [0]
		writer = MyStructuredVTKWriter(xcoords = x, ycoords = x, zcoords = z, data_dimension = 2)
		writer.filename = "mask"
		writer.casename = "mask_naca" + self.name
		writer.append_data(name = "object", SCALARS = True, data = self.mask)
		writer.output("mask")
		


def y_upper_test(x, kind="real"):
	if kind == "real":
		return np.where(x < 0.5, 0.5, x)
	elif kind == "int":
		y = np.linspace(start = 0, stop = 1, num = 100)
		return y[int(x)]

def y_upper(x, kind="real"):
	if kind == "real":
		return np.where(x < 0.5, 0.5, x)
	elif kind == "int":
		y = np.linspace(start = 0, stop = 1, num = 100)
		return y[int(x)]


def y_lower(x):
	return np.where(x < 0.5, 0.5, 0.5)


def main():
	naca_4digit_name = "0012"
	attack_angle_deg = 0
	np.set_printoptions(threshold = np.inf)
	grid_resolution = 2**9
	naca = N4d(naca_4digit_name, attack_angle_deg, grid_resolution)

	import matplotlib.pyplot as plt
	print(naca.y_u[-1] - naca.y_l[-1])
	plt.plot(naca.x_u, naca.y_u)
	plt.plot(naca.x_l, naca.y_l)
	plt.show()
	exit()
	aspectratio = 1.0  # y/x
	resize = False # int(grid_resolution / 8)
	s = Shape(naca.equidistant_y_u, naca.equidistant_y_l, grid_resolution, aspectratio, resize, naca_4digit_name)
	s.test_output()
	
	"""
	x = np.linspace(start = 0, stop = 1, num = 100)
	y_u = s.upper(x)
	y_l = s.lower(x)
	
	print(y_u)
	print(y_l)
	"""
	
if __name__ == '__main__':
    main()