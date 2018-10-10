# coding: utf-8
import scipy.linalg as linalg
import scipy.sparse.linalg as spla
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

# 複素数列z1, z2,...を, 2倍の長さの実数列x1, y1, x2, y2, ...に変換
def line_up_z2xy(z, two_rows=False):
	"""
	x = np.real(z).reshape(1, -1)
	y = np.imag(z).reshape(1, -1)
	xy = np.vstack((x, y)).T.reshape(-1)
	return xy
	"""
	if two_rows:
		return np.vstack((np.real(z).reshape(1, -1), np.imag(z).reshape(1, -1))).T
	else:
		return np.vstack((np.real(z).reshape(1, -1), np.imag(z).reshape(1, -1))).T.reshape(-1)

"""
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



def main():
	a = 1
	z1, size = get_complex_coords(type = 3, naca4 = "0012", size = 100)
	grid = np.zeros((size-1, size-1, 2))
	grid[0] = line_up_z2xy(z1[:size-1], two_rows = True)
	"""
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