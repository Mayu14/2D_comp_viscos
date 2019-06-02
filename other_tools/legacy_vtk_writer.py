# coding: utf-8
import os
import numpy as np
# pythonのvtkライブラリの使い方が分からなくてキレたため自分で書いた()
class MyStructuredVTKWriter(object):
	def erase_data(self):
		self.data_name = []
		self.data_kind = []
		self.data_value = []

	def __init__(self, xcoords, ycoords, zcoords, data_dimension):
		# xcoords等は点の座標であるため，セル数 + 1で作成すること
		self.xcoords = xcoords
		self.ycoords = ycoords
		self.zcoords = zcoords
		self.points = [len(xcoords), len(ycoords), len(zcoords)]
		
		self.data_dimension = data_dimension
		self.filename = "untitled"
		self.casename = "untitled"
		
		self.erase_data()
		if data_dimension == 1:
			self.cell = [len(xcoords) - 1, len(ycoords), len(zcoords)]
			self.totalcell = self.cell[0]
		elif data_dimension == 2:
			self.cell = [len(xcoords) - 1, len(ycoords) - 1, len(zcoords)]
			self.totalcell = (self.cell[0]) * (self.cell[1])
		else:
			self.cell = [len(xcoords) - 1, len(ycoords) - 1, len(zcoords) - 1]
			self.totalcell = (self.cell[0]) * (self.cell[1]) * (self.cell[2])
		
		self.header = "# vtk DataFile Version 3.0\n" + self.casename + "\nASCII\nDATASET RECTILINEAR_GRID\nDIMENSIONS " + str(
			self.points[0]) + " " + str(self.points[1]) + " " + str(self.points[2]) + "\n"
			
		
	def append_data(self, name, SCALARS, data):
		self.data_name.append(name)
		self.data_kind.append(SCALARS)
		if SCALARS == True: # kind = "SCALARS"
			data_dimension = 1
		else:   # kind = "VECTORS"
			data_dimension = 3

		if(self.data_dimension == 2):
			if (SCALARS == True):
				self.data_value.append(data.T.reshape(self.cell[0], self.cell[1], self.cell[2], 1))
			else:
				self.data_value.append(np.concatenate([data[:, :, 0].T.reshape(self.cell[0], self.cell[1], self.cell[2], 1),
				                      data[:, :, 1].T.reshape(self.cell[0], self.cell[1], self.cell[2], 1),
				                      np.zeros(self.cell[0] * self.cell[1]).reshape(self.cell[0], self.cell[1], self.cell[2], 1)], axis = 3))

		elif(self.data_dimension != 2):
			print("工事中")
			exit()
		"""
		if ((SCALARS == True) or (self.data_dimension == 3)):
			self.data_value.append(data.reshape(self.cell[0], self.cell[1], self.cell[2], data_dimension))    # points(i, j, k), direction(3)
		
		elif self.data_dimension == 2:
			data = data.reshape(self.cell[0], self.cell[1], self.cell[2], self.data_dimension)
			empty = np.zeros(self.cell[0]*self.cell[1]*self.cell[2]*1).reshape(self.cell[0], self.cell[1], self.cell[2], 3-self.data_dimension)
			self.data_value.append(np.concatenate([data, empty], axis = 3))
		"""
	
	def set_fname(self, fname):
		self.fname = fname
		
	def set_casename(self, casename):
		self.casename = casename

	
	def output(self, folder, TimeSeries=None):
		def auto_rename():
			flag = 0
			count = 0
			fname = self.filename
			self.filename = folder + fname
			while flag == 0:
				if os.path.exists(self.filename + ".vtk"):
					self.filename = folder + fname + str(count).zfill(4)
					count += 1
				else:
					flag = 1
		
		if TimeSeries == None:
			auto_rename()
		else:
			self.filename = folder + self.fname + str(TimeSeries).zfill(8) + ".vtk"
		
		with open(self.filename + ".vtk", 'w') as f:
			f.write(self.header)
			f.write("X_COORDINATES " + str(self.points[0]) + " float\n")
			[f.write(str(self.xcoords[i]) + "\n") for i in range(self.points[0])]
			f.write("Y_COORDINATES " + str(self.points[1]) + " float\n")
			[f.write(str(self.ycoords[i]) + "\n") for i in range(self.points[1])]
			f.write("Z_COORDINATES " + str(self.points[2]) + " float\n")
			[f.write(str(self.zcoords[i]) + "\n") for i in range(self.points[2])]
			
			f.write("CELL_DATA " + str(self.totalcell) + "\n")

			plot_scalar = lambda i, j, k: f.write(str(value[i, j, k, 0]) + "\n")
			plot_vector = lambda i, j, k: f.write(str(value[i, j, k, 0]) + " " + str(value[i, j, k, 1]) + " " + str(value[i, j, k, 2]) + "\n")
			
			for series in range(len(self.data_name)):
				name = self.data_name[series]
				kind = self.data_kind[series]
				value = self.data_value[series]
				
				if kind == True:
					f.write("SCALARS" + " " + name + " float\nLOOKUP_TABLE default\n" )
					# [plot_scalar(i, j ,k) for k in range(self.cell[2]) for j in range(self.cell[1]) for i in range(self.cell[0])]
					[plot_scalar(i, j, k) for i in range(self.cell[0]) for j in range(self.cell[1]) for k in range(self.cell[2])]
					"""
					for i in range(self.cell[0]):
						for j in range(self.cell[1]):
							for k in range(self.cell[2]):
								plot_scalar(i, j, k)
					"""
				else:   # kind == "VECTORS":
					f.write("VECTORS" + " " + name + " float\n")
					# [plot_vector(i, j, k) for k in range(self.cell[2]) for j in range(self.cell[1]) for i in range(self.cell[0])]
					[plot_vector(i, j, k) for i in range(self.cell[0]) for j in range(self.cell[1]) for k in range(self.cell[2])]
					"""
					for i in range(self.cell[0]):
						for j in range(self.cell[1]):
							for k in range(self.cell[2]):
								plot_vector(i, j, k)
					"""

		self.erase_data()

def main():
	pass
	
if __name__ == '__main__':
	main()
 