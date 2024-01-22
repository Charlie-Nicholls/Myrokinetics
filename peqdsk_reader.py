from numpy import nan, array
from scipy.interpolate import InterpolatedUnivariateSpline

class peqdsk(object):
	def __init__(self, filename = None, directory = "./"):
		import os.path
		self.filename = os.path.join(directory, filename)
		self.data = self.read(self.filename)
		
	def __getitem__(self, key):
        	return self.data[key]

	def read(self, filename):
		kinetics_file = open(filename)
		lines = kinetics_file.readlines()
		kinetics_file.close()
		kinetics_data = {}
		ni = 0
		while True:
			info = lines[ni].split(" ")
			n = eval(info[0])
			name1 = info[1]
			name2 = info[2].split("(")[0]
			name3 = info[3].strip("\n")
			list1 = []
			list2 = []
			list3 = []
			ni = ni+1
			for i in range(ni,n+ni):
				line = list(filter(None, lines[i].split(" ")))
				if 'nan' not in line[0]:
					list1.append(eval(line[0].strip(" \n")))
				else:
					list1.append(nan)
				if 'nan' not in line[1]:
					list2.append(eval(line[1].strip(" \n")))
				else:
					list2.append(nan)
				if 'nan' not in line[2]:
					list3.append(eval(line[2].strip(" \n")))
				else:
					list3.append(nan)
			list1 = array(list1)
			list2 = array(list2)
			list3 = array(list3)
			if name1 in kinetics_data.keys():
				f1 = InterpolatedUnivariateSpline(list1, list2)
				f2 = InterpolatedUnivariateSpline(list1, list3)
				kinetics_data[name2] = f1(kinetics_data[name1])
				kinetics_data[name3] = f2(kinetics_data[name1])
			else:
				kinetics_data[name1] = list1
				kinetics_data[name2] = list2
				kinetics_data[name3] = list3
			ni = ni+n
			#print(name1, name2, name3)
			if ni >= len(lines):
				break
		return kinetics_data
	
	def keys(self):
		return self.data.keys()

