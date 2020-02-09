import numpy as np
import math  as m
import argparse

def sin(x: np.ndarray, amp: float, w_num: float, loc: float, W: float):
    return loc + amp * np.sin(2 * m.pi * x / (W / w_num))


def cos(x: np.ndarray, amp: float, w_num: float, loc: float, W: float):
    return loc + amp * np.cos(2 * m.pi * x / (W / w_num))


def create_grid(curve, W: float, H: float):
	"""
	:params:
		curve: function, describing a curve
		W: width of the base box
		H: heigth of the base box
	:returns:
		zz: np.array
		curve_points: np.array
    The function creates description of a two 2d surfaces, created from 
    rectangular box with a curve given with by 'curve' argument.
    Description is parametrized with a regular mesh of size 128x128, with values {1, 0} for 
    surfaces respectfully. 
	"""

	curve_points = np.array([[a, b] for a, b in zip(
		np.linspace(0, W, 101), [curve(x) for x in np.linspace(0, W, 101)])
	])  

	xx, yy = np.meshgrid(np.linspace(0, W, 128),
	                     np.linspace(0, H, 128))
	zz = []
	for x, y in zip(xx.flatten(), yy.flatten()):
		zz.append(1. if (H - curve(x)) >= y else 0.)

	zz = np.array(zz)
	return zz, curve_points


def main():
	"""
	Aggregates grids descriptions into single csv file (see create_grid docstring)
	Each curve description is also written in separate csv file
	"""

	parser = argparse.ArgumentParser()
	parser.add_argument("batch_index", help="Index of the current batch"  , type=int)
	parser.add_argument("batch_size" , help="Number of grids in the batch", type=int)
	parser.add_argument("header"     , help="a prefix for the grid names" , type=str)
	args = parser.parse_args()
	
	header      = args.header
	batch_size  = args.batch_size
	batch_index = args.batch_index

	W           = 100.
	H           = 50.

	res = np.ones(128 * 128)

	for i in range(batch_size):

		loc   = np.random.uniform(10., 18.)
		w_num = np.random.uniform(3., 7.)
		amp   = np.random.uniform(0.6, 2.)
		
		curve = lambda x: sin(x, amp, w_num, H - loc, W)
		zz, curve_points = create_grid(curve, W, H)
		np.savetxt(header+f"_{i}.csv", curve_points, delimiter=",")
		res = np.vstack((res, zz))

	np.savetxt(f"{header}Batch{batch_index}.csv", np.delete(res, 0, 0), delimiter=",")


if __name__ == '__main__':
	main()