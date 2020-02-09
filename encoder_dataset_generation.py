import os
import argparse
import numpy as np
from grid_generation import sin, create_grid


def main():
	"""
	Creates a dataset for encoder, used to compress the description of the surfaces
	Dataset folder is called 'encoder' and consists of three subfolders:
	'train', 'test', and 'val'. Each of it contains csv files , numerated from 0 to 'dataset_size'

	"""

	parser = argparse.ArgumentParser()
	parser.add_argument("train_size" , help="size of the train dataset", type=int)
	parser.add_argument("val_size"   , help="size of the val dataset"  , type=int)
	parser.add_argument("test_size"  , help="size of the test dataset" , type=int)	
	
	args = parser.parse_args()

	root_dir  = os.path.join(os.getcwd(), 'encoder')
	W           = 100.
	H           = 50.

	try:
		os.mkdir(root_dir)
	except FileExistsError:
		print("dataset already exists")
		return

	def fill_dir(dir_name, dir_size):
		
		for i in range(dir_size):
		
			loc   = np.random.uniform(10., 18.)
			w_num = np.random.uniform(3., 7.)
			amp   = np.random.uniform(0.6, 2.)
			curve = lambda x: sin(x, amp, w_num, H - loc, W)

			zz, curve_points = create_grid(curve, W, H)
			np.savetxt(
				os.path.join(dir_name, f"{i}.csv"), 
				zz.reshape([128, 128]), delimiter=",")

	# train dir

	train_size = args.train_size 
	train_dir  = os.path.join(root_dir, 'train')
	os.mkdir(train_dir)
	fill_dir(train_dir, train_size)

	# val dir

	val_size = args.val_size 
	val_dir  = os.path.join(root_dir, 'val')
	os.mkdir(val_dir)
	fill_dir(val_dir, val_size)

	# test dir 

	test_size = args.test_size 
	test_dir  = os.path.join(root_dir, 'test')
	os.mkdir(test_dir)
	fill_dir(test_dir, test_size)


if __name__ == '__main__':
	main()