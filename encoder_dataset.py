from torch.utils.data import Dataset
from torchvision import transforms
import os
import re
import torch


class GridDataset(Dataset):
    """
    This class is implementing an interface to interact with pytorch Dataloaders

    Format for grid data is 'csv'
    Names of the files must contain numerals, due to the indexing issues
    File with the least numeral is considered to be a 'head' of the dataset
    """

    def __init__(self, root_dir, transform=None):
        """
        :params:
            root_dir (string): Directory with the images.
            transform (callable, optional): Optional transform to be applied
                on a sample.
        """
        self.root_dir    = os.path.join(os.getcwd(), root_dir)
        self.transform   = transform
        
        self.dir_content_data = [] # list with all grid filenames
        for file in os.listdir(self.root_dir):
            self.dir_content_data.append(
                (file, re.findall(r'(\d+)', file)[-1]))

        self.dir_content_data = [
            x[0] for x in list(sorted(self.dir_content_data, key = lambda x: x[1]))]
        
    def __len__(self):
        return len(self.dir_content_data)

    def __getitem__(self, idx):
        """
        Method allows to do 'online' loading from disk to save RAM
        """
        if torch.is_tensor(idx):
            idx = idx.tolist()

        grid_name = os.path.join(self.root_dir,
                                 self.dir_content_data[idx])

        data = np.genfromtxt(grid_name, delimiter=",")

        if self.transform:
            data = self.transform(data)

        return data


def create_encoder_dataset_loaders(batch_size: int):
 	"""
	:returns:
		train_batch_generator - pytorch Dataloader	
		val_batch_generator   - pytorch Dataloader
		test_batch_generator  - pytorch Dataloader	
	See https://pytorch.org/docs/stable/data.html?highlight=dataloader#torch.utils.data.DataLoader
	for more detailed description of dataloaders
	Function requires the dataset to be layed in the working directory in the 'encoder' folder
 	"""

 	transformer = transforms.Compose([  
    	transforms.Lambda(lambda x: torch.from_numpy(x).view(1, x.shape[0], x.shape[1]))
	])

 	root_dir  = os.path.join(os.getcwd(), 'encoder')
 	train_dir = os.path.join(root_dir, 'train')
 	val_dir   = os.path.join(root_dir, 'val')
 	test_dir  = os.path.join(root_dir, 'test')

 	if not os.direxists(root_dir):
 		raise Exception("Dataset does not exists")

 	train_dir = os.path.join(root_dir, 'train')
 	if not (os.direxists(train_dir) and os.direxists(val_dir) and os.direxists(test_dir)):
 		raise Exception("Corrupted dataset layout")


 	train_dataset = GridDataset(train_dir, transform=transformer)
 	val_data      = GridDataset(val_dir,   transform=transformer)
 	test_data     = GridDataset(test_dir,  transform=transformer)

 	train_batch_generator = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
 	val_batch_generator   = torch.utils.data.DataLoader(val_dataset,   batch_size=batch_size, shuffle=True)
 	test_batch_generator  = torch.utils.data.DataLoader(test_dataset,  batch_size=batch_size, shuffle=True)

 	return train_batch_generator, val_batch_generator, test_batch_generator 



