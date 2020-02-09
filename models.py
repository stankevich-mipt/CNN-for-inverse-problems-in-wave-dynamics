import torch.nn as nn
import torch.nn.functional as F


class UNetEncoder(nn.Module):
	"""
	This architecture is heavily based on the UNet, but simplified due to the
	peculiarity of the problem.
	See https://arxiv.org/pdf/1505.04597.pdf for the original article

	"""

	def __init__(self, activation='ReLU', pooling='max'):
	  
		super(UNet, self).__init__()

		if activation == 'ReLU':
			activation_function  = nn.ReLU()
		elif activation == 'LeakyReLU':
			activation_function  = nn.LeakyReLU()
		else:
			raise NotImplementedError
	
		if pooling == 'max':
			pooling_layer = nn.MaxPool2d(kernel_size=2, stride=2)
		else:
			raise NotImplementedError

		# ----------------- Encoder -------------------

		self.conv1 = nn.Sequential(
			nn.Conv2d(1, 32, kernel_size=3, padding=1),
			nn.BatchNorm2d(32),
			activation_function,
			nn.Conv2d(32, 32, kernel_size=3, padding=1),
			nn.BatchNorm2d(32),
			activation_function,
			pooling_layer
		)
		   
		self.conv2 = nn.Sequential(
			nn.Conv2d(32, 64, kernel_size=3, padding=1),
			nn.BatchNorm2d(64),
			activation_function,
			nn.Conv2d(64, 64, kernel_size=3, padding=1),
			nn.BatchNorm2d(64),
			activation_function,
			pooling_layer,
			nn.Dropout2d(p=0.3)
		)

		self.conv3 = nn.Sequential(
			nn.Conv2d(64, 128, kernel_size=3, padding=1),
			nn.BatchNorm2d(128),
			activation_function,
			nn.Conv2d(128, 128, kernel_size=3, padding=1),
			nn.BatchNorm2d(128),
			activation_function,
			pooling_layer,
			nn.Dropout2d(p=0.3)
		)

		self.conv4 = nn.Sequential(
		nn.Conv2d(128, 1, kernel_size=3, padding=1),
		)

		self.encoder = nn.Sequential(
			self.conv1,
			self.conv2,
			self.conv3,
			self.conv4
		)

		# ---------------------------------------------
		# ----------------- Decoder -------------------

		self.conv5 = nn.Sequential(
			nn.Conv2d(1, 128, kernel_size=3, padding=1),
			nn.BatchNorm2d(128),
			activation_function
		)

		self.conv6 = nn.Sequential(
			nn.Upsample(scale_factor=2, mode='bilinear'),
			nn.Conv2d(128, 64, kernel_size=3, padding=1),
			nn.BatchNorm2d(64),
			activation_function,
			nn.Conv2d(64, 64, kernel_size=3, padding=1),
			nn.BatchNorm2d(64),
			activation_function,
			nn.Dropout2d(p=0.3)
		)

		self.conv7 = nn.Sequential(
			nn.Upsample(scale_factor=2, mode='bilinear'),
			nn.Conv2d(64, 32, kernel_size=3, padding=1),
			nn.BatchNorm2d(32),
			activation_function,
			nn.Conv2d(32, 32, kernel_size=3, padding=1),
			nn.BatchNorm2d(32),
			activation_function,
			nn.Dropout2d(p=0.3)
		)

		self.conv8 = nn.Sequential(
			nn.Upsample(scale_factor=2, mode='bilinear'),
			nn.Conv2d(32, 32, kernel_size=3, padding=1),
			nn.BatchNorm2d(32),
			activation_function,
			nn.Conv2d(32, 1, kernel_size=3, padding=1)
		)

		self.decoder = nn.Sequential(
			self.conv5,
			self.conv6,
			self.conv7,
			self.conv8
		)

	def encode(self, x):
		h = self.encoder(x)
		return h

	def decode(self, h):
		x = self.decoder(h)
		return x
	  
	def forward(self, x):	
		return self.decode(self.encode(x))


class UNetEncoder2(nn.Module):
	"""
	This architecture is a lightweight version of the previous one
	See https://arxiv.org/pdf/1505.04597.pdf for the original article

	"""
	def __init__(self, activation='ReLU', pooling='max'):
	  
		super(UNet2, self).__init__()

		if activation == 'ReLU':
			activation_function  = nn.ReLU()
		elif activation == 'LeakyReLU':
			activation_function  = nn.LeakyReLU()
		else:
			raise NotImplementedError
		
		if pooling == 'max':
			pooling_layer = nn.MaxPool2d(kernel_size=2, stride=2)
		else:
			raise NotImplementedError
	
	# ----------------- Encoder -------------------

		self.conv1 = nn.Sequential(
			nn.Conv2d(1, 16, kernel_size=3, padding=1),
			nn.BatchNorm2d(16),
			activation_function,
			nn.Conv2d(16, 16, kernel_size=3, padding=1),
			nn.BatchNorm2d(16),
			activation_function,
			pooling_layer
		)
		   
		self.conv2 = nn.Sequential(
			nn.Conv2d(16, 32, kernel_size=3, padding=1),
			nn.BatchNorm2d(32),
			activation_function,
			nn.Conv2d(32, 32, kernel_size=3, padding=1),
			nn.BatchNorm2d(32),
			activation_function,
			pooling_layer,
			nn.Dropout2d(p=0.3)
		)

		self.conv3 = nn.Sequential(
			nn.Conv2d(32, 64, kernel_size=3, padding=1),
			nn.BatchNorm2d(64),
			activation_function,
			nn.Conv2d(64, 64, kernel_size=3, padding=1),
			nn.BatchNorm2d(64),
			activation_function,
			pooling_layer,
			nn.Dropout2d(p=0.3)
		)

		self.conv4 = nn.Sequential(
			nn.Conv2d(64, 1, kernel_size=3, padding=1),
		)

		self.encoder = nn.Sequential(
			self.conv1,
			self.conv2,
			self.conv3,
			self.conv4
		)
		
		# ---------------------------------------------
		# ----------------- Decoder -------------------

		self.conv5 = nn.Sequential(
			nn.Conv2d(1, 64, kernel_size=3, padding=1),
			nn.BatchNorm2d(64),
			activation_function
		)
		
		self.conv6 = nn.Sequential(
			nn.Upsample(scale_factor=2, mode='bilinear'),
			nn.Conv2d(64, 32, kernel_size=3, padding=1),
			nn.BatchNorm2d(32),
			activation_function,
			nn.Conv2d(32, 32, kernel_size=3, padding=1),
			nn.BatchNorm2d(32),
			activation_function,
			nn.Dropout2d(p=0.3)
		)
		
		self.conv7 = nn.Sequential(
			nn.Upsample(scale_factor=2, mode='bilinear'),
			nn.Conv2d(32, 16, kernel_size=3, padding=1),
			nn.BatchNorm2d(16),
			activation_function,
			nn.Conv2d(16, 16, kernel_size=3, padding=1),
			nn.BatchNorm2d(16),
			activation_function,
			nn.Dropout2d(p=0.3)
		)
		
		self.conv8 = nn.Sequential(
			nn.Upsample(scale_factor=2, mode='bilinear'),
			nn.Conv2d(16, 16, kernel_size=3, padding=1),
			nn.BatchNorm2d(16),
			activation_function,
			nn.Conv2d(16, 1, kernel_size=3, padding=1)
		)

		self.decoder = nn.Sequential(
			self.conv5,
			self.conv6,
			self.conv7,
			self.conv8
		)

	def encode(self, x):
		h = self.encoder(x)
		return h

	def decode(self, h):
		x = self.decoder(h)
		return x
	  
	def forward(self, x):	
		return self.decode(self.encode(x))