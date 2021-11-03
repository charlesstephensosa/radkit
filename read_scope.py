import glob
import numpy as np
import pandas as pd
from tqdm import tqdm
import lecroyutils.data
import warnings
#import numba

warnings.filterwarnings("ignore")

class data:

	def __init__(self):

		print('scope object created')

		self.lecroy = lecroyutils.data
		self.lecroy.warn('false')

	def setDataDirectory(self,directory,numFiles,channel):

		self.numFiles = numFiles

		self.directory = directory + channel + '*'

		print('directory set')

	def getDataInformation(self):

		self.paths = glob.glob(self.directory)

		self.paths = self.paths[0:self.numFiles]
		
		self.lecroy.warn('false')

		self.sampleFile = self.lecroy.LecroyScopeData.parse_file(self.paths[0])

		self.verticalRange = self.sampleFile.y_max - self.sampleFile.y_min

		self.segmentsize = self.sampleFile.y.shape[0]

		self.timeStep = self.sampleFile.Ts

		self.basePoints = int(0.10*self.segmentsize)

		self.baseMean = self.sampleFile.y[0:self.basePoints].mean()

		self.pulseStateID = 1

		if (self.baseMean>0):

			print('detected positive base')

			self.baseState = 'positive'

			self.baseMean = -1*self.baseMean

		else:

			print('detected negative base')

			self.baseState = 'negative'

		if (np.abs(self.sampleFile.y_min)<np.abs(self.sampleFile.y_max)):

			print('detected positive pulses')

			self.pulseState = 'positive'

			self.pulseStateID = 1

		else:

			print('detected negative pulses')

			self.pulseState = 'negative'

			self.pulseStateID = -1

		print('retrieved data information')

		del self.sampleFile


	def parseDataPulses(self):

		print('data reading and parsing started')

		framelist = []

		for i, p in (enumerate(tqdm(self.paths))):

			self.lecroy.warn('false')

			file = self.lecroy.LecroyScopeData.parse_file(p)

			framelist.append(pd.DataFrame(file.y))

		del file

		del self.lecroy

		print('data reading and parsing finished')

		self.numFiles = i + 1

		print('pandas dataframe construction started')

		self.yframe = pd.concat(framelist,axis=1,ignore_index=True)

		del framelist

		print('pandas dataframe construction finished')


	def adjustFrames(self):

		print('pandas dataframe adjustment started')

		if (self.baseState == 'positive'):

			print('correcting positive baseline')

			self.yframe = pd.eval('self.pulseStateID*(self.yframe+self.baseMean)')

		if (self.baseState == 'negative'):

			print('correcting negative baseline')

			self.yframe = pd.eval('self.pulseStateID*(self.yframe+self.baseMean)')

		print('pandas dataframe adjustment finished')


def filterLLD(f,threshold):

	o = f.shape[1]

	f = f.T[f.max()>threshold].T

	i = f.columns.tolist()

	return f, i

def filterULD(f,threshold):

	o = f.shape[1]

	f = f.T[f.max()<threshold].T

	i = f.columns.tolist()

	return f, i

def matchIndex(f2,f3,i2,i3):

	i = np.intersect1d(i2,i3)

	f2 = f2[i]

	f3 = f3[i]

	k = i.size

	return f2, f3, i

		

