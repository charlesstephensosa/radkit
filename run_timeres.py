import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.stats import norm
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from scipy.stats import exponnorm
#import seaborn as sns
#from fitter import Fitter, get_common_distributions, get_distributions

def gauss_flat(xvalues,amplitude,mean,sigma,slope,yintercept,up):
	
	return up + amplitude*np.exp(-(xvalues-mean)**2/(2*sigma**2)) + slope*xvalues + yintercept

def gauss_exp(xvalues,amplitude,mean,sigma,EXAMP,BEXP,up):

	return up + amplitude*np.exp(-(xvalues-mean)**2/(2*sigma**2)) + EXAMP*np.exp(xvalues/BEXP)

class analysis:

	def __init__(self):

		print('time resolution analysis object created')

	def setFrames(self,frameA,frameB):

		print('time resolution analysis: frames are set')

		self.frameA = frameA

		self.frameB = frameB

		self.yframeA_mod = frameA.yframe

		self.yframeB_mod = frameB.yframe

	def filterLLD(self,thresholdA,thresholdB):

		self.yframeA_mod = self.yframeA_mod.T[self.yframeA_mod.max()>thresholdA].T

		self.yframeB_mod = self.yframeB_mod.T[self.yframeB_mod.max()>thresholdB].T

		self.matchIndex()

		print('applied LLD')

	def filterULD(self,thresholdA,thresholdB):

		self.yframeA_mod = self.yframeA_mod.T[self.yframeA_mod.max()<thresholdA].T

		self.yframeB_mod = self.yframeB_mod.T[self.yframeB_mod.max()<thresholdB].T

		self.matchIndex()

		print('applied ULD')
	
	def correlateTime(self):

		print('correlating time...')

		peakposA = self.yframeA_mod.idxmax()
		peakposB = self.yframeB_mod.idxmax()

		meanpeakposA = int(peakposA.mean())
		meanpeakposB = int(peakposB.mean())

		leftboundA = meanpeakposA - 20
		rightboundA = meanpeakposA + 20

		leftboundB = meanpeakposB - 20
		rightboundB = meanpeakposB + 20

		self.yframeA_mod = self.yframeA_mod.T[(peakposA>leftboundA)&(peakposA<rightboundA)].T

		self.yframeB_mod = self.yframeB_mod.T[(peakposB>leftboundB)&(peakposB<rightboundB)].T

		self.matchIndex()

		print('correlating time... done')

	def getIntegrals(self):

		print('calculating integrals...')

		timeStepNS = self.frameA.timeStep*(1e9)

		#sumA = self.yframeA_mod.sum()

		#sumB = self.yframeB_mod.sum()

		self.trackindexA = np.array(self.yframeA_mod.index.tolist())

		self.trackindexB = np.array(self.yframeB_mod.index.tolist())

		self.integralsA = pd.eval('self.yframeA_mod.sum()*timeStepNS',engine='python')

		self.integralsB = pd.eval('self.yframeB_mod.sum()*timeStepNS',engine='python')

		#del sumA

		#del sumB

		print('calculating integrals... done')

	def getIntegrals2(self):

		print('calculating integrals 2...')

		timeStepNS = self.frameA.timeStep*(1e9)

		sumA = self.yframeA_mod.query('self.yframeA_mod.sum()',engine='python')

		sumB = self.yframeA_mod.query('self.yframeB_mod.sum()',engine='python')

		self.integralsA = pd.eval('sumA*timeStepNS')

		self.integralsB = pd.eval('sumB*timeStepNS')

		#del sumA

		#del sumB

		print('calculating integrals... done')

	def getIntegrals3(self):

		print('calculating integrals 3...')

		timeStepNS = self.frameA.timeStep*(1e9)

		self.trackindexA = np.array(self.yframeA_mod.columns.tolist())

		self.trackindexB = np.array(self.yframeB_mod.columns.tolist())

		self.integralsA = pd.eval('self.yframeA_mod.to_numpy().sum(axis=0)*timeStepNS',engine='python')

		self.integralsB = pd.eval('self.yframeB_mod.to_numpy().sum(axis=0)*timeStepNS',engine='python')

		self.integralsA = pd.DataFrame(self.integralsA).set_index(self.trackindexA).T

		self.integralsB = pd.DataFrame(self.integralsB).set_index(self.trackindexB).T

		del timeStepNS

		del self.trackindexA

		del self.trackindexB

		print('calculating integrals... done')

	def correlateEnergy(self,rangeA,rangeB):

		print('coorelating energies...')

		self.corrIntA = self.integralsA.T[np.logical_and((self.integralsA.T>rangeA[0]),(self.integralsA.T<rangeA[1]))].T

		self.corrIntB = self.integralsB.T[np.logical_and((self.integralsB.T>rangeB[0]),(self.integralsB.T<rangeB[1]))].T

		self.correlateEnergyIndex = np.intersect1d(self.corrIntA.T.index,self.corrIntB.T.index)

		self.yframeA_mod = self.yframeA_mod[self.correlateEnergyIndex]

		self.yframeB_mod = self.yframeB_mod[self.correlateEnergyIndex]

		#self.matchIndex()

		print('coorelating energies... done')

	def adjustFrameToCFD(self,CFD):

		self.yframeA_mod = self.yframeA_mod - self.yframeA_mod.max()*CFD

		self.yframeB_mod = self.yframeB_mod - self.yframeB_mod.max()*CFD

		self.matchIndex()

		print('adjusted vertical frame alignment to CFD value')

	def adjustFrameToCFD2(self,CFD):

		#self.CFD_A = self.yframeA_mod.max()*CFD

		#self.CFD_B = self.yframeB_mod.max()*CFD

		


		

		index = self.latestIndex

		peaklocA = self.yframeA_mod.idxmax().to_numpy()
		peaklocB = self.yframeB_mod.idxmax().to_numpy()

		psA = []
		psB = []

		for i in np.arange(0,len(index),1):
			peaksmoothA = self.yframeA_mod[index[i]][peaklocA[i]+0:peaklocA[i]+1].mean()*CFD
			peaksmoothB = self.yframeB_mod[index[i]][peaklocB[i]+0:peaklocB[i]+1].mean()*CFD
			psA.append(peaksmoothA)
			psB.append(peaksmoothB)

		self.CFD_A = np.array(psA)
		self.CFD_B = np.array(psB)

		self.yframeA_cfd = self.yframeA_mod - self.CFD_A

		self.yframeB_cfd = self.yframeB_mod - self.CFD_B

		fixedThreshold = 0.3



		xposA=[]
		xposB=[]

		startXA = []
		startXB = []

		for i in np.arange(0,len(index),1):
			startXA.append(self.yframeA_cfd[index[i]][self.yframeA_cfd[index[i]]>0].index[0])
			startXB.append(self.yframeB_cfd[index[i]][self.yframeB_cfd[index[i]]>0].index[0])

			rangeA = [startXA[i]-1,startXA[i]+1]
			rangeB = [startXB[i]-1,startXB[i]+1]

			xa1 = self.xframeA_mod[index[i]][rangeA[0]:rangeA[1]].to_numpy()

			ya1 = self.yframeA_mod[index[i]][rangeA[0]:rangeA[1]].to_numpy()

			xb1 = self.xframeB_mod[index[i]][rangeB[0]:rangeB[1]].to_numpy()

			yb1 = self.yframeB_mod[index[i]][rangeB[0]:rangeB[1]].to_numpy()

			#interpFA = interpolate.UnivariateSpline(ya1,xa1,s=0)
			#interpFB = interpolate.UnivariateSpline(yb1,xb1,s=0)

			interpFA = interpolate.interp1d(ya1, xa1)
			interpFB = interpolate.interp1d(yb1, xb1)

			interpA = interpFA(self.CFD_A[i])
			interpB = interpFB(self.CFD_B[i])

			#interpA = interpFA(fixedThreshold)
			#interpB = interpFB(fixedThreshold)

			#interpA = interpolate.interp1d(ya, xa)
			#interpB = interpolate.interp1d(yb, xb)

			xposA.append(interpA)
			xposB.append(interpB)



		#self.startXA = startXA
		#self.startXB = startXB

		self.xposA = np.array(xposA).flatten()
		self.xposB = np.array(xposB).flatten()

		print('adjusted vertical frame alignment to CFD value')

	def adjustFrameToCFD3(self,CFD):

		print('processing CFD...')

		index = self.latestIndex

		self.CFD_A = self.yframeA_mod.max().to_numpy()*CFD

		self.CFD_B = self.yframeB_mod.max().to_numpy()*CFD

		FA = self.yframeA_mod

		FB = self.yframeB_mod

		peakAloc = FA.idxmax().to_numpy()

		peakBloc = FB.idxmax().to_numpy()

		timestep = 0.4

		timedomain = np.arange(0,timestep*len(self.yframeA_mod),timestep)

		xposA=[]
		xposB=[]

		for i in np.arange(0,len(index),1):

			self.place = index[i]

			self.FAX = FA[index[i]][peakAloc[i]-10:peakAloc[i]+1]
			self.FAY = timedomain[peakAloc[i]-10:peakAloc[i]+1]

			self.FBX = FB[index[i]][peakBloc[i]-10:peakBloc[i]+1]
			self.FBY = timedomain[peakBloc[i]-10:peakBloc[i]+1]

			#self.fineTimeA = np.arange(0,self.FAY.max()[:-1],0.1)
			#self.fineTimeB = np.arange(0,self.FBY.max()[:-1],0.1)

			self.interpFA = interpolate.interp1d(self.FAX, self.FAY, kind='linear')
			self.interpFB = interpolate.interp1d(self.FBX, self.FBY, kind='linear')

			plt.plot(self.FAX,self.FAY);plt.show()

			self.interpA = self.interpFA(self.CFD_A[i])
			self.interpB = self.interpFB(self.CFD_B[i])

			xposA.append(self.interpA)
			xposB.append(self.interpB)



		self.xposA = np.array(xposA).flatten()
		self.xposB = np.array(xposB).flatten()

		print('processing CFD... done')

	def adjustFrameToCFD4(self,CFD):

		print('processing CFD...')

		FA = self.yframeA_mod

		FB = self.yframeB_mod

		peakAloc = FA.idxmax().to_numpy()

		peakBloc = FB.idxmax().to_numpy()

		index = self.latestIndex

		timestep = 0.4

		timedomain = np.arange(0,timestep*len(self.yframeA_mod),timestep)

		timedomainFine = np.arange(0,timestep*len(self.yframeA_mod),timestep)

		self.yframeA_mod = self.yframeA_mod - self.yframeA_mod.max()*CFD

		self.yframeB_mod = self.yframeB_mod - self.yframeB_mod.max()*CFD

		for i in np.arange(0,len(index),1):

			#interpFA = interpolate.UnivariateSpline(ya1,xa1,s=0)
			#interpFB = interpolate.UnivariateSpline(yb1,xb1,s=0)

			self.FAY = FA[index[i]][peakAloc[i]-10:peakAloc[i]+0].to_numpy()
			self.FAX = timedomain[peakAloc[i]-10:peakAloc[i]+0]

			self.FBY = FB[index[i]][peakBloc[i]-10:peakBloc[i]+50].to_numpy()
			self.FBX = timedomain[peakBloc[i]-10:peakBloc[i]+50]

			self.FAY = np.cumsum(self.FAY)
			self.FBY = np.cumsum(self.FBY)

			self.FAY = self.FAY/self.FAY.max()
			self.FBY = self.FBY/self.FBY.max()

			self.interpFA = interpolate.interp1d(self.FAX, self.FAY)
			self.interpFB = interpolate.interp1d(self.FBX, self.FBY)

			#plt.plot(timedomainFine, self.interpFA(timedomainFine));plt.show()

			#x = self.FAX

			plt.plot(self.FAX,self.FAY);plt.show()

			#self.RA = self.interpFA.roots[0]
			#self.RB = self.interpFB.roots[0]


		self.matchIndex()

		print('adjusted vertical frame alignment to CFD value')

		print('processing CFD... done')

	def getRootsXY(self):

		print('calculating roots')

		roots2 = []
		roots3 = []

		

		self.rootsA = roots2
		self.rootsB = roots3

		print('found roots')

	def getRootsXY2(self):

		print('calculating roots')

		roots2 = []
		roots3 = []

		index = self.latestIndex

		xA = self.frameA.xframe
		xB = self.frameB.xframe

		for i in np.arange(0,len(index),1):
			roots2.append(xA[self.latestIndex[i]][self.startXA[i]])
			roots3.append(xB[self.latestIndex[i]][self.startXB[i]])


		self.rootsA = roots2
		self.rootsB = roots3

		print('found roots')

	def getDelta(self,binWidth):

		print('getting delta...')

		rootsA = self.xposA
		rootsB = self.xposB

		deltaA = rootsA - rootsB
		deltaB = rootsB - rootsA

		meanA,stdA=norm.fit(deltaA)
		meanB,stdB=norm.fit(deltaB)

		if (meanA>0):
			thedelta = deltaA

		if (meanB>0):
			thedelta = deltaB

		self.delta = thedelta

		meanDelta,stdDelta=norm.fit(self.delta)

		self.deltaRange = np.arange(meanDelta-3.0*stdDelta,meanDelta+3.0*stdDelta,binWidth)

		print('getting delta... done')




	def matchIndex(self):

		self.yframeA_mod_index = self.yframeA_mod.columns.tolist()

		self.yframeB_mod_index = self.yframeB_mod.columns.tolist()

		self.latestIndex = np.intersect1d(self.yframeA_mod_index,self.yframeB_mod_index)

		self.yframeA_mod = self.yframeA_mod[self.latestIndex]

		self.yframeB_mod = self.yframeB_mod[self.latestIndex]

	
	def openplot(self):

		plt.figure()
		plt.grid('on')
		plt.xlabel('channels')
		plt.ylabel('counts')
		plt.plot(self.xaxis,self.yaxis,'--k')

		figManager = plt.get_current_fig_manager()
		figManager.window.state('zoomed')

		peaks = plt.ginput(n=-1,show_clicks=True,mouse_add=None,mouse_pop=None,mouse_stop=None,timeout=0)

		self.peaks = np.array(peaks)

		plt.scatter(self.peaks[:,0],self.peaks[:,1],c='blue',marker='o',linewidths=5)

		bases = plt.ginput(n=-1,show_clicks=True,mouse_add=None,mouse_pop=None,mouse_stop=None,timeout=0)

		self.bases = bases

	def guessParameters(self):

		self.centroids = self.peaks[:,0][0]
		self.amplitudes = self.peaks[:,1][0]
		self.expheight = self.bases[0][1]
		self.sigmas = (self.bases[2][0]-self.bases[1][0])/2.355
		self.run = self.bases[3][0]-self.bases[0][0]
		self.decay = -1*self.run
		self.rise = self.bases[0][1]-self.bases[3][1]
		self.slope = -1*self.rise/self.run
		self.yintercept = self.bases[0][1] + self.slope*self.bases[0][0]
		self.up = 0
	
		if (self.fitType == 0):

			self.guess = [self.amplitudes,self.centroids,self.sigmas,self.slope,self.yintercept,self.up]

		if (self.fitType == 1):

			self.guess = [self.amplitudes,self.centroids,self.sigmas,self.expheight,self.decay,self.up]

		print(self.guess)

		self.xnew = np.arange(self.bases[0][0],self.bases[3][0],1)
		self.ynew = self.interp(self.xnew)

		plt.plot(self.xnew,self.ynew);plt.show()

	def doFits(self):

		if (self.fitType == 0):

			self.fitFunction = curve_fit(gauss_flat,self.xnew,self.ynew,p0=self.guess)

			self.yfit = gauss_flat(self.xnew,*self.fitFunction[0])

		if (self.fitType == 1):

			self.fitFunction = curve_fit(gauss_exp,self.xnew,self.ynew,p0=self.guess)

			self.yfit = gauss_exp(self.xnew,*self.fitFunction[0])

	def plotDataWithFit(self):

		plt.plot(self.xnew,self.ynew)
		plt.plot(self.xnew,self.yfit)
		plt.show()



	

		







