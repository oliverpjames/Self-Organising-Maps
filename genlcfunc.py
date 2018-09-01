import numpy as np
import matplotlib.pyplot as plt 
import scipy.signal as sg
from astropy.stats import LombScargle 
# program that defines functions that can be used in other scripts to create, fold, bin, and normalise light curves

def genlc(t, functype, freq=1./20000, a=10000, n=100, err=0.05, stt=86400, stop=28800, sec=0, interval=100, day=86400*10, window=True, haserr=True):
	"""generate a light curve of form functype
	if functype is binary1, binary2 or transit then require a to be chosen, must be int and divisible by 10, default is 100000
	else frequency must be chosen, must be float, default is 1/20000
	
	returns time (t), and flux (y) data
	"""
	
	if functype == 'sin':
		y = np.sin(freq*t*(2*np.pi))																												#sine wave
		
	elif functype == 'sawtooth':
		y = sg.sawtooth(2*np.pi*freq*t)																												#sawtooth function
		
	elif functype == 'rsawtooth':
		y = sg.sawtooth(2*np.pi*freq*t, 0)																											#reverse sawtooth function
		
	elif functype == 'triangle':
		y = sg.sawtooth(2*np.pi*freq*t, 0.5)																										#triangle function
		
	elif functype == 'square':
		y = sg.square(2*np.pi*freq*t)																												#square wave
		
		
	elif functype == 'binary1':
		b = a/10																																	#choose feature start point (seconds)
		y = t*0 + 1 																																#create an array of length of t where each number is equal to 1
		for i in range(n):
			y[(a*i + b)/interval : (a*i + 2*b)/interval] = y[(a*i + b)/interval : (a*i + 2*b)/interval]	  -	np.linspace(0,0.02,b/interval)			#give a v shape dip of a binary (2% dip in light)
			y[(a*i + 2*b)/interval : (a*i + 3*b)/interval] = y[(a*i + 2*b)/interval : (a*i + 3*b)/interval] -	np.linspace(0.02,0,b/interval)	
	
	elif functype == 'binary2':
		b = a/10																																	#choose feature start point (seconds)
		y = t*0 + 1 																																#create an array of length of t where each number is equal to 1
		for i in range(n):
			y[(a*i + b)/interval : (a*i + 2*b)/interval] = y[(a*i + b)/interval : (a*i + 2*b)/interval]	  -	np.linspace(0,0.02,b/interval)			#give a v shape dip of a binary (2% dip in light)
			y[(a*i + 2*b)/interval : (a*i + 3*b)/interval] = y[(a*i + 2*b)/interval : (a*i + 3*b)/interval] -	np.linspace(0.02,0,b/interval)	
			
			y[(a*i + 6*b)/interval : (a*i + 7*b)/interval] = y[(a*i + 6*b)/interval : (a*i + 7*b)/interval]	  -	np.linspace(0,0.01,b/interval)		#give a 2nd v shape dip of a binary (1% dip in light)
			y[(a*i + 7*b)/interval : (a*i + 8*b)/interval] = y[(a*i + 7*b)/interval : (a*i + 8*b)/interval] -	np.linspace(0.01,0,b/interval)	
	
	elif functype == 'transit':
		b = a/10																																	#choose feature start point (seconds)
		y = t*0 + 1 																																#create an array of length of t where each number is equal to 1
		for i in range(n):
			y[(a*i + b)/interval : (a*i + 2*b)/interval] = y[(a*i + b)/interval : (a*i + 2*b)/interval]	  -	np.linspace(0.02,0.02,b/interval)		#give a u shape dip of a transit (2% dip in light)
			y[(a*i + 2*b)/interval : (a*i + 3*b)/interval] = y[(a*i + 2*b)/interval : (a*i + 3*b)/interval] -	np.linspace(0.02,0.02,b/interval)	

	

	if window == True:
		dellist = []																		#create empty list
	
		for i in range(len(t)*interval/day+1):												#for each set of 24 hours
			dellist.append(range((stop + sec)/interval, (stt + sec)/interval))				#add a list of hours to be deleted
			sec = sec + day																	#adds time (eg a day) to the total time passed 
		dellist = np.reshape(dellist, len(dellist)*(stt-stop)/interval)						#reshape dellist to be a list of integers rather than a list of lists
  	
		t = np.delete(t, dellist)
		y = np.delete(y, dellist)

	if haserr == True:
		for i in range(len(y)):
			y[i] = y[i] + np.random.normal(0, err)											#add normally distributed errors to y

	return t, y
	
	


def fold(t, y, periodogram=False):
	"""Folds data on pattern frequency
	
	if periodogram = True then returns T, frequency, power
	otherwise by default
	returns T only"""
	
	frequency, power = LombScargle(t, y).autopower()				#find frequencies that contribute to describe data 
	fpat = frequency[np.argmax(power)] 								#find the most important frequency, this will likely be the period of the orbit
	

	#create array of folded times T
	T = t*0															#initialize array for folded times T 
	c = 0															#initialize c the multiplication factor to choose the bin in which each datapoint goes

	for i in range (len(t)):										#iterate over all numbers in t
		while (t[i] >= c/fpat):										#if the time t[i] is larger than the current multiple of the period  (=1/pattern frequency) 
			c = c + 1												#	then add 1 to the multiplication factor c until t[i] is smaller than c*period
		T[i] = t[i] - c/fpat										#remove c-1 periods from the time 
		
	if periodogram == True:
		return T, frequency, power
	else:
		return T
	

def norm(T, y):
	"""normalises T and y
	
	Returns Tnorm, ynorm"""	
	
	Tnorm = T/max(T)
	
	ynorm = y - min(y)
	ynorm = ynorm/max(ynorm)

	return Tnorm, ynorm




def bindata(T, y, bins=10):
	"""puts data into bins
	default is number of bins is 10
	
	returns arrays of bin times, means, and standard deviations
	"""
	
	tbin = np.linspace(min(T), max(T), bins+1) 					#split the bins into equal steps (+1 allows for drawing of all bins)								
	binindex = np.zeros(len(T), dtype=int)						#initialize an array of zeroes the same length as T

	for i in range (len(T)):	
		for  j in range(bins):									#for each bin
			if (T[i] >= tbin[j]) and (T[i] < tbin[j + 1]):		#	if the time is in bin i
				binindex[i] = j									#		note the index of the bin
				break											#		stop testing if it is in other bins



	#Find bin errors
	powlst = [[] for i in xrange(bins)]							#create list of bins lists (for the number of times bins, powlst = [] and list them)


	for i in range(len(binindex)):
		powlst[binindex[i]].append(y[i])						#append binindex[i] list of powers with the values of y that should be in that bin


	meanpow = np.zeros(bins)									#create an array of the same length of tbins in which to hold the mean of the power in those bins
	stdpow = np.zeros(bins)

	for i in range(bins):
		meanpow[i] = np.mean(powlst[i])							#calculate the mean of each bin (sublist)
		stdpow[i] = np.std(powlst[i])							#calculate the standard deviation of each bin

	binwidth = (max(tbin) - min(tbin))/bins						#find the width of the bins
	tbin = tbin[0:-1]+binwidth/2
	
	return tbin, meanpow, stdpow
	
	
#an example set of constants 
"""
#np.random.seed(5)	
a = 100000#											#choose the length of the period (seconds)
n = 10												#choose number of periods recorded 
interval = 100										#choose time between recordings recorded (seconds)
t = np.arange(float(a*n/interval))*interval			#create array of a meassurement times for each period
err = 0.1											#choose error to add to y values

sec = 3600 * 0										#initialize number of passed seconds (3600 converts number of hours to seconds)
stop = 3600 * 8										#choose second to stop recording data (3600 converts number of hours to seconds)
stt = 3600 * 24 *1									#choose second to start recrding	(3600 converts number of hours to seconds)
day = 86400		*1									#number of seconds in a day (this depends on how often sessions of measurement are doesn't necessarily have to be a day)

freq = 1./(20000)									#choose frequency
freq2 = 1./(4000*0.4)								#choose frequency
"""

# an example call of the functions
"""t =np.arange(1000000/100.)*100
t, y = genlc(t, 'binary1')                           #creates a square light curve
T = fold(t, y)                                      #folds it on its period
T, y = norm(T, y)                                   #normalises the time and magnitude
tbin, weights, stdpow = bindata(T, y, bins=10)      #puts data into bins
plt.plot(T, y, 'x'), plt.show()
"""

