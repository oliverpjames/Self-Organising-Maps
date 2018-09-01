"""generate light curves"""

import numpy as np
import matplotlib.pyplot as plt 
import  scipy.signal as sg
from astropy.stats import LombScargle 

plt.rcParams.update({'font.size': 25})				#choose font size for figures
#np.random.seed(5)									#uncomment if it is desired to refind exact light curves (number in brackets can be changed)

a = 100000											#choose the length of the period (seconds)
n = 10												#choose number of periods recorded 
interval = 100										#choose time between recordings recorded (seconds)
t = np.arange(float(a*n/interval))*interval			#create array of a meassurement times for each period
err = 0.1											#choose error to add to y values

#create light curves
#b = a/10											#choose feature start point (seconds)
#y = t*0 + 1 										#create an array of length of t where each number is equal to 1

#for i in range(n):
#	y[(a*i + b)/interval : (a*i + 2*b)/interval] = y[(a*i + b)/interval : (a*i + 2*b)/interval]	  -	np.linspace(0,0.02,b/interval)			#give a v shape dip of a binary (2% dip in light)
#	y[(a*i + 2*b)/interval : (a*i + 3*b)/interval] = y[(a*i + 2*b)/interval : (a*i + 3*b)/interval] -	np.linspace(0.02,0,b/interval)		#give a v shape dip of a binary (2% dip in light)

#	y[a*i + b + 60 : a*i + b + 70] = y[a*i + b + 60 : a*i + b + 70] - 	np.linspace(0,0.01,10)												#give a 2nd v shape dip of a binary (1% dip in light)
#	y[a*i + b + 70 : a*i + b + 80] = y[a*i + b + 70 : a*i + b + 80] -	np.linspace(0.01,0,10)												#give a 2nd v shape dip of a binary (1% dip in light)

#	y[a*i + b : a*i + b + 10] = y[a*i + b : a*i + b + 10] 			-	np.linspace(0.02,0.02,10)											#give a u shape dip of a binary (2% dip in light)
#	y[a*i + b + 10 : a*i + b + 20] = y[a*i + b + 10 : a*i + b + 20] - 	np.linspace(0.02,0.02,10)											#give a u shape dip of a binary (2% dip in light)

freq1 = 1./(20000)													#choose frequency
#freq2 = 1./(4000*0.4)												#choose frequency
#y = np.sin( freq1*t*(2*np.pi) ) + np.sin( freq2*t*(2*np.pi) )		#create sinusoids

y = sg.sawtooth(2*np.pi*freq1*t)				#sawtooth function
y = sg.sawtooth(2*np.pi*freq1*t, 0)				#reverse sawtooth function
y = sg.sawtooth(2*np.pi*freq1*t, 0.5)			#triangle function
y = sg.square(2*np.pi*freq1*t)					#square wave
#plt.step(t,y), plt.show()



#delete data to create a windoww function
dellist = []															#create empty list
sec = 3600 * 0															#initialize number of passed seconds (3600 converts number of hours to seconds)
stop = 3600 * 8															#choose second to stop recording data (3600 converts number of hours to seconds)
stt = 3600 * 24 *1														#choose second to start recrding	(3600 converts number of hours to seconds)
day = 86400		*1														#number of seconds in a day (this depends on how often sessions of measurement are doesn't necessarily have to be a day)

for i in range(len(t)*interval/day+1):									#for each set of 24 hours
	dellist.append(range((stop + sec)/interval, (stt + sec)/interval))	#add a list of hours to be deleted
	sec = sec + day														#adds time (eg a day) to the total time passed 
dellist = np.reshape(dellist, len(dellist)*(stt-stop)/interval)			#reshape dellist to be a list of integers rather than a list of lists
  
t = np.delete(t, dellist)												#delete time data
y = np.delete(y, dellist)												#delete flux data

	
#add normally distributed errors to y
for i in range(len(y)):
	y[i] = y[i] + np.random.normal(0, err)


frequency, power = LombScargle(t, y).autopower()				#find frequencies that contribute to describe data 
f1 = frequency[np.argmax(power)] 								#find the most important frequency, this will likely be the period of the orbit

power[np.argmax(power)*1-250:np.argmax(power)*1+250] = -0.1		#lower the power of the data around the highest frequency
f2 = frequency[np.argmax(power)]								#find second most important frequency

#fbeat = abs(f1-f2)												#find beat frequency
fpat = f1


#create array of folded times T
T = t*0														#initialize array for folded times T 
c = 1														#initialize c the multiplication factor to choose the bin in which each datapoint goes
clst = []
tlst = []
for i in range (len(t)):									#iterate over all numbers in t
	while (t[i] >= c/fpat):									#if the time in t[i] is larger than the bin time
		c = c + 1											#	move to the next bin
	T[i] = t[i] - (c-1)/fpat								#remove c periods from the time 

T = T/max(T)												#normailse length of T
ynorm = y-min(y)											#normalise power
ynorm = ynorm/max(ynorm)

bins = 10													#choose the number of bins
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
	powlst[binindex[i]].append(ynorm[i])					#append binindex[i] list of powers with the values of y that should be in that bin


meanpow = np.zeros(bins+1)									#create an array of the same length of tbins in which to hold the mean of the power in those bins
stdpow = np.zeros(bins+1)

for i in range(bins):
	meanpow[i] = np.mean(powlst[i])							#calculate the mean of each bin (sublist)
	stdpow[i] = np.std(powlst[i])							#calculate the standard deviation of each bin

meanpow[-1] = meanpow[-2]											
stdpow[-1] = stdpow[-2]


binwidth = (max(tbin) - min(tbin))/bins						#find the width of the bins


#plt.figure(1, figsize=(30,5)), plt.plot(t,y, 'x'), plt.xlabel('Time (s)'), plt.ylabel('Luminosity')
#plt.savefig('data_orig.png', bbox_inches='tight')
#plt.show()											#plot original data

#plt.figure(1, figsize=(30,5)), plt.plot(frequency[0:len(frequency)/25], power[0:len(frequency)/25]),  plt.xlabel('Frequency (Hz)'), plt.ylabel('Power')
#plt.savefig('periodogram.png', bbox_inches='tight')	
#plt.show()															#plot periodogram
#plt.plot([abs(freq1-freq2)/7], [max(power)], 'go'), 
#plt.plot([f1, f2], [max(power), max(power)], 'rx'), 
#plt.plot(1/day, max(power), 'cx'), 
#plt.plot(fbeat, max(power)-0.01, 'ko'), 
#plt.plot(t, env, 'gx')
#plt.plot(t, np.sin(2*np.pi*fbeat*t), 'r')
#plt.plot(t, np.sin(2*np.pi*f2*t), 'k')


plt.figure(1, figsize=(30,5)), plt.plot(T,ynorm, 'x')#, plt.show()	
#plt.plot(t*8/3, sin(t, 100))
#plt.plot(t, np.sin(2*np.pi*100*fbeat*t))
#plt.plot(x, a, 'g')#, '.')
#plt.plot(x, -b, 'r')#, '.')
#plt.plot(x, d, '--r')
#plt.plot(x, e, '--g')
#plt.grid()																						#plot folded data   

plt.step(tbin, meanpow, 'r', where='post',linewidth=3.0),  plt.xlabel('Fraction of Period'), plt.ylabel('Luminosity')
	#plot binned data

plt.errorbar(tbin[0:-1]+binwidth/2, meanpow[0:-1], stdpow[0:-1], linestyle='None', ecolor='r', linewidth=2.0, capsize=4)#, plt.show()
plt.savefig('data_folded.png', bbox_inches='tight')		
plt.show() 	

