"""create SOMs and classify light curves"""

import numpy as np
import matplotlib.pyplot as plt
import random
import matplotlib.animation as animation
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import genlcfunc as func
from time import time as tim

def Nnew(N, c, r0, t, d, L0, a, b):
	"""adjust the weights of each node"""
	#adjust the weights of each element in the neighbourhood of the BMU
	N1 = N							#create a coppy of array N
	
	delta = N[:,:,] - c 			# array of differences between chosen weight and each node weight
	dist = np.zeros((a, b)) 		# array of zeros of dimensions of the grid
	
	for i in range(a):
		for j in range(b):
			dist[i,j] = np.sqrt( np.sum( np.square(delta[i,j,] )))	#find Euclidean distance between input RGB value and node RGB value)
	#
								
	BMU 		= np.argwhere( dist == np.min(dist) ) 				#find position of best matching unit returned as 2D array
	
	neibr 		= r0*np.exp(-t/d)									#find neibourhood radius at particular timestep
				
	learnrate 	= L0*np.exp(-t/d)									#find learning rate at particular timestep
					
	adjweight 	= np.exp(-np.square(dist)/(2*neibr**2))				#find weight adjustment with distance
	
	#iterate through each node
	for i in range(a):
		for j in range(b):
			
			if np.sqrt( (float(i)-(BMU[0,0]))**2 + (float(j)-(BMU[0,1]))**2 ) <= neibr: 		# if node is in neibourhood
				N1[i,j,] = N[i,j,] + adjweight[i,j]*learnrate* (c - N[i,j,])					# adjust the weight of the node, according to corresponding weights

	return N1

def BMU(weights, N=N, a=a, b=b):
	"""function that finds the location of the BMU of a light curve"""
	dist = np.zeros((a, b)) 												# array of zeros of dimensions of the grid

	for i in range(a):
		for j in range(b):
			dist[i,j] = np.sqrt( np.sum( np.square(N[i,j,] - weights )))	#find Euclidean distance between input RGB value and node RGB value)
			
	closest = np.argwhere( dist == np.min(dist) )							#find nearest neibour
	return np.array([closest[0,0], closest[0,1]])
	
#change line widths
mpl.rcParams['lines.linewidth'] = 1

# define variables
t = int(0)							#initialise timestep
t_max = int(90) 					#number of iterations							(Choose)
a = int(40) 						#width of grid									(Choose)
b = int(40) 						#height of grid									(Choose)
L0 = float(0.2) 					#initial learning rate							(Choose)
r0 = max(float(a),float(b))/1		#initial neighbourhood radius					(Choose)
d = t_max/np.log(r0)				#decay constant
# reserve i,j,k for indexing arrays

np.random.seed(5)

#weights = np.array([[1,1,1,1,1],[1,0,0,1,1],[1,1,0.5,1,1],[0,0,1,1,1],[1,0,1,0,1],[1,1,0,1,1],[0,1,0,1,0],[0,0,0,0,0]])		#	(Choose) weights to be organised (here primary, secondary, and white light)	
#weights = np.array([[1,1,1,1,1],[0,0,1,1,1],[1,1,1,0,0],[1,0,1,0.5,1],[1,0.8,0.6,0.4,0.2]])	
	
p = 10								#number of points on curve (change when changing weights)

#create training data from genlcfunc
time = np.arange(1000000/100.)*100

weights = np.zeros((5, 10))
err=0.1

time, y = func.genlc(time, 'sin', err=err)
T = func.fold(time, y)
tbin, weights[0], stdpow = func.bindata(T, y, bins=10) 
tbin, weights[0] = func.norm(tbin, weights[0])			

time, y = func.genlc(time, 'rsawtooth', err=err)
T = func.fold(time, y)
tbin, weights[1], stdpow = func.bindata(T, y, bins=10) 
tbin, weights[1] = func.norm(tbin, weights[1])

time, y = func.genlc(time, 'triangle', err=err)
T = func.fold(time, y)
tbin, weights[2], stdpow = func.bindata(T, y, bins=10) 
tbin, weights[2] = func.norm(tbin, weights[2])
weights[2] = -weights[2]

time, y = func.genlc(time, 'square', err=err)
T = func.fold(time, y)
tbin, weights[3], stdpow = func.bindata(T, y, bins=10) 
tbin, weights[3] = func.norm(tbin, weights[3])

time, y = func.genlc(time, 'triangle', err=err)
T = func.fold(time, y)
tbin, weights[4], stdpow = func.bindata(T, y, bins=10) 
tbin, weights[4] = func.norm(tbin, weights[4])
	
#np.random.seed(100)							
N = np.random.rand(a, b, p)					#create grid of random values final number must be the same length as input values 



#iterate T_max times for each weight
for t in range(t_max):
	for i in range( len(weights[:,])):
		N = Nnew(N, weights[i], r0, t, d, L0, a, b)										#find new node weights at each step
		
for i in range (a):
	for j in range (b):
		N[i,j,:] = ( N[i,j,:] - min(N[i,j,:]) )/( max(N[i,j,:]) - min(N[i,j,:]) ) 		#normailse the nodes
		


	



#create light curves from combinations of group centres
#0: sin, 1: rsawtooth, 2: upside down triangle, 3: Square, 4: triangle

boundweights = np.zeros((8, 10))

tbin, boundweights[0] = func.norm(tbin, weights[1] +  0.5*weights[0])
tbin, boundweights[1] = func.norm(tbin, weights[1] + 0.75*weights[0])
tbin, boundweights[2] = func.norm(tbin, weights[1] +      weights[0])
tbin, boundweights[3] = func.norm(tbin, weights[1]*0.75 + weights[0])
tbin, boundweights[4] = func.norm(tbin, weights[1]*0.5  + weights[0])

tbin, boundweights[5] = func.norm(tbin, weights[0] + weights[3])
tbin, boundweights[6] = func.norm(tbin, weights[0] + weights[4])
tbin, boundweights[7] = func.norm(tbin, weights[0] + weights[2])

#create array of BMU positions of combined group centres
bound = np.array([ BMU(boundweights[0]), BMU(boundweights[1]), BMU(boundweights[2]), BMU(boundweights[3]),	BMU(boundweights[4]),
					 BMU(boundweights[5]),	BMU(boundweights[6]), BMU(boundweights[7]) ])	
					 
#create array of group centre position
neibcent = np.array([ BMU(weights[0]), BMU(weights[1]), BMU(-weights[4]), BMU(weights[3]),	BMU(weights[4]) ])		s

print '\n', bound
print '\n', neibcent


plt.figure(figsize=(18,18))			     		  		#create a figure
plt.subplots_adjust(wspace=0, hspace=0)	     		  	#remove whitespace between figures

axes = []								        		#create empty list in which the figure data can be put
x=np.arange(0,p)										#set x-values must be the same length as input 

n=1										       			#start the figure number at 1

for i in range (a):
	for j in range(b):
		axes.append(plt.subplot(b,a,n))		            					#append subplot number n on a grid of dimensions a by b
		
		for k in range(8):													
			if (i == bound[k,0]) & (j == bound[k,1]):						#plot locations of combined light curves as a thicker red line
				mpl.rcParams['lines.linewidth'] = 2
				plt.step(x, N[i,j,:], 'r', where='mid')
				mpl.rcParams['lines.linewidth'] = 1
				break
			
			elif (k < 5):
				if (i == neibcent[k,0]) & (j == neibcent[k,1]):				#plot locations of neighbourhood centres as a thicker green line
					mpl.rcParams['lines.linewidth'] = 2
					plt.step(x, N[i,j,:] , 'g', where='mid') 
					mpl.rcParams['lines.linewidth'] = 1
					break
			
			elif k == 5:
				plt.step(x, N[i,j,:], where='mid')							#draw a step line on the subplot (if not classifying data this is the only line in the for loop that is needed and the loop itself can be commented out)

		plt.ylim([-0.5,1.5])				           						#set y-limits to keep the graphs the same shape
		plt.xlim([-0.5,p-0.5])			            						#as above remember to change when changing number of points in weights 
		n = n + 1							            					#increase figure number


plt.setp(axes, xticks=[], yticks=[])					#remove ticks and numbers from the axes of the subplots in list axes


#plt.savefig('SOlcM test.png')							#save SOM
#plt.show()												#show SOM
#np.save("lcmapnorm", N)								#save the SOM array so it can be used in other programs

print 'fin SOlcM'
