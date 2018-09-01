"""Time self-organising map of coulours, and write data to file on the qulity of the SOM"""

import numpy as np
import matplotlib.pyplot as plt
import random
import matplotlib.animation as animation
import time #import time as tim

def Nnew(N, c, r0, t, d, L0, a, b):
    """function that adjusts the weights of each node"""
    #adjust the weights of each element in the neighbourhood of the BMU
    N1 = N							                                #create a coppy of array N

    delta = N[:,:,] - c 			                                # array of differences between chosen colour and each node weight
    dist = np.zeros((a, b)) 		                                # array of zeros of dimensions of the grid

    for i in range(a):
	    for j in range(b):
		    dist[i,j] = np.sqrt( np.sum( np.square(delta[i,j,] )))	#find Euclidean distance between input RGB value and node RGB value)
							
    BMU1 		= np.argwhere( dist == np.min(dist) ) 				#find position of best matching unit returned as 2D array

    neibr1 		= r0*np.exp(-t/d)									#find neibourhood radius at particular timestep
			
    learnrate1 	= L0*np.exp(-t/d)									#find learning rate at particular timestep
				
    adjweight1 	= np.exp(-np.square(dist)/(2*neibr1**2))			#find weight adjustment with distance

    #iterate through each node
    for i in range(a):
	    for j in range(b):
		
		    if np.sqrt( (float(i)-(BMU1[0,0]))**2 + (float(j)-(BMU1[0,1]))**2 ) <= neibr1: 		#if node is in neibourhood
			    N1[i,j,] = N[i,j,] + adjweight1[i,j]*learnrate1* (c - N[i,j,])					#adjust the weight of the node, according to corresponding weights

    return N1


# define variables
t = int(0)									#initialise timestep
t_max = int(500) 							#number of iterations							(Choose)
a = int(40) 								#width of grid									(Choose)
b = int(40) 								#height of grid									(Choose)
L0 = float(0.1) 							#initial learning rate							(Choose)
r0 = max(float(a),float(b))*0.5     		#initial neighbourhood radius					(Choose)
d = t_max/np.log(r0)						#decay constant

cols = np.array([[1,1,1],[1,1,0],[1,0,1],[0,1,1],[0,0,1],[0,1,0],[1,0,0]])		#	(Choose) colours to be organised (here primary, secondary, and white light)


t0 = np.zeros(10)
t1 = np.zeros(10)
for num in range (10):
    np.random.seed(5)				                                #set the seed for the random numbers, to keep them the map the same each time SOCM is run, increases reproducability

    t0[num] = time.time()  
    N = np.random.rand(a, b, 3)		                                #create grid of random RGB values

	
    for t in range(t_max):
	    for i in range( len(cols[:,])):
		    N = Nnew(N, cols[i], r0, t, d, L0, a, b)				#find new node weights

    t1[num] = time.time()

ttot = t1-t0                                                        #create array of times to complete organisating process


fh = open("rules test L0.txt", "a+")                                                                                                        #open a file for appending, if one doesn't exist then create and open one

fh.write("%f\t\t" %L0 + "%f\t\t" %np.mean(ttot) + "%f\t\t" %np.std(ttot))                                                                   #append data about the learning rate mean time and standard deviation of time

for k in range(7):                                                                                                                          #for each of the colours
    dist = np.zeros((a, b)) 															                                                    #create array of zeros of dimensions of the grid

    for i in range(a):
	    for j in range(b):
		    dist[i,j] = np.sqrt( np.sum( np.square(N[i,j,] - cols[k] )))				                                                    #find Euclidean distance between input RGB value and all node RGB values
			
    closest = np.argwhere( dist == np.min(dist) )										                                                    #find location of BMU
    nneib = N[closest[0,0], closest[0,1]]												                                                    #find weights of BMU
    fh.write("%f\t\t" %cols[k,0] +"%f\t\t" %cols[k,1] +"%f\t\t" %cols[k,2] + "%f\t\t" %nneib[0] + "%f\t\t" %nneib[1] + "%f\t\t" %nneib[2] ) #append the colour, and BMU's weights to the file

fh.write('\n')                                                                                                                              #return line in file
fh.close()                                                                                                                                  #close file

print 'finished SOCM timing test for L0 =', L0			                                  


