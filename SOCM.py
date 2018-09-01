"""create self-organising map of coulours"""

import numpy as np
import matplotlib.pyplot as plt
import random
import matplotlib.animation as animation


def Nnew(N, c, r0, t, d, L0, a, b):
    """function to adjust the weights of each node"""
    #adjust the weights of each element in the neighbourhood of the BMU
    N1 = N															#create a coppy of array N

    delta = N[:,:,] - c 											# array of differences between chosen colour and each node weight
    dist = np.zeros((a, b)) 										# array of zeros of dimensions of the grid

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
		
		    if np.sqrt( (float(i)-(BMU1[0,0]))**2 + (float(j)-(BMU1[0,1]))**2 ) <= neibr1:      #if node is in neibourhood
			    N1[i,j,] = N[i,j,] + adjweight1[i,j]*learnrate1* (c - N[i,j,])					#adjust the weight of the node, according to corresponding weights

    return N1


# define variables
t = int(0)							        #initialise timestep
t_max = int(500) 					        #number of iterations							(Choose)
a = int(40) 						        #width of grid									(Choose)
b = int(40) 						        #height of grid									(Choose)
L0 = float(0.1) 					        #initial learning rate							(Choose)
r0 = max(float(a),float(b))*0.5             #initial neighbourhood radius					(Choose)
d = t_max/np.log(r0)				        #decay constant

cols = np.array([[1,1,1],[1,1,0],[1,0,1],[0,1,1],[0,0,1],[0,1,0],[1,0,0]])		#	(Choose) colours to be organised (here primary, secondary, and white light)

#np.random.seed(5)				#set the seed for the random numbers, to keep them the map the same each time SOCM is run, helpful in testing
N = np.random.rand(a, b, 3)		# create grid of random RGB values

#create a list of figures of the map organising itself
ims = []						                                #create an empty list of images
fig = plt.figure(figsize=(14,14))		                        #create a figure on which to plot sequence of images

#im = plt.imshow(N, interpolation='nearest', animated=True)	    #create an image of RGB values, viewed as 'closest' colours that can be displayed
#ims.append([im])											    #add im to the list of images
                                                                #uncomment if gif of organisation process is desired instead ofstill image of completed map

for t in range(t_max):                                          #iterate from 0 to t_max
    for i in range( len(cols[:,])):                             #show all colours at each iteration
	    N = Nnew(N, cols[i], r0, t, d, L0, a, b)				#find new node weights

im = plt.imshow(N, interpolation='nearest', animated=True)	    #create an image of closest RGB values
ims.append([im])											    #add im to the list of images
                                                                #if gif of the organisation process is desired then these lines should be indented to fit into one or both of the for loops
                                                                #if still image of the resulting plot is desired then lave unindented

anim = animation.ArtistAnimation(fig, ims, interval=1)			#put images in ims into a sequence
#anim.save('Rules test L0=test.png')			            	#save animation 
plt.show()	                                                    #show animation

print 'fin SOCM'			                                  


