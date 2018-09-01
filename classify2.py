"""create a U-matrix, and boundary map"""

import numpy as np
import matplotlib.pyplot as plt
import genlcfunc as func            # load in the functions to create and fold light curves
N = np.load("lcmapnorm.npy")		# load in a SOM

#define size of map base on N
a = len(N[:,0,0])
b = len(N[0,:,0])
p = len(N[0,0,:])
Ua = 2*a-1
Ub = 2*b-1


#create blanck U-matrix
U = np.zeros([2*a-1, 2*b-1])	

#create U-matrix
for i in range (a-1):
	for j in range (b-1):
	
		U[2*i+1, 2*j] = np.sqrt( np.sum( np.square(N[i,j,] - N[i+1,j,] )))			                                                           # Find column
		U[2*i, 2*j+1] = np.sqrt( np.sum( np.square(N[i,j,] - N[i,j+1,] )))			                                                           # Find rows
		U[2*i+1, 2*j+1] = np.mean([np.sqrt( np.sum( np.square(N[i,j,] - N[i+1,j+1,] ))),np.sqrt( np.sum( np.square(N[i,j+1,] - N[i+1,j,] )))]) # diagonal distance

for i in range (a-1):
	U[2*a-2, 2*i+1] = np.sqrt( np.sum( np.square(N[a-1,i+1,] - N[a-1,i,] )))		# fill last row

for j in range (b-1):
	U[2*j+1, 2*b-2] = np.sqrt( np.sum( np.square(N[j+1,b-1,] - N[j,b-1,] )))		# fill last column

"""# fill SOM squares
for i in range (0, a-1):
	for j in range (1, b-1):
		U[2*i, 2*j] = np.mean([U[2*i-1, 2*j-1], U[2*i-1, 2*j], U[2*i, 2*j+1], U[2*i, 2*j-1], U[2*i, 2*j+1], U[2*i+1, 2*j-1], U[2*i+1, 2*j], U[2*i+1, 2*j+1]])
"""

# Create boundary map

# test if local max with comparable cells
B = U*0
for i in range (6,Ua-6): 				#test the inner cells
	for j in range (6,Ub-6):
		# if U[i,j] is the largest of the values in the vertical of i+-2,+-4, and +-6 then add it ot the bounds matrix
		if (U[i-2, j] < U[i, j]) & (U[i+2, j] < U[i, j]) & (U[i-4, j] < U[i, j]) & (U[i+4, j] < U[i, j]) & (U[i-6, j] < U[i, j]) & (U[i+6, j] < U[i, j]):
			B[i, j] = U[i, j]
		# if U[i,j] is the largest of the values in the horizontal of i+-2,+-4, and +-6 then add it ot the bounds matrix
		if (U[i, j-2] < U[i, j]) & (U[i, j+2] < U[i, j]) & (U[i, j-4] < U[i, j]) & (U[i, j+4] < U[i, j]) & (U[i, j-6] < U[i, j]) & (U[i, j+6] < U[i, j]) :
			B[i, j] = U[i, j]
			
# testing local max at edges
for i in range(6,Ua-6):  				# test upper rows
	for j in range (6):
		# only test horisontal 
		if (U[i-2, j] < U[i, j]) & (U[i+2, j] < U[i, j]) & (U[i-4, j] < U[i, j]) & (U[i+4, j] < U[i, j]) & (U[i-6, j] < U[i, j]) & (U[i+6, j] < U[i, j]):
			B[i, j] = U[i, j]		

for i in range(6):						#test left columns
	for j in range (6, Ub-6):
		# only test vertical 
		if  (U[i, j-2] < U[i, j]) & (U[i, j+2] < U[i, j])& (U[i, j-4] < U[i, j]) & (U[i, j+4] < U[i, j]) & (U[i, j-6] < U[i, j]) & (U[i, j+6] < U[i, j]) :
			B[i, j] = U[i, j]
			
for i in range(6,Ua-6):			        # test bottom rows
	for j in range (Ub-6,Ub): 
		# only test horisontal 
		if (U[i-2, j] < U[i, j]) & (U[i+2, j] < U[i, j]) & (U[i-4, j] < U[i, j]) & (U[i+4, j] < U[i, j]) & (U[i-6, j] < U[i, j]) & (U[i+6, j] < U[i, j]):
			B[i, j] = U[i, j]		

for i in range(Ua-6,Ub):				# test right columns
	for j in range (6, Ub-6):
		# only test vertical 
		if  (U[i, j-2] < U[i, j]) & (U[i, j+2] < U[i, j])& (U[i, j-4] < U[i, j]) & (U[i, j+4] < U[i, j]) & (U[i, j-6] < U[i, j]) & (U[i, j+6] < U[i, j]) :
			B[i, j] = U[i, j]


		
# boolean array that is true if there is no point and false if there is one (wrong way round but helpful for next bit)
npoint = np.zeros((Ua,Ub), dtype=bool)	
for i in range (Ua):
	for j in range (Ub):
		if B[i,j] != 0:
			npoint[i,j] = False
			
		else:
			npoint[i,j] = True
			
# remove all datapoints that are on their own (to clean data)
for i in range (1,Ua-1):
	for j in range (1,Ub-1):
		if npoint[i,j] == False:
			
			if all( np.delete( np.reshape(npoint[i-1:i+2,j-1:j+2], (9,1) ), 4 ) ):	
				B[i,j] = 0		# set bounnds matrix = 0

# set points on boundaries = 1
B = B/B			#create array of 1's where there are values and nan's elsewhere
for i in range(Ua):
	for j in range (Ub):
		# set all nan's = 0
		if np.isnan(B[i,j]) == True:
			B[i,j] = 0

"""
#make the locations of classified light curves clear to see in the map
U[2*21,2*  8] = 1.75
U[2*20,2*  8] = 1.75
U[2*20,2*  9] = 1.75
U[2*20,2* 10] = 1.75
U[2*19,2* 11] = 1.75
U[2* 5,2*  8] = 1.75
U[2* 0,2* 25] = 1.75
U[2*20,2* 25] = 1.75

U[2*14,2* 23] = 1.00
U[2*31,  0] = 1.00
U[2*34,2* 13] = 1.00
U[ 0, 2* 0] = 1.00
U[2* 8,2* 38] = 1.00
"""

# Draw
plt.figure(figsize=(10,10))	                                    #choose size of Figure
plt.imshow(B, cmap='gist_earth', interpolation='nearest')       #draw plot if B then boundry map will be plotted, if U then U-matrix will be plotted
                                                                #cmap=gist_earth for colour plot, cmap=gray for greyscal plot
#plt.colorbar()#scale=0.8)
#plt.savefig('U-matrix test.png')
plt.show()

print 'fin classify'



