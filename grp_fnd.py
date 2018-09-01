"""example  of how to define group locations"""

import numpy as np
import matplotlib.pyplot as plt
import random


a = 5
u = np.array([(0,1,0,0,0),(5,1,10,1,0),(1,1,0,1,0),(0,1,1,1,0),(0,15,1,0,0)])	#create example array in which to find groups
print '\n', u
y = np.arange(a)

x = np.arange(a)


n=0													# number of attempts at numbering the nodes
while any(u.reshape(a*a) == 0):						# while any of the nodes are still 0 continue
	n = n+1		
	random.shuffle(x)								# shuffle x so that the nodes are iterated on in a random order									
	random.shuffle(y)								# shuffle y so that the nodes are iterated on in a random order
	for j in y:
		for i in x:	
			gn = 0									# reset the group number
			bound = np.array([0, 0, 0, 0])			# reset the boundry positions
													# array state boundries for ymin, ymax, xmin, xmax in that order
			if  u[j,i] == 1:						# if u in this position is 1 then exit this iteration as it is a boundry
				#print  0
				#print  j, i, bound
				#print  u
				continue
		
			for k in range(a):											# find ymin bound
				if (j-k) < 0:											# if the index is out of the boundry of the matrix
					bound[0] = 0					
					#print  1
					break
				if u[j-k, i] == 1:										# if the label is a bound then label the bound location
					bound[0] = j - k + 1			
					#print  2
					break
				elif (u[j-k, i] != 0) & (u[j-k, i] >= gn):				# if the number is a group number
					gn = u[j-k, i]										# change the group number
					#print  gn*10

			for k in range(a):
				if (j+k) >= a:											# equivilent for ymax bound
					bound[1] = a
					#print  3
					break
				if u[j+k, i] == 1:
					bound[1] = j+k 
					#print  4
					break
				elif (u[j+k, i] != 0) & (u[j-k, i] >= gn):
					gn = u[j+k, i]
					#print  gn*10
				
			for k in range(a):											# equivilent for xmin bound
				if (i-k) < 0:
					bound[2] = 0
					#print  6
					break
				if u[j, i-k] == 1:
					bound[2] = i-k + 1
					#print  7
					break
				elif (u[j, i-k] != 0) & (u[j-k, i] >= gn):
					gn = u[j, i-k]
					#print  gn*10
				
			for k in range(a):											# equivilent for xmax bound
				if (i+k) >= a:
					bound[3] = a
					#print  8
					break
				if u[j, i+k] == 1:
					bound[3] = i+k 
					#print  9
					break
				elif (u[j, i+k] != 0) & (u[j-k, i] >= gn):
					gn = u[j, i+k]
					#print  gn*10
				
			u[bound[0]:bound[1], i] = gn								# change values in ybound to the group number
			u[j, bound[2]:bound[3]] = gn								# change values in xbpund to the group number

	if n > 4:															# break incase zeroes aren't dissappearing
		break
		
print '\n', u, '\n'
print  'fin'









