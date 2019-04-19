"""BF_ClosestPair: Given a set of n points in the plane, describe by their Cartesian
coordinates, pi(xi, yi). Using the Euclidean distance metric, d(pi, pj) = sqrt[(xi − xj)^2 + (yi − yj)^2],
find the two closest points in the set, using the brute force approach. Write a program
that does the computation and prints the indexes and the coordinates of the closest pair/s.
"""

from math import *

def removeSpaces(xypairs):
    for pair in xypairs:
        if '' in pair:
            del pair[1]

file = open("data6.txt")

n = int(file.readline())
file.readline()

xypairs = []      #create list to hold points

for line in file:               #read points into list of lists
   xypairs.append(line.strip().split(' '))

removeSpaces(xypairs)

CPdist = 1000000
reserve = {}
#do the math
for i in range (0, n):
    for j in range (1, n):
        if i != j:
            dist = sqrt((int(xypairs[i][0]) - int(xypairs[j][0]))**2 + (int(xypairs[i][1]) - int(xypairs[j][1]))**2)
            if dist < CPdist:
                CPdist = dist
                CPindex1 = i
                CPindex2 = j
                c1 = xypairs[i][0] + "," + xypairs[i][1]
                c2 = xypairs[j][0] + "," + xypairs[j][1]
                reserve.clear()
            elif dist == CPdist:
                r1 = xypairs[i][0] + "," + xypairs[i][1]
                r2 = xypairs[j][0] + "," + xypairs[j][1]
                if dist not in reserve:
                    reserve[dist] = [i, j, r1, r2]
                else:
                    reserve[dist].append([i, j, r1, r2])

print("\nClosest Pair Indexes:\n")
print(CPindex1, ", ", CPindex2, ". coordinates: ", c1, ", ", c2, "\n")
if reserve:
    for i in reserve:
        print(reserve[i][0], ", ", reserve[i][1], ". coordinates: ", reserve[i][2], ", ", reserve[i][3], "\n")

file.close()
