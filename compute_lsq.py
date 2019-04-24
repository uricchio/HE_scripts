# 
import sys
import numpy as np

def lsq(arr,arr2):
    diff = np.subtract(arr,arr2)
    tot = np.sum(np.multiply(diff,diff))
    return tot

fh = open(sys.argv[1], 'r')

obs = []
for line in fh:
    obs =  [float(x) for x in line.strip().split()]

fh2 = open(sys.argv[2],'r')
    
comp = {}
for line in fh2:
    data =  [float(x) for x in line.strip().split()]
    tot = lsq(data[6:],obs)
    comp[tot] = data
 
i = 0
for tot in sorted(comp):
    #print tot, 
    for thing in comp[tot][0:6]:
        print thing,
    print
    i += 1
    if i >= 5000:
        break

