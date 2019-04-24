# 
import sys
import numpy as np

def lsq(arr,arr2):
    diff = np.subtract(arr,arr2)
    tot = np.sum(np.multiply(diff,diff))
    return tot

obs = []

# prior file 
prior = open(sys.argv[1],'r')

all_priors = []
for line in prior:
    data = [float(x) for x in line.strip().split()]
    all_priors.append(data)

# observed data
od = open(sys.argv[2],'r')

count = 0
for line in od:

    obs =  [float(x) for x in line.strip().split()]
    
    comp = {}
    for p in all_priors:
        tot = lsq(p[7:],obs)
        comp[tot] = p[1:7]
 
    i = 0
    of = '/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/est/estimates'+sys.argv[3]+'.'+str(count)+'.post.txt'
    ofh = open(of, 'w')
    for tot in sorted(comp):
        #print >> ofh, tot, 
        for thing in comp[tot]:
            print >> ofh, thing,
        print >> ofh 
        i += 1
        if i >= 1000:
            break

    count += 1
