#!/usr/bin/python
#$ -e sim.abc.log
#$ -o sim.abc.log
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-100
#$ -l arch=linux-x64
#$ -l mem_free=4G
#$ -l netapp=1G

import sys
import os
import subprocess
import math
import numpy as np
import re
import gzip
import shutil
from collections import defaultdict
from scipy import interpolate
from scipy.stats import gaussian_kde
from scipy.optimize import minimize

#from scipy.interpolate import interp1d

num_params = 6
fac = 1000
id = 1
if 'SGE_TASK_ID' in os.environ:
    id = int(os.environ['SGE_TASK_ID'])
mnum = fac*int(id)-(fac-1)

# this is the master file with all the params and summary stats
all_priorfile = sys.argv[1]

# typically 0.0005 and 25
tol = sys.argv[3]
spaces = sys.argv[4]

fid = str(id)+'.t'+tol+'.s'+spaces

# this the path prefix to the out files
outfile = sys.argv[2]+'.'+fid+'.txt'
f = open(all_priorfile,'r')
out = open(outfile, 'w')

# get prior data
data = {}
te = re.compile('nan')
cur = 0
for line in f:
    if(te.search(line)):
        #print line.strip()
        continue
    line = line.strip().split()

    #for i in xrange(num_params+7,len(line)-1):
    #    line[i] = 1-newline[i-num_params]
    #if float(line[4]) < 0.6:
    #    continue
    #a[1./nsam,2./nsam,3./nsam,4./nsam,5./nsam,7./nsam,10./nsam,15./nsam,20./nsam,50./nsam,100./nsam,150./nsam,200./nsam,250./nsam,300./nsam,1.]
    #data[cur] = [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8]+line[9]+line[10]+line[11]+line[12],line[13]+line[14]+line[15]+line[16],line[17]+line[18]+line[19]+line[20],line[21]]
    data[cur] = line
    cur += 1
f.close()

# 

# infer peak from hist
def infer(post,wid,mval):
    if len(post) == 0:
        return float('nan')
 
    # comment these 3 lines to use old histogram based estimation
    dens = gaussian_kde(post) #bw_method=0.01)
    est = minimize(lambda x: -1.*dens(x), np.mean(post)) 
    if est.x[0] > 1:
        return 1
    if est.x[0] < 0:
        return 0
    return est.x[0]

    post = sorted(post)
    #temp = np.histogram(post,bins=np.arange(0,mval,wid))
    #post.insert(0,0.)
    tot = range(1,len(post)+1)
    xvals = np.arange(0,mval,wid)
    
    try:
        myinter = interpolate.interp1d(post,tot,bounds_error=False)
    except:
        return float('nan')

    try:
        interpout = myinter(xvals)
    except:
        return float('nan')
    
    for i in xrange(len(interpout)):
        if math.isnan(interpout[i]):
            interpout[i] = 0.
    
    check = np.gradient(interpout)
    #myinter2 = interpolate.UnivariateSpline(xvals,check)    
    #check = myinter2(xvals)

    k = 0
    for thing in check:
        if thing == max(check):
            break 
        k+=1
  
    # figure out how to do gaussian kde instead to fix artifacts

    #dens = gaussian_kde(post) #bw_method=0.01)
    #est = minimize(lambda x: -1.*dens(x), mean(post)) 

    return xvals[k]+wid/2.
  

lnum = 0
i = 0

data2 = {}
pri = {}


for cur in sorted(data):

    if not (cur >= (mnum-1) and cur < (mnum + fac -1)):
        pri[cur] = data[cur]
    else:
        data2[cur] = data[cur]

data = {}
priorfile = '/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/priors/out.prior.'+fid+'.txt'
datafile = '/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/data/out.data.'+fid+'.txt'

p = open(priorfile, 'w')
d = open(datafile, 'w')

for m in sorted(pri):

    print >> p, m,
    for suma in pri[m]:
        print >> p, suma,
    print >> p
 
p.close()

tr = []
for cur in sorted(data2):

    tr.append(data2[cur][0:num_params])
    for suma in data2[cur][num_params:]:
        print >> d, suma,
    print >> d

d.close()

exec_array = ['python','/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/compute_lsq_sim.py', priorfile, datafile, fid]
for thing in exec_array:
    print thing,
print

p1 = subprocess.Popen(exec_array)
p1.wait()

for m in sorted(data2):

    # open posterior file
    t = open('/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/est/estimates'+fid+'.'+str(lnum)+'.post.txt','r')
    a = [[] for i in xrange(num_params)]
    for line in t:
        line = line.strip().split()
        place = 0
        for thing in line:
            a[place].append(float(line[place]))
            place +=1
    t.close()
    
    rhotrue = float(tr[lnum][0])
    tautrue = float(tr[lnum][1])
    sptrue = float(tr[lnum][2])
    rhosdtrue = float(tr[lnum][3])
    tausdtrue = float(tr[lnum][4])
    spsdtrue = float(tr[lnum][5])
    
    #Ltrue = float(tr[lnum][4])
    rho = infer(a[0],1./(float(spaces)),1.)     
    #if rhotrue > 0.95 and tautrue > 0.95:
    #    print rhotrue
    #    print  
    #    for thing in rho:
    #        print thing
    #    print 
    #    for thing in a[0]:
    #        print thing
    #    exit()
    tau = infer(a[1],1./(float(spaces)),1.)     
    sp = infer(a[2],1./(float(spaces)),1.)     
    rhosd = infer(a[3],0.5/(float(spaces)),0.5)     
    tausd = infer(a[4],0.5/(float(spaces)),0.5)     
    spsd = infer(a[5],0.5/(float(spaces)),0.5)     
    #L = infer(a[4],0.05,1)

    #print >> out, mtrue,sdvtrue,tautrue,rhotrue,Ltrue,me,sdv,tau,rho,L,np.mean(a[0]),np.mean(a[1]),np.mean(a[2]),np.mean(a[3]),np.mean(a[4])
    print >> out, rhotrue,tautrue, sptrue,rhosdtrue,tausdtrue,spsdtrue,rho,tau,sp,rhosd,tausd,spsd,np.mean(a[0]),np.mean(a[1]),np.mean(a[2]),np.mean(a[3]),np.mean(a[4]),np.mean(a[5])
     
    f_in = open('/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/est/estimates'+fid+'.'+str(lnum)+'.post.txt','rb')
    f_out = gzip.open('/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/est/estimates'+fid+'.'+str(lnum)+'.post.txt.gz','wb')

    shutil.copyfileobj(f_in, f_out)
    f_in.close()
    os.remove('/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/est/estimates'+fid+'.'+str(lnum)+'.post.txt')
    
    lnum += 1

out.close()
os.remove(priorfile)
