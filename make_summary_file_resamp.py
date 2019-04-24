#!/usr/bin/python
#$ -e HE.abc.log
#$ -o HE.abc.log
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-1
#$ -l arch=linux-x64
#$ -l mem_free=2G
#$ -l netapp=1G

import sys
import numpy as np
import math
import random
import gzip
import os
from scipy.stats import beta
from collections import defaultdict

# job id
jid = 1
if 'SGE_TASK_ID' in os.environ:
    jid = int(os.environ['SGE_TASK_ID'])

# params
bw_rho = 0.01
bw_tau = 0.05
bw_p_s = 0.1

# functions
def next_m(m,x,k):
    return m+(x-m)/k

def next_s(s,x,mn,m):
    return s+(x-mn)*(x-m)

def run_var(xs):
    m = xs[0]+0.
    s = 0.
    k = 1
    for x in xs[1:]:
        k += 1
        mn = m
        m = next_m(m,x,k)
        s = next_s(s,x,mn,m)
    return s/k

def get_beta_params(m,s):
    a = (m**2 - m**3 -m*(s**2))/(s**2)
    b = (m-1)*(m**2-m+s**2)/(s**2)
    return ((a,b))

def dist(x0,x1,x2,y0,y1,y2):
    return ((x0-y0)**2+(x1-y1)**2+(x2-y2)**2)**0.5

def get_index(val, bins):
    return int(math.floor(val*len(bins)))
# this is the old version
#def apply_var(data, bins):
#    ret = np.zeros(360)
#    for i in xrange(len(bins)):
#        if i == 0:
#            mylen = bins[0]+0.
#            for j in xrange(0,bins[0]):
#                ret[j] += data[0]/mylen
#        else:
#            mylen = bins[i]-bins[i-1]+0.
#            for j in xrange(bins[i-1],bins[i]):
#               ret[j] += data[i]/mylen   
#    return ret

def apply_var(data, bins):
    ret = np.zeros(360)
    for i in xrange(len(bins)):
        if bins[i] == 0:
            continue
        if i == 0:
            fact = data[i]/np.sum(np.divide(1.,np.arange(1,bins[i]+1)))
            var_in_bins = np.multiply(fact, np.divide(1.,np.arange(1,bins[i]+1)))
            for j in xrange(0,bins[i]):
                ret[j] += var_in_bins[j]
        elif bins[i]-bins[i-1] == 1:
                ret[bins[i]-1] += data[i]
        else:
            if np.sum(np.divide(1.,np.arange(bins[i-1]+1,bins[i]+1))) == 0:
                #print bins[i], bins[i-1], bins
                continue
            fact = data[i]/np.sum(np.divide(1.,np.arange(bins[i-1]+1,bins[i]+1)))
            var_in_bins = np.multiply(fact, np.divide(1.,np.arange(bins[i-1]+1,bins[i]+1)))
            for j in xrange(bins[i-1],bins[i]):
                ret[j] += var_in_bins[j-bins[i-1]-1]
    return ret

# first, get all the data 

pyhash = lambda: defaultdict(pyhash)
#bins_rho = np.arange(0,1,bw_rho)
#bins_tau = np.arange(0,1,bw_tau)

# now adjusting for oversampling of high rho/ high tau
bins_tau = beta.ppf(np.arange(bw_tau,1+bw_tau,bw_tau),0.35,0.35)
bins_rho = beta.ppf(np.arange(bw_rho,1+bw_rho,bw_rho),0.35,0.35)
bins_s_p = np.arange(0,1,bw_p_s)

params = pyhash()

for br in xrange(len(bins_rho)):
    for bt in xrange(len(bins_tau)):
        for bs in xrange(len(bins_s_p)):
            params[br][bt][bs] = []
       
params_file  = sys.argv[1]

f = open(params_file, 'r')
   
for line in f:
    data = line.strip().split()
    data[2:] = [float(x) for x in data[2:]] 
    if beta.cdf(data[2],0.35,0.35) == 1.0:
        data[2] = 0.9999
    if beta.cdf(data[3],0.35,0.35) == 1.0:
        data[3] = 0.9999
    if int(math.floor(beta.cdf(data[2],0.35,0.35)*len(bins_rho))) not in params:
        print  "R", int(math.floor(beta.cdf(data[2],0.35,0.35)*len(bins_rho)))
    elif int(math.floor(beta.cdf(data[3],0.35,0.35)*len(bins_tau))) not in  params[int(math.floor(beta.cdf(data[2],0.35,0.35)*len(bins_rho)))]:
        print "T", int(math.floor(beta.cdf(data[3],0.35,0.35)*len(bins_tau))), beta.cdf(data[3],0.35,0.35)
    elif int(math.floor(data[4]*len(bins_s_p))) not in params[int(math.floor(beta.cdf(data[2],0.35,0.35)*len(bins_rho)))][int(math.floor(beta.cdf(data[3],0.35,0.35)*len(bins_tau)))]:
        print "P", int(math.floor(data[4]*len(bins_s_p)))
    #print int(math.floor(data[2]*len(bins_rho))), int(math.floor(data[3]*len(bins_tau))), int(math.floor(data[4]*len(bins_s_p)))
    #print params[int(math.floor(data[2]*len(bins_rho)))][int(math.floor(data[3]*len(bins_tau)))][int(math.floor(data[4]*len(bins_s_p)))]
    params[int(math.floor(beta.cdf(data[2],0.35,0.35)*len(bins_rho)))][int(math.floor(beta.cdf(data[3],0.35,0.35)*len(bins_tau)))][int(math.floor(data[4]*len(bins_s_p)))].append(data)

#for br in params:
#    for bt in params[br]:
#        for bs in params[br][bt]:
#            if len( params[br][bt][bs]) != 0:
#                print br, bt, bs, params[br][bt][bs]

# now just need to select according to beta dist
ngenes = 10000
nrows = 1000

all_params = []
fh = open(sys.argv[2],'r')
for line in fh:
    data = [float(x) for x in line.strip().split()]
    all_params.append(data)
fh.close()

# open output file
of ="/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/resamp_sims.txt"
ofh = open(of,'w')

nrow = 0

while nrow < nrows:
    # draw parameters for beta distribution of tau
    mean_tau  =all_params[nrow][1] 
    sd_tau =all_params[nrow][4]
    (a_tau,b_tau) = get_beta_params(mean_tau,sd_tau)

    # draw parameters for beta distribution of rho
    mean_rho = all_params[nrow][0]
    sd_rho =all_params[nrow][3]
    (a_rho,b_rho) = get_beta_params(mean_rho,sd_rho)

    mean_p_s = all_params[nrow][2]
    sd_p_s =all_params[nrow][5]
    (a_p_s,b_p_s) = get_beta_params(mean_p_s,sd_p_s)

    
    taus = beta.rvs(a_tau,b_tau, size = ngenes)
    rhos = beta.rvs(a_rho,b_rho, size = ngenes)
    p_ss = beta.rvs(a_p_s,b_p_s, size = ngenes)

    while np.isnan(np.sum(taus)):
        mean_tau  = random.random()
        sd_tau = (((1-mean_tau)*mean_tau)**0.5)*random.random()
        (a_tau,b_tau) = get_beta_params(mean_tau,sd_tau)
        taus = beta.rvs(a_tau,b_tau, size = ngenes)
    while np.isnan(np.sum(rhos)):
        mean_rho = random.random()
        sd_rho = (((1-mean_rho)*mean_rho)**0.5)*random.random()
        (a_rho,b_rho) = get_beta_params(mean_rho,sd_rho)
        rhos = beta.rvs(a_rho,b_rho, size = ngenes)
    while np.isnan(np.sum(p_ss)):
        mean_p_s = random.random()
        sd_p_s = (((1-mean_p_s)*mean_p_s)**0.5)*random.random()
        (a_p_s,b_p_s) = get_beta_params(mean_p_s,sd_p_s)
        p_ss = beta.rvs(a_p_s,b_p_s, size = ngenes)

    rows = [] 

    #all_in = defaultdict(dict)

    for i in xrange(ngenes):
        index_rho = get_index(rhos[i], bins_rho)
        index_tau = get_index(taus[i], bins_tau)
        index_p_s = get_index(p_ss[i], bins_s_p)
    
        if len(params[index_rho][index_tau][index_p_s]) > 0:
            newparam = params[index_rho][index_tau][index_p_s][np.random.choice(xrange(len(params[index_rho][index_tau][index_p_s])))]
            #if newparam[0] not in all_in or newparam[1] not in all_in[newparam[0]]:
            rows.append(newparam)
            #    all_in[newparam[0]][newparam[1]] = 0   
        else:
            while len(params[index_rho][index_tau][index_p_s]) == 0:
                taus[i] = beta.rvs(a_tau,b_tau, size=1)[0]
                rhos[i] = beta.rvs(a_rho,b_rho, size =1)[0]
                p_ss[i] = beta.rvs(a_p_s,b_p_s, size=1)[0]            
    
                while np.isnan(taus[i]):
                    taus[i] = beta.rvs(a_tau,b_tau, size=1)[0]
                while np.isnan(rhos[i]):
                    rhos[i] = beta.rvs(a_rho,b_rho, size =1)[0]
                while np.isnan(p_ss[i]):
                    p_ss[i] = beta.rvs(a_p_s,b_p_s, size=1)[0]            

                index_rho = get_index(rhos[i], bins_rho)
                index_tau = get_index(taus[i], bins_tau)
                index_p_s = get_index(p_ss[i], bins_s_p)
                
                if len(params[index_rho][index_tau][index_p_s]) > 0:
                    newparam = params[index_rho][index_tau][index_p_s][np.random.choice(xrange(len(params[index_rho][index_tau][index_p_s])))]
                    #if newparam[0] not in all_in or newparam[1] not in all_in[newparam[0]]:
                    rows.append(newparam)
    #print mean_rho, mean_tau, mean_p_s
    #for row in rows:
    #    print row
    #exit()
    #print np.mean(taus), np.mean(rhos), np.mean(p_ss), mean_tau, mean_rho, mean_p_s
    #print np.var(taus)**0.5, np.var(rhos)**0.5, np.var(p_ss)**0.5, sd_tau, sd_rho, sd_p_s

    # now, for each row in rows, get the data summary data
    all_vare = []
    tot_var = 0    
    mk = 0.
    sk = 0.

    num_row = 0
    for row in rows:
    
        # open file
        hef = "/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/1000000/K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1/"+row[0]+"/"+row[0]
        hef += "_1000000_PCadj_freqBin_K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1_NJOB="+row[1]+"_he.h2.gz"
        
        try:
            f = gzip.open(hef,'r')
        except:
            continue
        num_row += 1
        for line in f:
            data = line.strip().split()
            mylen = 20
            vare = data[1:(mylen+1)]
            mybins = data[(2*mylen+1):(3*mylen+1)]
            if "NA" in vare:
                while vare[len(vare)-1] == "NA":
                    vare.pop()
                    mybins.pop()
            vare = [float(x) for x in vare]
            mybins = [int(720*float(x)) for x in mybins]
            myvar = apply_var(vare, mybins)
            tot_var+=np.sum(vare)
            all_vare.append(myvar)

    if num_row < ngenes/2.:
        continue
    
    num_all = len(all_vare)    
    all_vare = np.array(all_vare)
    all_vare = all_vare.T
    all_vare_av = np.zeros(360)
    all_vare_std = np.zeros(360)
    for k in xrange(len(all_vare)):
        all_vare_av[k] = np.sum(all_vare[k])/tot_var
        all_vare_std[k] = np.std(np.multiply(all_vare[k],num_all/tot_var))
    
    #print mean_rho, sd_rho, np.sum(all_vare_av)
    #print all_vare_av[0:10]
    #print all_vare_std[0], all_vare_std[1],all_vare_std[2],all_vare_std[5],all_vare_std[20],all_vare_std[60],all_vare_std[120],all_vare_std[240],all_vare_std[359]
    #exit()

    ofh.write(str(mean_rho)+' '+str(mean_tau)+' '+str(mean_p_s)+' '+str(sd_rho)+' '+str(sd_tau)+' '+str(sd_p_s)+' ')
    tot = 0
    cumu_var = np.zeros(360)
    for i in xrange(len(all_vare)):
        cumu_var[i] = tot+all_vare_av[i]
        tot += all_vare_av[i]
    out_freqs = [0,1,2,4,9,19,59,119,179,239,359]
    for fr in xrange(len(out_freqs)):
        if fr == 0:
            #ofh.write(str(cumu_var[out_freqs[fr]])+' '+str(all_vare_std[out_freqs[fr]])+' '+str(all_vare_std[out_freqs[fr+1]])+' ')
            ofh.write(str(cumu_var[out_freqs[fr]])+' ')
        else:
            ofh.write(str(cumu_var[out_freqs[fr]]-cumu_var[out_freqs[fr-1]])+' ')
     
    ofh.write('\n')
    ofh.flush()

    nrow += 1

ofh.close()
