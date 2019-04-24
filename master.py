#!/usr/bin/python
#$ -e HE.abc.log
#$ -o HE.abc.log
#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-10560
#$ -l arch=linux-x64
#$ -l mem_free=1G
#$ -l netapp=1G

import gzip
import sys
import os
import numpy as np
import random
import subprocess
from scipy.stats import beta

# job id
jid = 1
if 'SGE_TASK_ID' in os.environ:
    jid = int(os.environ['SGE_TASK_ID'])

# params
linemax = 10**5

# functions 
def maf(ac,n):
    return np.min([2*n-ac,ac])

# first, get a gene
gene_num = jid
gg = open('/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/goodExp.txt','r')

i = 0
gene_name = ''
for line in gg:
    if i == (gene_num-1):
        gene_name = line.strip().split()[0] 
        break
    i += 1    

#gene_file_list = sys.argv[1]
#gl = open(gene_file_list,'r')

#gene_list = []
#for line in gl:
#    line = line.strip()
#    if line[-6:] == "012.gz":
#        gene_list.append(line.strip())
#gl.close()

mygene = '/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/1000000/'+gene_name+'.1000000.txt.012.gz'
parfile= "/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/params/"+gene_name+".txt"
parf = open(parfile, 'w')
# get actual name instead of file name
#gene_name = mygene.split("/")
#gene_name = gene_name[1].split(".")[0]

for sim in xrange(1,201):
    try:
        gl = gzip.open(mygene,'r')
    except:
        print "no such gene in data"
        exit(-1)

    genos = []
    freqs = []
    for line in gl:
        data = [int(x) for x in line.strip().split()]
        genos.append(data)
        freqs = np.zeros(len(data))
        freqs = np.add(freqs,data)
        break

    for line in gl:
        data = [int(x) for x in line.strip().split()]
        genos.append(data)
        freqs = np.add(freqs,data)
    gl.close()

    # now make a list of lists with all the positions of the freqs
    # first get mafs
    for i in xrange(len(freqs)):
        freqs[i] = maf(freqs[i],360)

    freq_pos = [set() for i in xrange(360)]
    j = 0
    for freq in freqs:
        freq_pos[int(freq)-1].add(j)
        j += 1

    # now sample params from priors
    # note that these priors have nothing to do with the ABC
    tau = np.random.beta(0.35,0.35)  # tau for pheno model
    rho = np.random.beta(0.35,0.35)  # rho for pheno model
    p_b = 0.02+0.18*np.random.random()  # proportion of sites that are causal, no more than 50%, no less than 10%
    #p_b = 0.1
    p_s = np.random.random()  # proportion of sites drawn from Torg selection dist
    # now get effect sizes for each site in file
    num_sites = len(genos[0])
    num_causal = 0

    if num_sites < 100:
        print "not enough sites"
        exit(-1)

    while num_causal == 0:
        num_causal = np.random.binomial(num_sites,p_b)

    # now get all the effect sizes for num_causal variants
    num_torg = 0
    #while num_torg == 0:
    num_torg = np.random.binomial(num_causal,p_s)

    # get sel coeffs
    torg_file = os.path.join(os.path.expanduser('~'),'forRyan', 'tenTorgSel.txt')
    boyk_file = os.path.join(os.path.expanduser('~'),'forRyan', 'tennessenBoykoSel.txt')

    torf = open(torg_file,'r')
    boyf = open(boyk_file,'r')

    boy_s = []
    tor_s = []

    ln = 0
    for line in boyf:
        data = [float(x) for x in line.strip().split()]
        boy_s.append(data)    
        ln+=1
        if ln > linemax:
            break
    ln = 0
    for line in torf:
        data = [float(x) for x in line.strip().split()]
        tor_s.append(data)
        ln +=1 
        if ln > linemax:
            break

    torf.close()
    boyf.close()

    # map effect sizes to causal sites in the gene by frequency (greedy approach to deal with freq mismatch should be fine)
    # first get frequencies of all the sites
    boyk_coeffs_I = np.random.choice(len(boy_s), size=num_causal-num_torg)
    boyk_coeffs = [boy_s[i] for i in boyk_coeffs_I]

    torg_coeffs_I = np.random.choice(len(tor_s), size=num_torg)
    torg_coeffs = [tor_s[i] for i in torg_coeffs_I]

    # now find a variant in the actual data with same or as close as possible frequency for each variant
    coeffs = torg_coeffs+boyk_coeffs

    # transform s to effect sizes
    for var in coeffs:
        if random.random() > rho:
            var[1] = coeffs[np.random.randint(low=0,high=len(coeffs))][1]
        if random.random() > 0.5:
            var[1] = -1.*(abs(var[1])**tau)
        else:
            var[1] = abs(var[1])**tau

    # get a list of the number of sites of each frequency that are in the simulated data
    # also get a list of corresponding effect sizes
    coeff_sfs = np.zeros(360)
    coeff_effs = [[] for i in xrange(360)]
    for var in coeffs:
        ac = int(maf(var[0],360)-1)
        coeff_sfs[ac]+=1
        coeff_effs[ac].append(var[1])

    gen_positions = []
    gp_dict = {}
    j = 0
    for num_vars in coeff_sfs:
        if num_vars == 0:
            gen_positions.append([])
            j += 1
            continue
        if num_vars <= len(freq_pos[j]):
            genpos = random.sample(freq_pos[j],int(num_vars))
            gen_positions.append(genpos)
            #print genpos, len(genpos)
            #print coeff_effs[j], len(coeff_effs[j])
        # the below code block is for when there is no exact match in frequency ... then just add the extra variant(s) to the next frequency bin
        else:
            num_vars_t = len(freq_pos[j])
            genpos = random.sample(freq_pos[j],int(num_vars_t)) 
            gen_positions.append(genpos)
            if j < 359:
                while num_vars > len(freq_pos[j]):
                    num_vars -=1
                    coeff_sfs[j] -=1
                    coeff_sfs[j+1] += 1
                    coeff_effs[j+1].append(coeff_effs[j].pop())
        j += 1

    # calculate arch exactly and store somewhere
    # don't actually need arch but do need gp_dict
    var_exp  = np.zeros(360)
    for fr in xrange(len(gen_positions)):
        var_num = 0
        var_exp_fr = 0
        for site in gen_positions[fr]:
            gp_dict[site] = coeff_effs[fr][var_num]
            var_exp_fr += (coeff_effs[fr][var_num]**2)*(1-(fr+1.)/720.)*((fr+1.)/720.)
            var_num += 1
        var_exp[fr] = var_exp_fr
    
    var_exp = np.multiply(var_exp,1./np.sum(var_exp))
    
    #print rho, tau, p_s,
    #for thing in var_exp[0:10]:
    #    print thing,
    #print 
    #exit()

    # need to actually store all these parameters and data somewhere
    #print rho,tau, p_b, p_s, np.sum(var_exp)
    #for thing in var_exp:
    #    print thing

    # Simulate phenotypes, fix h^2 at 0.4 (observed mean), then rescale to be mean 0 and variance 1
    h2 = 0.09
    phenos = np.zeros(360)
    for dip in xrange(len(genos)):
        for var in gp_dict:
            phenos[dip] += genos[dip][var]*gp_dict[var]
        
    tot_var = np.var(phenos)

    # tot_var/(E + tot_var) = h2
    phenos = np.add(phenos, np.random.normal(scale = ((1-h2)*tot_var/h2)**0.5,size=len(phenos)))
    phenos = np.add(phenos,-np.mean(phenos)) 
    phenos = np.multiply(phenos, 1./np.var(phenos)**0.5)

    # write out to pheno file
    pid = sim
    pfile= "/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/phenos/"+gene_name+"/"+gene_name+"."+str(pid)+".txt"
    pdir= "/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/phenos/"+gene_name+"/"
    try:
        os.stat(pdir)
    except:
        os.mkdir(pdir)
    pf = open(pfile,'w')
    for thing in phenos:
        pf.write(str(thing)+"\n")
    pf.close()

    # Run HE on the phenotypes
    odir = "/netapp/home/lawrence.uricchio/ABC_herit/infer_group_genes/HE_regress_infer/1000000/K=20_MAC=1_ngPC=10_npPC=10_QN=1_DSAMP=1/"+gene_name
    try:
        os.stat(odir)
    except:
        os.mkdir(odir)
    myargs = ['Rscript', 'HEvLMM_gene_range_SNPlist_noResid.R', 'runHE=1', 'runLMM=0', 'SIM=0', 'ngPC=10', 'npPC=10','QN=1','NJOB='+str(pid),'GENE='+gene_name]
    p = subprocess.Popen(myargs)
    p.wait()

    # Store parameters for this gene
    parf.write(gene_name+' '+str(pid)+' '+str(rho)+' '+str(tau)+' '+str(p_s)+' '+str(p_b)+'\n')

parf.close()
