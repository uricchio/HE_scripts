import sys
import numpy as np

# real.data is the input file

fh = open(sys.argv[1], 'r')
non_cumu = []
for line in fh:

    data = [float(x) for x in line.strip().split()]
    non_cumu.append(data[0])
    for i in range(1,len(data)):
        non_cumu.append(data[i]-data[i-1])
    for i in range(len(non_cumu)):
        non_cumu[i] /= data[-1]


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


my_bins = [1,2,3,5,8,11,17,25,36,50,68,89,114,143,174,207,242,280,319,360]

fin_data= apply_var(non_cumu,my_bins)

cumus = [fin_data[0]]
for i in range(1,len(fin_data)):
    cumus.append(fin_data[i]+cumus[i-1])

out_freqs = [0,1,2,4,9,19,59,119,179,239,359]

tot = cumus[0]
print cumus[0],

for i in range(1,len(out_freqs)):
    
    tot += cumus[out_freqs[i]]-cumus[out_freqs[i-1]]
    print cumus[out_freqs[i]]-cumus[out_freqs[i-1]],

print 

