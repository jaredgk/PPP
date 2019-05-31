import sys
import argparse
import logging
import numpy as np

from pgpipe.model import Model, read_single_model

def readLocus(ima_f,pop_count):
    line = ima_f.readline()
    seqs = []
    seqs_per_pop = [int(i) for i in line.strip().split()[1:1+pop_count]]
    total_seqs = sum(seqs_per_pop)
    for i in range(total_seqs):
        line = ima_f.readline()
        seqs.append(line.strip().split()[1])
    return seqs,seqs_per_pop

def adjCoef(n):
    a = 0
    for i in range(1,max(2,int(n))):
        a += float(1)/float(i)
    return a

def wattersonPerPop(seqs,seqs_per_pop):
    seq_idx = 0
    watt_pop = []
    for i in range(len(seqs_per_pop)):
        varsites = 0
        ref_seq = seqs[seq_idx]
        for j in range(len(ref_seq)):
            for k in range(1,seqs_per_pop[i]):
                if ref_seq[j] != seqs[k+seq_idx][j]:
                    varsites += 1
                    break
        a = adjCoef(seqs_per_pop[i]/2)
        #print (varsites,a,seqs_per_pop[i])
        watt_pop.append(float(varsites)/float(a))
        seq_idx += seqs_per_pop[i]
    #print (watt_pop)
    #exit()
    return watt_pop

def statGeomean(locus_stats):
    means = []
    for i in range(len(locus_stats[0])):
        ll = np.array([float(1) if locus_stats[j][i] == 0 else locus_stats[j][i] for j in range(len(locus_stats))])
        logl = np.log(ll)
        means.append(np.exp(logl.sum()/len(logl)))
        #lmean = 0
        #for j in range(len(locus_stats)):
    return means


locus_stats = []

ima_f = open(str(sys.argv[1]),'r')

ima_f.readline()
line = ima_f.readline()
while line[0] == '#':
    line = ima_f.readline()
pop_count = int(line.strip())
line = ima_f.readline()
pop_names = line.strip().split()
if len(pop_names) != pop_count:
    raise Exception("Error parsing, length of population list is different from population count")
line = ima_f.readline()
if ':' in line:
    line = ima_f.readline()

locus_count = int(line.strip())
for i in range(locus_count):
    locus_data,seqs_per_pop = readLocus(ima_f,pop_count)
    locus_stats.append(wattersonPerPop(locus_data,seqs_per_pop))

gms = statGeomean(locus_stats)
print (gms)
ne = max(gms)
print (max(gms))
q = ne*5
t = ne*2
m = float(2)/float(ne)
print ("Q: ",q)
print ("T: ",t)
print ("M: ",m)