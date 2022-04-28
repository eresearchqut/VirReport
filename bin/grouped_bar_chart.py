#!/usr/bin/env python
import pandas as pd
from functools import reduce
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

run_data = pd.DataFrame()
iterator=int(0)

for fl in glob.glob("*_read_length_dist.txt"):
    sample = (fl.replace('_read_length_dist.txt', ''))
    sample_data = pd.read_csv(fl, header=None, sep='\t',index_col=None)
    sample_data.columns = ["length", sample]
    if iterator == 0: 
        run_data = run_data.append(sample_data)
    else:
        run_data = pd.merge(run_data, sample_data, how="outer", on=["length"])
    iterator += 1
print(iterator)

#length = ()
length = int(round(iterator/4))
print(length)
run_data = run_data.reindex(sorted(run_data.columns), axis=1)
run_data = run_data.set_index('length')
print(run_data)
dim=len(run_data.columns)
print(dim)

#if the number of samples is not divisible by 4, add dummy columns
if (iterator % 4 != 0):
    leftover = int(((length+1)*4)-iterator)
    print(leftover)

    for i in range(1, leftover+1):
        run_data['XX'+ str(i)] = 1
    dim=len(run_data.columns)
    print(dim)
    new_length=length+1
    print(new_length)
    fig, a = plt.subplots(new_length, 4, figsize=(10, dim), tight_layout=True)
    #delete the dummy subplots
    for i in range(1, leftover+1):
        fig.delaxes(a[new_length-1][4-i])
    run_data.plot.barh(ax=a, subplots=True, fontsize=7)
#run_data = run_data.reindex(sorted(run_data.columns), axis=1)
#run_data = run_data.set_index('length')
#print(run_data)
#derive the height of the final PDF based on columns


elif iterator < 4:
    fig, a = plt.subplots(1, 4, figsize=(10, dim), tight_layout=True)
    print(leftover)
    print(length)
    for i in range(1, leftover+1):
        fig.delaxes(a[4-i])
    run_data.plot.barh(ax=a, subplots=True, fontsize=7)
#else:
    #if (iterator % 4 != 0):
    #    new_length=length+1
    #    print(new_length)
    #    fig, a = plt.subplots(new_length, 4, figsize=(10, dim), tight_layout=True)
    #    #delete the dummy subplots
    #    for i in range(1, leftover+1):
    #        fig.delaxes(a[length-1][4-i])
else:
    fig, a = plt.subplots(length, 4, figsize=(10, dim), tight_layout=True)
    run_data.plot.barh(ax=a, subplots=True, fontsize=7)

fig.savefig('run_read_size_distribution.pdf', format='pdf')
