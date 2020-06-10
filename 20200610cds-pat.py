#!/usr/bin/env python
# coding: utf-8



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import bisect as bi
import argparse 

chromosome_ls =['I','II','III','IV','V','X']


parser = argparse.ArgumentParser(description='To mention sequence')
parser.add_argument('location_file',type=str)
parser.add_argument('counts_file',type=str)
parser.add_argument('cds',type=str)
parser.add_argument('pattern',type=str)
args = parser.parse_args()

def main():
    cds_ls = from_cds_file_to_cds_list(args.cds)
    hits_per_chromosome = from_genome_wide_hits_to_each_chromome(args.location_file,args.counts_file,chromosmenumber=6)
    distance_ls = distance_list(hits_per_chromosome,cds_ls,6)
    plotting(distance_ls)

def from_cds_file_to_cds_list(cdsfile): # change the file to read to argparse later
    df_cds = pd.read_table(cdsfile,usecols=[0,3,4,6,8],names=['chromosome','start','end','strand','annotation'])
    df_posi= df_cds.loc[df_cds['strand']=='+'].drop_duplicates(subset ='annotation',keep='first').drop(columns=['end','strand','annotation'])
    df_minus = df_cds.loc[df_cds['strand']=='-'].drop_duplicates(subset ='annotation',keep='last').drop(columns=['start','strand','annotation'])
    df_posi.rename(columns={'start':'position'},inplace=True,)
    df_minus.rename(columns={'end':'position'},inplace=True)
    cds = pd.concat([df_posi,df_minus],ignore_index=True)
    cds.sort_values(['chromosome','position'],ascending=True)
    cds_ls = [cds.loc[cds['chromosome']==k,'position'].copy().to_numpy() for k in chromosome_ls]
    return cds_ls

def from_genome_wide_hits_to_each_chromome(location,counts,chromosmenumber=6): 
    df_hits = pd.read_csv(location,sep='\t',names=['start','end','strand'])
    df_counts= pd.read_csv(counts,sep='\t',names=['hit_count'])
    total=0
    hits_per_chromosome=[0]*chromosmenumber
    for k in range(chromosmenumber):
        hits_per_chromosome[k] = df_hits[total:total+df_counts.loc[k]['hit_count']]
        total += df_counts.loc[k]['hit_count']
    return hits_per_chromosome

def distance_calculation (row,cds):
    if row['strand']=='+':
        hit_posi = row['start']
    elif row['strand']=='-':
        hit_posi = row['end']
    
    if np.amin(cds) < hit_posi < np.amax(cds):
        distance_from_downstream_cds = hit_posi - cds[bi.bisect_left(cds,hit_posi)]          
        distance_from_upstream_cds = hit_posi - cds[bi.bisect_left(cds,hit_posi)-1]
        if abs(distance_from_downstream_cds)<abs(distance_from_upstream_cds):
            distance = distance_from_downstream_cds
        else:
            distance = distance_from_upstream_cds
    elif hit_posi < np.amin(cds):
        distance = hit_posi - np.amin(cds)
    elif hit_posi > np.amax(cds):
        distance = hit_posi - np.amax(cds)
    elif hit_posi == np.amin(cds):
        distance=0
        
    if row['strand']=='+':
        return distance
    elif row['strand']=='-':
        return distance * -1
def distance_list(hits_ls,cds_ls,chromosomenumber=6):
    distance_ls = [list([distance_calculation(row[1],cds_ls[n]) for row in hits_ls[n].iterrows()]) for n  in range(chromosomenumber)]
    return distance_ls

def plotting(distance_ls):
    fig, axes = plt.subplots(3,2)
    for k,theaxis in enumerate(axes.flat): 
        theaxis.hist(distance_ls[k],100)
        theaxis.set_xlabel('Chromosome '+ chromosome_ls[k]) 
    plt.tight_layout()
    plt.suptitle('distance_from_{}_to_closest_ORF'.format(args.pattern))
    plt.subplots_adjust(top=0.9)
    plt.savefig('distance_from_{}_to_closest_ORF.pdf'.format(args.pattern))
    plt.show()




if __name__ == '__main__':
    main()

