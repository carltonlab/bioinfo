#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 15:25:22 2020

@author: taka
"""

import numpy as np
import pandas as pd
import argparse 
import bisect as bi
import matplotlib.pyplot as plt

#all the global variable
plus_strand_hits = ['chromosomeI_hit_plus','chromosomeII_hit_plus','chromosomeIII_hit_plus','chromosomeIV_hit_plus','chromosomeV_hit_plus','chromosomeX_hit_plus']
minus_strand_hits = ['chromosomeI_hit_minus','chromosomeII_hit_minus','chromosomeIII_hit_minus','chromosomeIV_hit_minus','chromosomeV_hit_minus','chromosomeX_hit_minus']
plus_strand_cds = ['chromosomeI_plus_cds','chromosomeII_plus_cds','chromosomeIII_plus_cds','chromosomeIV_plus_cds','chromosomeV_plus_cds','chromosomeX_plus_cds']
minus_strand_cds = ['chromosomeI_minus_cds','chromosomeII_minus_cds','chromosomeIII_minus_cds','chromosomeIV_minus_cds','chromosomeV_minus_cds','chromosomeX_minus_cds']
chromosome_list = ['I','II','III','IV','V','X']


# this program 
#input: files that contains matches posi strand etc,Hitcounts file, processed cds file, searched pattern
#output: histogram to show the distane between searched pattern and ORF

def main():
    args = file_input()
    df_loc, df_hits, df_csv = file_read(args.location_file,args.counts_file,args.cds)
    d = from_genome_wide_hits_to_each_chromome(df_hits,df_loc) # d has 6keys for each chromosome, value dtype = pandas.DataFrame
    d2 = from_chrome_hits_to_chrome_plus_minus(d) # d2 has 12 keys; +/- for each chromosome, value dtype = ndarray
    d3 = from_cds_ls_to_chrome_plus_minus(df_csv) # d3 has 12 keys; +/- for each chromosome, value dtype = ndarray
    d5 = from_hits_to_hist_data(d2,d3) # d4 12 keys; +/- for each chromosome
    from_hist_data_to_histogram(d5,args.pattern)
    print(d5)

def file_input():
    parser = argparse.ArgumentParser(description='To mention sequence')
    parser.add_argument('location_file',type=str)
    parser.add_argument('counts_file',type=str)
    parser.add_argument('cds',type=str)
    parser.add_argument('pattern',type=str)
    args = parser.parse_args()
    return args

def file_read(hit_loc,hit_counts,cds): #args.location_file, args.counts_file
    df_loc = pd.read_csv(hit_loc,header = None, sep = '\t')
    df_hits = pd.read_csv(hit_counts,header = None, sep = '\t')
    df_csv = pd.read_csv(cds,header = None, sep = '\t')
    return df_loc, df_hits, df_csv

# To separate Hits into each chromosome
def from_genome_wide_hits_to_each_chromome(df_hits,df_loc): 
    d={}
    total=0
    for k,i in enumerate(chromosome_list):
        d["chromosome{}_hit".format(i)] = df_loc.iloc[total:total + df_hits.iloc[k][0]]
        total += df_hits.iloc[k][0]
    return d

#to separate +/- of pattern search results and remove end point of Hit accordinf to +/- value, leaving only one column
def from_chrome_hits_to_chrome_plus_minus(d): 
    d2 = {}
    for k, value in enumerate(d.values()):
        dataframe_posi = d["chromosome{}_hit".format(chromosome_list[k])].groupby(2).get_group('+').drop([1,2],axis=1)
        d2["chromosome{}_hit_plus".format(chromosome_list[k])] = dataframe_posi.iloc[:,0].to_numpy()
        dataframe_minus = d["chromosome{}_hit".format(chromosome_list[k])].groupby(2).get_group('+').drop([1,2],axis=1)
        d2["chromosome{}_hit_minus".format(chromosome_list[k])] = dataframe_minus.iloc[:,0].to_numpy()
    return d2

#to separate +/- of the cds file and remove end point of Hit accordinf to +/- value, leaving only one column
def from_cds_ls_to_chrome_plus_minus(df_csv):
    d3={}
    for k,i in enumerate(chromosome_list):
        cds_plus = df_csv.groupby(0).get_group(i).groupby(2).get_group('+').drop([0,2],axis=1)
        d3["chromosome{}_plus_cds".format(i)] = cds_plus.iloc[:,0].to_numpy()
        cds_minus = df_csv.groupby(0).get_group(i).groupby(2).get_group('-').drop([0,2],axis=1)
        d3["chromosome{}_minus_cds".format(i)] = cds_minus.iloc[:,0].to_numpy()
    return d3

#input:iterator(hit_loc) and chromo_cds output: distance corresponding to hit_loc
def from_hitposi_to_distance(hit_loc,chromo_cds):
    distance_ls=np.array([])
    for hit_posi in hit_loc:
        if hit_posi < np.amin(chromo_cds):
            distance = hit_posi - np.amin(chromo_cds)
            distance_ls = np.append(distance_ls, distance)
        
        if hit_posi == np.amin(chromo_cds):
            distance_ls = np.append(distance_ls, 0)
        
        if np.amin(chromo_cds) < hit_posi < np.amax(chromo_cds):
            distance_from_downstream_cds = hit_posi - chromo_cds[bi.bisect_left(chromo_cds,hit_posi)]          
            distance_from_upstream_cds = hit_posi - chromo_cds[bi.bisect_left(chromo_cds,hit_posi)-1] 
            if  distance_from_upstream_cds < distance_from_downstream_cds * -1:
                distance = distance_from_upstream_cds
                distance_ls = np.append(distance_ls, distance)
            elif distance_from_upstream_cds > distance_from_downstream_cds * -1: 
                distance = distance_from_downstream_cds
                distance_ls = np.append(distance_ls, distance)
            elif  distance_from_upstream_cds == distance_from_downstream_cds * -1: 
                distance = distance_from_upstream_cds
                distance_ls = np.append(distance_ls, [distance,-1*distance])
        
        if hit_posi == np.amax(chromo_cds):
            distance_ls = np.append(distance_ls, 0)
            
        if hit_posi > np.amax(chromo_cds):
            distance = hit_posi - np.amax(chromo_cds) 
            distance_ls = np.append(distance_ls, distance)
    return distance_ls

def from_hits_to_hist_data(d2,d3):
    d4={}
    d5={}
    for i in range(6):
        d4['minimumdistance{}_plus'.format(chromosome_list[i])] = from_hitposi_to_distance(d2[plus_strand_hits[i]],d3[plus_strand_cds[i]])
        d4['minimumdistance{}_minus'.format(chromosome_list[i])] = from_hitposi_to_distance(d2[minus_strand_hits[i]],d3[minus_strand_cds[i]]) * -1 # because upstream/downstream is reversed in minus strands
        d5['minimumdistance{}'.format(chromosome_list[i])] = np.append(d4['minimumdistance{}_plus'.format(chromosome_list[i])],d4['minimumdistance{}_minus'.format(chromosome_list[i])])
    return d5

def from_hist_data_to_histogram(d5,pattern):
    fig, axes = plt.subplots(3,2)
    
    for k,theaxis in enumerate(axes.flat): 
        theaxis.hist(d5['minimumdistance{}'.format(chromosome_list[k])],200)
        theaxis.set_xlabel('Chromosome '+ chromosome_list[k]) 
    fig.suptitle("closest startpoint of exon from {}".format(pattern))
    plt.tight_layout()
    path = 'distance_{}_to_startcodon.pdf'.format(pattern)
    plt.subplots_adjust(top=0.9)
    plt.savefig(path)
    plt.show()

if __name__ == '__main__':
    main()
    
    
