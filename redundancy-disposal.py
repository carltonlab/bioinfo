#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 14:23:18 2020

@author: taka
"""

import pandas as pd 
import argparse


def main():
    args = get_input_from_bash()
    df_positive , df_negative = duplicates_removal(args.plus_file, args.minus_file)
    exons = concatenation(df_positive,df_negative)
    sorting(exons)
    exons.to_csv('c_elegans_ws274_startcodons.txt',index=False,sep='\t', header=False)

def get_input_from_bash ():
    parser = argparse.ArgumentParser(description='getting plus strand file and minus strand file')
    parser.add_argument('plus_file',type=str)
    parser.add_argument('minus_file',type= str)
    args = parser.parse_args()
    return args

def duplicates_removal(plus_file,minus_file):
    df_positive = pd.read_csv(plus_file, header = None, sep='\t')
    df_negative = pd.read_csv(minus_file, header = None, sep='\t')
    
    df_positive.drop_duplicates(subset=4, keep='first',inplace=True)
    df_negative.drop_duplicates(subset=4, keep='last', inplace=True) 
    return df_positive, df_negative

def concatenation(df_positive,df_negative):
    exon_dict = {'positive_starand':df_positive[[0,1,2]],'negative_strand':df_negative[[0,1,2]]}
    exons = pd.concat(exon_dict)
    return exons

def sorting(exons):
    exons.sort_values([0,1],ascending=[1,1],inplace=True)
    return exons

if __name__=="__main__":
    main()
