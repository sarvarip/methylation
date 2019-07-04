# Copyright (C) 2019 Peter Sarvari
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import sys, argparse
import pandas as pd
import numpy as np

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-m1', '--methylation1',
                        dest='meth1', required=True,
                        help='first mehylation dataset that is transposed')
    parser.add_argument('-m2', '--methylation2',
                        dest='meth2', required=True,
                        help='reduced mehylation dataset')
    parser.add_argument('-col', '--columns',
                        dest='col', required=True,
                        help='Produced by extract.py')
    parser.add_argument('-o', '--output',
                        dest='out', required=True,
                        help='Output for sample correlation data')
    args = parser.parse_args()

    df = pd.read_csv(args.meth1,sep='\t', header=0, index_col = 0)
    df = df.iloc[:,:-1]
    df = df.T
    found_cols = []
    with open(args.col) as f:
          for line in f:
                found_cols.append(line.strip())

    found_samples = []
    with open(args.meth2) as met:
          for line in met:
                parts = line.strip().split()
                found_samples.append(parts[0][1:])
    found_samples = found_samples[1:] #to get rid of header

    retained_sample = []
    for sample in found_samples:
          if sample in df.index:
                retained_sample.append(sample)

    df2 = pd.read_csv(args.meth2,sep='\t', header=0, index_col = 0)

    idx = df2.index
    newidx = []
    for i in idx:
          i = i[1:]
          #wat
          newidx.append(i)
    df2.index = newidx

    df2_redsample = df2.loc[retained_sample,:]
    df_red = df.loc[retained_sample, found_cols]
    df2_redsample = df2_redsample.T
    df_red = df_red.T

    correlations = []
    for i in np.arange(len(retained_sample)):
          arr1 = df2_redsample.iloc[:,i].values
          arr2 = df_red.iloc[:,i].values
          corr = np.corrcoef(arr1, arr2)
          correl = corr[0][1]
          correlations.append(correl)

    out = pd.DataFrame(data=correlations, index=retained_sample, columns = ['Sample correlations between the two methods'])
    print(out)
    out.to_csv(args.out, header=None, sep='\t')

if __name__ == '__main__':
      main()

