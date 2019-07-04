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
from operator import itemgetter
import numpy as np

def counting_sort(arr):
    mn = min(arr)
    mx = max(arr)
    shift = mn
    length = len(arr)
    count = np.zeros(mx-mn+1, dtype=int)
    sorted_arr = np.zeros(length, dtype=int)
    ix = np.zeros(length, dtype=int)
    for i in np.arange(length):
        num = arr[i]
        count[num-shift] += 1
    print('Counting sort bins: ')
    print(count)
    sumCount = np.cumsum(count)

    for i in np.arange(length-1,-1,-1):
        num = arr[i]
        sorted_arr[sumCount[num-shift]-1] = num
        ix[sumCount[num-shift]-1] = i
        sumCount[num-shift] -= 1


    return ix, sorted_arr

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='infile',
                        required=True, help='input file to be sorted')
    parser.add_argument('-i2', '--input2', dest='aux',
                        required=False, help='auxiliary file to be aligned')
    parser.add_argument('-o', '--outfile', dest='outfile',
                        required=True, help='output file name where sorted input is stored')
    parser.add_argument('-o2', '--outfile2', dest='aligned_aux',
                        required=False, help='output file name where aligned auxiliary input is stored')
    parser.add_argument('-idx', '--sortedidx', dest='idx',
                        required=False, help= 'file to output the indices of the sorted array')
    parser.add_argument('-v', '--verbose', action='store_true', dest='VERBOSE',
                        required=False, help= 'print extra run info')
    args = parser.parse_args()

    chromarr = []
    sitearr = []

    with open(args.infile) as f:
      for line in f:
        (chrom, site) = line.strip().split(':')
        chromnum = chrom[3:]
        if chromnum == 'X':
            chromnum = 23
        elif chromnum == 'Y':
            chromnum = 24
        chromnum = int(chromnum)
        sitenum = int(site)
        chromarr.append(chromnum)
        sitearr.append(sitenum)
    chromarr = np.array(chromarr)
    sitearr = np.array(sitearr)
    indices = np.argsort(sitearr)
    site_sorted = sitearr[indices]
    chromarr = chromarr[indices]
    indx, chrom_sorted = counting_sort(chromarr)
    site_sorted = site_sorted[indx]
    indices = indices[indx]

    if args.aux is not None and args.aligned_aux is not None:
        aux_list = []
        with open(args.aux) as au:
            for line in au:
                coord = line.strip()
                aux_list.append(line)
        aux_list_aligned = list(np.array(aux_list)[indices])

    with open(args.outfile, 'w+') as o:
        for i in np.arange(len(site_sorted)):
            if chrom_sorted[i] == 23:
                curr_chrom = 'X'
            elif chrom_sorted[i] == 24:
                curr_chrom = 'Y'
            else:
                curr_chrom = str(chrom_sorted[i])
            o.write('chr' + curr_chrom + ':' + str(site_sorted[i]) + '\n')

    if args.idx is not None:
        with open(args.idx, 'w+') as idxfile:
            for index in indices:
                idxfile.write(str(index) + '\n')

    if args.aux is not None and args.aligned_aux is not None:
        with open(args.aligned_aux, 'w+') as al_au:
            for coord in aux_list_aligned:
                al_au.write(coord)


if __name__ == '__main__':
    main()
