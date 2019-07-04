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

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--good-cols',
                        dest='good_cols_filename', required=True,
                        help='file with columns to extract (one column)')
    parser.add_argument('-p', '--probe-mapping',
                        dest='probe_map_filename', required=True,
                        help='file of probes to coordinates (two columns)')
    parser.add_argument('-m', '--methylation', dest='meth_filename',
                        required=True, help='table of methylation levels')
    parser.add_argument('-o', '--outfile', dest='outfile',
                        required=True, help='output file name')
    parser.add_argument('-v', '--verbose', action='store_true', dest='VERBOSE',
                        required=False, help= 'print extra run info')
    args = parser.parse_args()

    good_cols = []
    with open(args.good_cols_filename) as good_cols_file:
        for line in good_cols_file:
            (chrom, site) = line.strip().split()
            good_cols.append(chrom + ':' + str(int(site) + 1)) #because mapping in -m uses 1-based indexing

    found_cols = open('found_cols.txt', 'w+')
    probe_map = {}
    with open(args.probe_map_filename) as probe_map_file:
        for line in probe_map_file:
            (colname, site) = line.split()
            probe_map[colname] = site
            if site in good_cols:
                [chrom, num] = site.split(':')
                newnum = str(int(num)-1) 
                original_col = chrom + ':' + newnum
                found_cols.write(original_col + '\n')
    found_cols.close()

    out = open(args.outfile, 'w')
    meth = open(args.meth_filename)
    header = meth.readline()
    parts = header.strip().split()
    good_colnum = []
    count = 0
    for part in parts:
        try:
            if probe_map[part] in good_cols:
                good_colnum.append(count)
        except KeyError:
            continue
        count += 1
    good_values = []
    for j in good_colnum:
        good_values.append(parts[j])
    out.write('\t' + '\t'.join(good_values) + '\n')
    
    for line in meth:
        parts = line.strip().split()
        good_values = []
        for j in good_colnum:
            good_values.append(parts[j+1])
        out.write(parts[0] + '\t' + '\t'.join(good_values) + '\n')

if __name__ == '__main__':
    main()
