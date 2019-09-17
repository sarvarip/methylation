import pandas as pd
import numpy as np

ref = pd.read_csv('proms_merged_sorted', header=None, sep='\t')
result = []
gene_array = []

def isbigger(str1, str2):
    # is str1 bigger (in chromosomal order)
    # than str2
    if str1 == 'chrY':
        return True
    if str2 == 'chrY':
        return False
    if str1 == 'chrX': # we know that str2 is not chrY
        return True
    if str2 == 'chrX': # we know that str1 is not chrY
        return False
    if str1 == 'chrM': # we know that str2 is neither chrY or chrX
        return True
    if str2 == 'chrM': # we know that str1 is neither chrY or chrX
        return False
    if len(str2) > len(str1):
        return False
    elif len(str1) > len(str2):
        return True
    else:
        return (str1>str2)

with open('sorted_EPIC_probe_coords') as f:
    start_compare = 0
    found = 0
    for line in f:
        (coord, location) = line.strip().split('\t')
        (chrom, site) = location.split(':')
        print('Chromosome and site for probe: ' + chrom + ', ' + site)
        site = int(site)
        ref_chrom = ref.iloc[start_compare, 0]
        ref_val_lower = int(ref.iloc[start_compare, 1])
        ref_val_upper = int(ref.iloc[start_compare, 2])
        if site >= ref_val_lower and site <= ref_val_upper and chrom == ref_chrom:
            # we might think it has to be bigger than ref_val_lower,
            # because lines are sorted, but we need to check lower limit
            # if continue was called from while loop, which only check
            # for upper limit!
            # need to check chromosome because site might be smaller
            # than last end of gene on previous chromosome
            # E.g. it might happen that the last probe coord
            # on a chromosome matches a gene and then the next
            # probe is from a new chromosome but the ref chrom is
            # still the previous one
            if found == 1:
                result[-1].append(coord)
                gene_array[-1].append(start_compare)
                print('Probe belongs to same gene as probe before')
                print(start_compare)
                found = 1
                continue
            else:
                result.append([coord])
                gene_array.append([start_compare])
                print('New gene found for probe')
                print(start_compare)
                found = 1
                continue
        while site > ref_val_upper and chrom == ref_chrom:
            # need to check chrom in case site is bigger than the end
            # of the last gene on that chromosome (in which case
            # it would get stuck in loop if we did not compare chroms)
            start_compare = start_compare + 1
            ref_chrom = ref.iloc[start_compare, 0]
            while len(ref_chrom) > 5:
                start_compare = start_compare + 1
                ref_chrom = ref.iloc[start_compare, 0]
            ref_val_lower = int(ref.iloc[start_compare, 1])
            ref_val_upper = int(ref.iloc[start_compare, 2])
        if site >= ref_val_lower and chrom == ref_chrom:
            result.append([coord])
            gene_array.append([start_compare])
            print('New gene found for probe')
            print(start_compare)
            found = 1
            continue
        else:
            if chrom != ref_chrom:
                if isbigger(ref_chrom, chrom):
                    print('Probe was after the last gene in the chromosome')
                    print('Reference chrom (gene): ' + str(ref_chrom))
                    print('Chrom of probe: ' + str(chrom))
                    found = 0
                    continue
                else:
                    while site > ref_val_upper or chrom != ref_chrom:
                        start_compare = start_compare + 1
                        ref_chrom = ref.iloc[start_compare, 0]
                        while len(ref_chrom) > 5:
                            start_compare = start_compare + 1
                            ref_chrom = ref.iloc[start_compare, 0]
                        ref_val_lower = int(ref.iloc[start_compare, 1])
                        ref_val_upper = int(ref.iloc[start_compare, 2])
                    if site >= ref_val_lower:
                        result.append([coord])
                        gene_array([start_compare])
                        print('New gene on a new chromosome found for probe')
                        found = 1
                        continue
                    else:
                        print('First probe on a new chromosome has no corresponding genes')
                        found = 0
                        continue
            else:
                # was in between genes
                print("Probe is in between genes")
                found = 0
                continue

with open('coord_lists_per_gene', 'w+') as o:
    for element in result:
        o.write('\t'.join(element) + '\n')

with open('gene_array', 'w+') as o:
    for genes in gene_array:
        o.write('\t'.join([str(x) for x in genes]) + '\n')










