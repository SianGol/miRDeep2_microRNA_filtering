import re
import collections
from collections import defaultdict


def check_read_stack(sam_file, coords_file):
    # first make nested dict using defaultdict from collections, values are stored as list
    precursor_dict = defaultdict(lambda: defaultdict(list))
    # regions_dict function creates nested dict
    # chr: {strand: {precursor_id_half: [precursor half start, precursor half end]}}
    split_regions = regions_dict(coords_file)

    with open(sam_file) as sam:
        for s in sam:
            # skip any headers in sam file
            if '@' in s:
                continue
            chrom = s.strip().split('\t')[2]
            strand = s.strip().split('\t')[1]
            # determine if read on + or - strand as this determines co-ordinate of 5' start
            if strand == "+":
                read_start = int(s.strip().split('\t')[3])
                # check each read against any precursors on the same chr and strand as read
                # id_coords is {precursor_id_half : [precursor half start, precursor half end]}
                try:
                    id_coords = split_regions[chrom][strand]
                except KeyError:
                    continue

                for ID, loc in id_coords.items():
                    # check if read maps within precursor half and upstream by 10 nts
                    if (loc[0] - 10) <= read_start <= loc[1]:
                        precursor_dict[strand][ID].append(read_start)

            else:
                read_start = int(s.strip().split('\t')[3]) + (len(str(s.strip().split('\t')[4])) - 1)
                try:
                    id_coords = split_regions[chrom][strand]
                except KeyError:
                    continue

                for ID, loc in id_coords.items():
                    # check if read maps within precursor half and upstream by 10 nts (on - strand so upstream is +)
                    if loc[0] <= read_start <= (loc[1] + 10):
                        precursor_dict[strand][ID].append(read_start)

    precursor_dict.default_factory = None

    # make a list for precursors to append to which pass 50% homogeneity test
    precursors = []
    for strnd, nested_items in precursor_dict.items():
        for precursor_id_half, read_starts in nested_items.items():
            # count_list counts frequencies of read_start in list, are more than 50% the same?
            x = count_list(strnd, precursor_id_half, read_starts)
            precursors.append(x)

    # check if precursor exists in the above list in duplicate as this means stacked reads on both 5p and 3p arm
    precursors_filt = list(filter(None, (item for item, count in collections.Counter(precursors).items() if count > 1)))

    return precursors_filt


# create dictionary containing chr: {strand: {precursor_half: [start, end]}} for each half of predicted precursor
def regions_dict(infile):
    # first make nested dictionary using default dict
    dict_split_precursor = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    # open miRDeep2 results file (make sure tab separated)
    with open(infile) as file:
        for f in file:
            # only do for predicted precursors
            if 'chr' in f and 'mmu' not in f:
                split_line = f.strip().split('\t')

                name = split_line[0]
                chr = name.split('_')[0]
                strand = split_line[16].split(':')[2]
                precursor_start = int(split_line[16].split(':')[1].split('..')[0])
                precursor_end = int(split_line[16].split(':')[1].split('..')[1])
                length = int(precursor_end - precursor_start)
                start_first_half = precursor_start
                start_second_half = precursor_start + round(length/2)

                dict_split_precursor[chr][strand][name + "_half_1"] = [start_first_half, start_second_half - 1]
                dict_split_precursor[chr][strand][name + "_half_2"] = [start_second_half, precursor_end]

    dict_split_precursor.default_factory = None
    return dict_split_precursor


# count frequencies of read_start/+1 in list, are more than 50% the same?
def count_list(strand, precursor_id, lst):
    # check if list length is odd or even
    # if odd need to round up by adding 1 ('round' rounds length of an odd list of numbers down not up)
    if (len(lst) % 2) == 0:
        x = round(len(lst)/2)
    else:
        x = round(len(lst)/2) + 1

    for i in lst:
        if strand == "+":
            upper = i + 1
        else:
            upper = i - 1
        # does a read start account for >50% of reads?
        if (lst.count(i) + lst.count(upper)) >= x:
            full_id = str(precursor_id)
            # remove '_half_1/_half_2' from ID
            return full_id[:-7]


# keep precursors with folding free energy <= -0.2 kcal/mol/nt
def ffe_check(precursor_string_file):
    with open(precursor_string_file) as precursor_str:
        # make dict with precursor ID:[seq, bracket-dot notation] for every 3 lines of precursors.str
        pass_test = {}
        count = 0
        for p in precursor_str:
            count += 1
            if count == 1:
                p_id = p.strip().split(">")[1]
            if count == 2:
                sequence = p.strip()
            if count == 3:
                notation = p.strip().split(" ")[0]
                folding_energy = float(re.search(r"[+-]?\d+(?:\.\d+)?", p).group())
                folding_energy_nt = folding_energy/len(notation)
                count = 0

                if float(folding_energy_nt) <= -0.2:
                    pass_test[p_id] = [sequence, notation]

        return pass_test


# check for 0-4nt overhang
def find_overhang(pre, mat_fp, mat_tp):
    # first part of this function finds out if 3' of 3p arm has 0-4nt overhang
    # get paired brackets dict
    pairing = pair_dict(pre)
    # get position of first nt of 5p arm and check if paired, if not, find the next paired nt and get position
    if mat_fp.index("(") == 0:
        mat_fp_start = pre.index(mat_fp)
        # get position of last nt of 3p arm
        mat_tp_end = pre.index(mat_tp) + len(mat_tp)
    else:
        # adjust by adding the position of first bracket in precursor to the start position of mature 5p arm
        mat_fp_start = pre.index(mat_fp) + pre.index("(")
        # get position of last nt of 3p arm but adjust it by same the distance as above
        mat_tp_end = pre.index(mat_tp) + len(mat_tp) - pre.index("(")

    if mat_fp_start in pairing.keys():
        # -4 staggered end, blunt end
        if mat_tp_end - 4 <= pairing[mat_fp_start] <= mat_tp_end:
            # second part of function finds out if 3' of 5p arm has 0-4nt overhang
            # have to flip ) in index check
            if mat_tp.index(")") == 0:
                mat_tp_start = pre.index(mat_tp)
                mat_fp_end = pre.index(mat_fp) + len(mat_fp)
            else:
                mat_tp_start = pre.index(mat_tp) + mat_tp.index(")")
                mat_fp_end = pre.index(mat_fp) + len(mat_fp) - mat_tp.index(")")

            for k, v in pairing.items():
                if mat_tp_start == v:
                    match2 = k
                    if mat_fp_end - 4 <= match2 <= mat_fp_end:
                        return True


# number each forward/reverse paired bracket
def find(s, c):
    return (i for i, ltr in enumerate(s) if ltr == c)


# find which nucleotide pairs with which in hairpin
def pair_dict(str):
    x = find(str, "(")
    y = find(str, ")")
    # reverse x list (5p arm) so in same order as 3p
    xr = x[::-1]
    # make dictionary for complementary positions
    dict = {}
    for m, n in list(zip(xr, y)):
        dict[m] = n
    return dict


if __name__ == "__main__":

    from pathlib import Path
    #import os.path

    home = str(Path.home())
    # INPUT FILES
    # samfile of fasta reads (combined from all samples) mapped to genome with 0 mismatches in first 18nt and then 1 mismatch after that
    # generated using bowtie –n 0 –e 40 –l 18 –a –m 5 –best –strata (adjusted from mirdeep2 pipeline)
    samfile = f"{home}/bowtieoutput.sam"
    # coords file from miRDeep2 output file filtered for predicted precursors with >= 10 reads on both precursor arms
    coords = f"{home}/mirdeep2_filtered_min_10reads.txt"
    # precursors bracket-dot notation from mirdeep2 output
    precursors_str = f"{home}/mirdeep2_precursors.str"
    outfile = open(f"{home}/pass_filters_microRNAs.txt", "w")

    # get list of precursors which have stacked reads
    stacked = check_read_stack(samfile, coords)
    # get dict with { precursors: [seq, notation] } that pass ffe energy check
    ffe_pass = ffe_check(precursors_str)

    with open(coords) as coords:
        overhang_list = []
        for c in coords:
            name = str(c.strip().split('\t')[0])
            mature = c.strip().split('\t')[13].upper()
            star = c.strip().split('\t')[14].upper()
            precursor = c.strip().split('\t')[15].upper()

            # check if precursor has stacked reads
            if name in stacked:
                # get full sequence and bracket-dot notation from precursors_str file
                # this removes any that didnt pass ffe check
                try:
                    seq, notation = ffe_pass[name]
                except KeyError:
                    continue

                # remove any dodgy precursors made of repeated seqs that have got past stacked test
                if seq.count(mature) > 1 or seq.count(star) > 1:
                    continue

                # get dot-bracket notation of actual precursor in coords file
                precursor_not = notation[seq.index(precursor):seq.index(precursor) + len(precursor)]

                # find out if 'mature' defined by mirdeep2 is actually on 5p or 3p arm of precursor
                if precursor.index(mature) == 0:
                    fp_not = notation[seq.index(mature):seq.index(mature) + len(mature)]
                    tp_not = notation[seq.index(star):seq.index(star) + len(star)]
                else:
                    fp_not = notation[seq.index(star):seq.index(star) + len(star)]
                    tp_not = notation[seq.index(mature):seq.index(mature) + len(mature)]

                overhang = find_overhang(precursor_not, fp_not, tp_not)

                # check if overhang == True
                if overhang:
                    outfile.write(name + '\n')

    outfile.close()


