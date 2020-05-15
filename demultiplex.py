# split-seq demultiplexing

from __future__ import print_function

from Bio import Seq, SeqIO

import argparse

import pandas as pd

import gzip

import os, sys, time, traceback


######################################################### utility #########################################################
# version
def version():
    return '0.1.0'

# get current time
def time_now():
    return time.strftime('[%a %b-%d-%Y %X %Z]')
###########################################################################################################################

######################################################### barcode #########################################################
# dna alphabet
def alphabet():
    return ['A', 'C', 'G', 'T', 'N']

# load barcode file
# barcode file has 3 columns: barcode_to, barcode_from, time_point
# header is True and columns are separated by \t
# no missing values
# 
# barcode_file: a string of barcode file path
# returns a pandas dataframe of barcode info
def load_barcode_file(barcode_file):
    barcode = pd.read_csv(barcode_file, sep='\t')
    barcode.columns = ['barcode_to', 'barcode_from', 'time_point']

    return barcode

# generate barcode sequences within 1 mismatch with input barcode
# barcode: a string of barcode
# returns a list of string containing new barcode sequences within 1 mismatch
# with input barcode (including original barcode)
def generate_permutations(barcode):
    base_alphabet = alphabet()
    barcode_list = list(barcode) # string is immutable, so change to list

    # check for invalid character in barcode
    for i, base in enumerate(barcode_list):
        if base not in base_alphabet:
            print('Invalid character %s found in barcode: %s. Converting to N.' %(base, barcode), 
                file=sys.stderr)
            barcode_list[i] = 'N'

    barcode = ''.join(barcode_list)
    perms = [barcode] # include original barcode

    # permutation
    for i, base in enumerate(barcode_list):
        for entry in base_alphabet:
            if entry != base: # permutation to other bases
                barcode_list[i] = entry
                perms.append(''.join(barcode_list))

                barcode_list = list(barcode) # reset barcode list

    return perms

# wrapper
# barcode_file: a string of barcode file path
# returns 3 dictionaries: 
# barcode_all_dict: key is permutated sequence with <= 1 mismatch of barcode, value is original barcode
# barcode_replace_dict: key is barcode_from seq, value is barcode_to seq
# time_point_dict: key is barcode_to_seq, value is time_point
def get_barcode_dict(barcode_file):
    try:
        barcode_all_dict = {}
        barcode_replace_dict = {}
        time_point_dict = {}

        barcode = load_barcode_file(barcode_file) # load barcode file
        barcode = barcode.values.tolist() # to list

        for barcode_to, barcode_from, time_point in barcode:
            # permutation
            perms_barcode_to = generate_permutations(barcode_to)
            perms_barcode_from = generate_permutations(barcode_from)

            # barcode_to_dict
            for perm_barcode_to in perms_barcode_to:
                barcode_all_dict[perm_barcode_to] = barcode_to # key is permutation, val is original seq

            # barcode_from_dict
            for perm_barcode_from in perms_barcode_from:
                barcode_all_dict[perm_barcode_from] = barcode_from # key is permutation, val is original seq

            # barcode_replace_dict
            barcode_replace_dict[barcode_from] = barcode_to

            # time_point_dict
            time_point_dict[barcode_to] = time_point

        return barcode_all_dict, barcode_replace_dict, time_point_dict, 0

    except: # fail
        traceback.print_exc()
        return None, None, None, None, -1
###########################################################################################################################

##################################################### quality control #####################################################
# calcualte mean of a list
# in_list: input list of ints/floats
# returns a float of mean of list
def average(in_list):
    return 1.0 * sum(in_list) / len(in_list)

# count the occurrence of elements in a list less than a value
# in_list: input list of ints/floats
# val: an int/float cutoff
# returns an integer of counts where elements in list less than given value
def less_than(in_list, val):
    count = 0

    for entry in in_list:
        if entry <= val:
            count += 1

    return count

# open file for read or write
# infile_path: a string of input file path
# mode: file open mode
# returns a file handler
def safe_open(infile_path, mode='r'):
    if infile_path.endswith('.gz'):
        return gzip.open(infile_path, mode)
    else:
        return open(infile_path, mode)

# fastq quality format
# args: input command line arguments
# returns a string of fastq quality format
def fq_qual(args):
    qual = 'fastq'

    if args.phredQualtiy == 64:
        qual = 'fastq-illumina'

    return qual

# check fastq quality
def check_qual(phred_qual_list, args):
    # filter if average quality score is too low
    check1 = average(phred_qual_list) >= args.averageQuality
    # filter if too many low quality bases
    check2 = less_than(phred_qual_list, args.minQuality) <= args.minPercent * len(phred_qual_list) / 100.0

    return check1 and check2

# check barcode
# barcode: input string of barcode sequence
# barcode_dict: all possible barcode sequences within 1 mismatch
# returns 0 if the original barcode, 1 if 1 mismatch, 2 otherwise
def check_barcode(barcode, barcode_dict):
    if barcode in barcode_dict:
        if barcode == barcode_dict[barcode]: # original barcode
            return 0
        else: # 1 mismatch
            return 1
    else:
        return 2

# wrapper of quality control for paired-end sequencing
# args: input command line arguments
# barcode_all_dict: a dictionary of sequences with less than or equal to 1 mismatch
# barcode_replace_dict: a dictionary of barcode to replace (from key to val)
# compared with the original barcodes
# returns 0 if success and -1 otherwise
def run_qc_pe(args, barcode_all_dict, barcode_replace_dict, time_point_dict):
    try:
        fq_format = fq_qual(args) # get fastq format

        fq1_in = safe_open(args.fastq1, 'r')
        fq2_in = safe_open(args.fastq2, 'r')

        fq1 = SeqIO.parse(fq1_in, fq_format)
        fq2 = SeqIO.parse(fq2_in, fq_format)

        time_points = time_point_dict.values() # values are time points
        out_fqs = {}
        for entry in time_points:
            out_fqs[entry] = [
                safe_open('%s.%s.fastq.gz' %(args.outFile1Prefix, entry), 'a'), 
                safe_open('%s.%s.fastq.gz' %(args.outFile2Prefix, entry), 'a')
            ]

        entry_count = 0
        pass_count = 0

        while True:
            try:
                fq1_record = fq1.next()
                fq2_record = fq2.next()
            except:
                break

            entry_count += 1
            if entry_count % 1000000 == 0:
                print('%s %d fastq entries processed' %(time_now(), entry_count),
                      file=sys.stderr)

            phred_qual_list1 = fq1_record.letter_annotations['phred_quality']
            phred_qual_list2 = fq2_record.letter_annotations['phred_quality']

            # umi and barcode
            # 1:10 UMI (1-based)
            umi_qual = phred_qual_list2[0:10]
            # 11:18 barcode3 (1-based)
            # 49:56 barcode2 (1-based)
            # 87:94 barcode1 (1-based)
            barcode3 = fq2_record.seq[10:18]
            barcode2 = fq2_record.seq[48:56]
            barcode1 = fq2_record.seq[86:94]

            # quality check
            qual_check1 = check_qual(phred_qual_list1, args)
            qual_check2 = check_qual(phred_qual_list2, args)
            # barcode check
            barcode_check3 = check_barcode(barcode3, barcode_all_dict)
            barcode_check2 = check_barcode(barcode2, barcode_all_dict)
            barcode_check1 = check_barcode(barcode1, barcode_all_dict)
            barcode_check = barcode_check3 + barcode_check2 + barcode_check1
            # umi check
            umi_check = average(umi_qual) > 10

            if (qual_check1 and qual_check2 and 
                barcode_check1 <= 1 and barcode_check2 <= 1 and barcode_check3 <=1 and 
                umi_check):
                # alter sequences and quality scores in barcode regions
                if barcode_all_dict[barcode1] in barcode_replace_dict:
                    fq2_new_barcode_region = (fq2_record.seq[0:10] + 
                        barcode_all_dict[barcode3] + 
                        barcode_all_dict[barcode2] + 
                        barcode_replace_dict[barcode_all_dict[barcode1]])

                    time_point = time_point_dict[barcode_replace_dict[barcode_all_dict[barcode1]]]
                else:
                    fq2_new_barcode_region = (fq2_record.seq[0:10] + 
                        barcode_all_dict[barcode3] + 
                        barcode_all_dict[barcode2] + 
                        barcode_all_dict[barcode1])

                    time_point = time_point_dict[barcode_all_dict[barcode1]]
                #fq2_seq_region = fq2_record.seq[94:]

                fq2_qual = fq2_record.letter_annotations['phred_quality']
                fq2_qual_barcode_region = fq2_qual[0:10] + fq2_qual[10:18] + fq2_qual[48:56] + fq2_qual[86:94]
                #fq2_qual_seq_region = fq2_qual[94:]

                # change sequences
                fq2_record.letter_annotations = {} # clear dict, will reconstruct later
                fq2_record.seq = fq2_new_barcode_region

                # reconstruct quality scores
                fq2_record.letter_annotations = {'phred_quality': fq2_qual_barcode_region}

                # output file
                fq1_out = out_fqs[time_point][0]
                fq2_out = out_fqs[time_point][1]

                fq1_out.write(fq1_record.format(fq_format))
                fq2_out.write(fq2_record.format(fq_format))

                pass_count += 1

        print('%s %d fastq entries processed in total, and %d entries passed filtering'
              %(time_now(), entry_count, pass_count), file=sys.stderr)

        fq1_in.close()
        fq2_in.close()

        for entry in time_points:
            out_fqs[entry][0].close()
            out_fqs[entry][1].close()

        return 0

    except:
        traceback.print_exc()
        return -1
###########################################################################################################################

########################################################## args ###########################################################
# command line arguments
def get_args():
    parser = argparse.ArgumentParser(description=' '.join([
        'split-seq paired-end quality control utility (ver. %s).' %(version()), 
        '<Shiyuan Chen, shc@stowers.org>'
    ]))

    parser.add_argument('-1', '--fastq1', type=str, required=True, 
        help='Input fastq file 1')
    parser.add_argument('-2', '--fastq2', type=str, required=True, 
        help='Input fastq file 2 (for paired-end)')
    parser.add_argument('-o1', '--outFile1Prefix', type=str, required=True, 
        help='Output fastq file 1 Prefix')
    parser.add_argument('-o2', '--outFile2Prefix', type=str, required=True, 
        help='Output fastq file 2 (for paired-end) Prefix')
    parser.add_argument('-b', '--barcode', type=str, required=True, 
        help='SplitSeq barcode file')
    parser.add_argument('-q', '--minQuality', type=int, default=5, 
        help='Minimum quality score to keep')
    parser.add_argument('-p', '--minPercent', type=int, default=50, 
        help='Minimum percent of bases that must have MINQUALITY quality')
    parser.add_argument('-m', '--mismatch', type=int, default=1, 
        help='Barcode maximum mismatches allowed')
    parser.add_argument('-a', '--averageQuality', type=float, default=10, 
        help='Average quality of reads less than AVERGEQUALITY will be discarded')
    parser.add_argument('-f', '--phredQualtiy', type=int, default=33, 
        help='Phred quality format')

    return parser, parser.parse_args()

# check validity of input arguments
# 1. required options should be set
# 2. mismatch should be 1
def check_args(parser, args):
    check1 = (
        (args.fastq1 and args.fastq2) and # input fastq files
        (args.barcode) and # barcode file
        (args.outFile1Prefix and args.outFile2Prefix) # output fastq prefixes
    )

    if not check1:
        print('Not all required options are set: --fastq1, --fastq2, --barcode, --outFile1Prefix, --outFile2Prefix', 
            file=sys.stderr)
        parser.print_help()

    check2 = (args.mismatch == 1)
    if not check2:
        print('Currently, we only support 1 barcode mismatch (--mismatch 1)', file=sys.stderr)

    return check1 and check2
###########################################################################################################################


########################################################## main ###########################################################
def main():
    parser, args = get_args()

    print('%s Pipeline started' %(time_now()), file=sys.stderr)
    args_return = check_args(parser, args)

    if args_return:
        print('%s Loading barcode file' %(time_now()), file=sys.stderr)
        barcode_all_dict, barcode_replace_dict, time_point_dict, barcode_return = get_barcode_dict(args.barcode)

        if barcode_return == 0: # barcode loaded successfully
            print('%s Filtering paired-end reads...' %(time_now()), file=sys.stderr)
            qc_return = run_qc_pe(args, barcode_all_dict, barcode_replace_dict, time_point_dict)

            if qc_return == 0: # qc runs successfully
                print('%s Successfully completed' %(time_now()), file=sys.stderr)
            else: # qc fails
                print('%s Job failed' %(time_now()), file=sys.stderr)
        else: # barcode fails
            print('%s Barcode file loading failed' %(time_now()), file=sys.stderr)

    print('%s Pipeline finished' %(time_now()), file=sys.stderr)

if __name__ == '__main__':
    main()
