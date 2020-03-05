#!/usr/bin/env python3

#
# Description: Script to extract a subsequence from a fasta sequence
#
# Usage: python ./get_subsequence.py -i fasta_file -s start_position_of_wanted_sequence -e stop_position_of_wanted_sequence
#   -t title_of_fasta_if_multi_fasta-ed -r do_reverse_complement_on_sequence. I must recheck to see what the -b and -f do. It
#   may be even a subsequence of the subsequence
#
# Output location: standard out
#
# Modules required: Biopython must be available in python instance
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import sys
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input file', required=True, dest='input_file')
parser.add_argument('-s', '--start', help='start position', required=True, type=int, dest='start')
parser.add_argument('-e', '--end', help='end position', required=True, type=int, dest='end')
parser.add_argument('-t', '--title', help='header', required=False, type=str, dest='title')
parser.add_argument('-r', '--reverse', help='reverse complement strand?', required=False, type=bool, dest='reverse', default=False)
parser.add_argument('-b', '--begin', help='rev-comp begin position', required=False, type=int, dest='begin', default=0)
parser.add_argument('-f', '--finish', help='rev-comp finish position', required=False, type=int, dest='finish', default=0)
parser.add_argument('-o', '--header_out', help='header info to output', required=False, type=str, dest='header_out', default=0)
parameters=parser.parse_args()

# Counts the number of sequences in the file using > as the delimiter
with open(parameters.input_file, 'r') as myfile:
    data=myfile.read().replace('\n', ' ')
sequence_count=data.count(">")
#print("SC:", sequence_count)

#Sets the initial record to catch if it has not been found during the parsing
match_record=None

#Cahnges how to handle a single entry vs multi-entry fasta file
if sequence_count == 0:
    #print("No sequences found in input fasta file, exiting")
    exit()
elif sequence_count == 1:
    match_record = SeqIO.read(parameters.input_file, "fasta")
elif parameters.title == None:
    counter=0
    for record in SeqIO.parse(parameters.input_file, "fasta"):
        match_record = record
        #print("Using first sequence:", record.id)
        # Purposely left to break after first sequence is retrieved
        break
else:
    for record in SeqIO.parse(parameters.input_file, "fasta"):
        if parameters.title in record.id:
            match_record = record
            #print("Found match for:", parameters.title, "as", record.id)
            break
if match_record is None:
    print("No matching header found, exiting...")
#else:
    #print(match_record.id)
    #print(match_record.description)
    #print(match_record.seq)


sequence = ""
#match_record = SeqIO.read(sys.argv[1],"fasta")
#match_record = SeqIO.read(parameters.input_file, "fasta")
#start = int(sys.argv[2])-1
#end = int(sys.argv[3])
#print(":"+str(parameters.start)+":"+str(parameters.end)+":")
#print(match_record.id)
#print(match_record.description)

search_DNA_seq = match_record.seq[parameters.start:parameters.end]
reverse_record = SeqRecord(Seq(str(search_DNA_seq)), match_record.id, '', '')
#print(reverse_record.seq)
reverse_record=reverse_record.reverse_complement()
#print(reverse_record.seq)

if parameters.reverse:
    start_sub = parameters.begin
    if start_sub > 0:
        start_sub-=1
    if parameters.finish <= parameters.begin:
        print("Reverse finish is less than reverse begin???")
        exit()
    else:
        end_sub = parameters.finish
    print(reverse_record.seq)[start_sub:end_sub]
    exit()
if parameters.header_out is not None:
    print(">"+parameters.header_out+"\n"+search_DNA_seq)
else:
    print(search_DNA_seq)
