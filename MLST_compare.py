#!/usr/bin/env python3

#
# Description: Script to compare and group similar MLST profiles together to allow SNV analysis
#
# Usage: python3 ./MLST_compare.py -i input -o outbreak_folder_name
#
# Output location: parameter
#
# Modules required: None
#
# v1.0.2 (1/15/2020)
#
# Created by Rich Stanton (njr5@cdc.gov)
#

import sys
import glob
import itertools
import re
import argparse

#Usage: python MLST_Compare_SciComp_Exe.py Input_Samples_File Outbreak_Folder_Name
#Written by Rich Stanton njr5@cdc.gov

# Parse all arguments from command line
def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to group similar MLSTs together')
    parser.add_argument('-i', '--input', required=True, help='input samples file')
    parser.add_argument('-o', '--output', required=True, help='output/outbreak folder name')
    return parser.parse_args()


def List_Scorer(list1, list2):
    Similarity_Score = 0
    for positions in range(len(list1)):
        if list1[positions] == list2[positions]:
            Similarity_Score = Similarity_Score + 1
        else:
            continue
    return Similarity_Score

def MLST_List_Maker(MLST_file):
    """Makes a list of alleles from an MLST file"""
    f = open(MLST_file, 'r')
    String1 = f.readline()
    f.close()
    List1 = list(filter(None, re.split("[()\t\n]+", String1)))
    if len(List1) != 17:
        List1 = ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0']
    List2 = [[List1[4]], [List1[6]], [List1[8]], [List1[10]], [List1[12]], [List1[14]], [List1[16]]]
    List3 = []
    for items in List2:
        item = list(filter(None, re.split(',', items[0])))
        List3.append(item)
    Output_Lists = list(itertools.product(*List3))
    return Output_Lists

def MLST_List_Maker_Names(MLST_file):
    """Makes a list of alleles from an MLST file"""
    f = open(MLST_file, 'r')
    String1 = f.readline()
    f.close()
    List1 = list(filter(None, re.split("[()\t\n]+", String1)))
    if len(List1) != 17:
        List1 = ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0']
    List2 = [List1[3], List1[5], List1[7], List1[9], List1[11], List1[13], List1[15]]
    Allele_List = MLST_List_Maker(MLST_file)
    Out_List = []
    for Matches in Allele_List:
        Current_List = []
        for Alleles in range(0, 7):
            Current = List2[Alleles] + '_' + Matches[Alleles]
            Current_List.append(Current)
        Out_List.append(Current_List)
    return Out_List

def Max_MLST_Matcher(MLST1, MLST2):
    """Makes a max comparison score of 0-7 for # of overlapping MLST alleles even with multiple alleles"""
    MLST_List1 = MLST_List_Maker(MLST1)
    MLST_List2 = MLST_List_Maker(MLST2)
    Max = 0
    for mlst1 in MLST_List1:
        for mlst2 in MLST_List2:
            Score = List_Scorer(mlst1, mlst2)
            if Score > Max:
                Max = Score
    return Max

def Max_MLST_Matcher_Names(MLST1, MLST2):
    """Makes a max comparison score of 0-7 for # of overlapping MLST alleles even with multiple alleles"""
    MLST_List1 = MLST_List_Maker_Names(MLST1)
    MLST_List2 = MLST_List_Maker_Names(MLST2)
    Max = 0
    for mlst1 in MLST_List1:
        for mlst2 in MLST_List2:
            Score = List_Scorer(mlst1, mlst2)
            if Score > Max:
                Max = Score
    return Max

def Max_MLST_Entries(name1, name2):
    MLST_1 = glob.glob('/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/' + name1+ '/MLST/*.mlst')
    MLST_2 = glob.glob('/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/' + name2+ '/MLST/*.mlst')
    Max = 0
    for files1 in MLST_1:
        for files2 in MLST_2:
            Score = Max_MLST_Matcher(files1, files2)
            if Score > Max:
                Max = Score
    return Max

def MLST_Species(MLST_file):
    f = open(MLST_file, 'r')
    String1 = f.readline()
    f.close()
    Species = 'Unknown'
    List1 = list(filter(None, re.split("[()\t\n]+", String1)))
    if len(List1) == 17 or len(List1) == 19:
        Species = List1[1]
        Species = list(filter(None,re.split("_+", Species)))
        Species = Species[0]
    return Species

def Taxa_Stats(input_stats_file):
    f = open(input_stats_file, 'r')
    String1 = f.readline()
    Species = 'Unknown'
    while String1 != '':
        List1 = list(filter(None, re.split(" ", String1)))
        if List1[0] == 'Taxa':
            Species = List1[4][0].lower() + List1[5]
            if Species[-1] == '\n':
                Species = Species[0:-1]
            String1 = f.readline()
        else:
            String1 = f.readline()
    f.close()
    return Species

def MLST_ST(MLST_file):
    f = open(MLST_file, 'r')
    String1 = f.readline()
    f.close()
    ST = ['None']
    #print(String1)
    List1 = list(filter(None, re.split("[()\t\n]+", String1)))
    if len(List1) == 17 or len(List1) == 19:
        ST = List1[2]
        ST = list(filter(None,re.split("[/,]+", ST)))
        #print(ST)
    return ST

def Species_Cluster_Maker(input_list):
    """Reads in an input file list and returns a list of species overlaps"""
    Species_List = []
    Out_List = []
    for isos in input_list:
        Name = isos.split('/')[1]
        MLST = glob.glob('/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/' + isos + '/MLST/*.mlst')
        #print("IS:",isos, MLST[0])
        Species = MLST_Species(MLST[0])
        if Species == 'Unknown':
            Stats = '/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/' + isos + '/' + Name + '_pipeline_stats.txt'
            Species = Taxa_Stats(Stats)
        Iso_Info = [isos, Species]
        Species_List.append(Iso_Info)
    Out_List = []
    All_Species = []
    for items in Species_List:
        if (items[1] in All_Species) == False:
            All_Species.append(items[1])
    for species in All_Species:
        Current_Species = []
        for entries in Species_List:
            if entries[1] == species:
                Current_Species.append(entries)
        Out_List.append(Current_Species)
    return Out_List

def ST_Cluster_Maker(input_list):
    """Reads in an input file list and returns a list of ST overlaps"""
    ST_List = []
    Out_List = []
    for isos in input_list:
        MLST = glob.glob('/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/' + isos + '/MLST/*.mlst')
        ST = []
        for files in MLST:
            ST = ST + MLST_ST(files)
        Iso_Info = [isos, ST]
        ST_List.append(Iso_Info)
    Out_List = []
    Out_List.append([ST_List[0]])
    for items in ST_List[1:]:
        Add = 0
        for Matches in Out_List:
            for STs in items[1]:
                for entries in Matches:
                    if (STs in entries[1]) == True:
                        Matches.append(items)
                        Add = 1
                        break
        if Add == 0:
            Out_List.append([items])
    return Out_List

def Close_Cluster_Maker(input_isolate_list):
    """Reads in an input file list and returns a list of ST overlaps"""
    Output_Clusters = []
    Output_Clusters.append([input_isolate_list[0]])
    for isolates in input_isolate_list[1:]:
        Add = 0
        for clusters in Output_Clusters:
            Name1 = clusters[0]
            Name2 = isolates
            Max = Max_MLST_Entries(Name1, Name2)
            if Max >= 6:
                Add = 1
                clusters.append(isolates)
                break
        if Add == 0:
            Output_Clusters.append([isolates])
    return Output_Clusters

def Close_STs(input_clusters_list):
    """Expands cluster lists to include closely related isolates"""
    Output_Clusters = []
    Output_Clusters.append(input_clusters_list[0])
    for cluster1 in input_clusters_list[1:]:
        Add = 0
        for cluster2 in Output_Clusters:
            Name1 = cluster1[0][0]
            Name2 = cluster2[0][0]
            Max = Max_MLST_Entries(Name1, Name2)
            if Max >= 6:
                Add = 1
                New_Cluster = cluster1 + cluster2
                Output_Clusters.remove(cluster2)
                Output_Clusters.append(New_Cluster)
                break
        if Add == 0:
            Output_Clusters.append(cluster1)
    return Output_Clusters

def Same_STs(input_clusters_list):
    """Expands cluster lists to include isolates from the same ST"""
    Output_Clusters = []
    Output_Clusters.append(input_clusters_list[0])
    for cluster1 in input_clusters_list[1:]:
        Add = 0
        for cluster2 in Output_Clusters:
            Name1 = cluster1[0][0]
            Name2 = cluster2[0][0]
            Max = Max_MLST_Entries(Name1, Name2)
            if Max == 7:
                Add = 1
                New_Cluster = cluster1 + cluster2
                Output_Clusters.remove(cluster2)
                Output_Clusters.append(New_Cluster)
                break
        if Add == 0:
            Output_Clusters.append(cluster1)
    return Output_Clusters

def Sample_Maker(input_list, output_file):
    """Makes an output file of the isolates in a list"""
    Output = open(output_file, 'w')
    for items in input_list:
        Output.write(items[0] + '\n')
    Output.close()

def Repeat_Remover(any_list):
    """Removes repeats for any list"""
    new_list = []
    for items in any_list:
        if (items in new_list) == False:
            new_list.append(items)
    return new_list

def ST_List_Maker(input_isolate_list):
    """Makes a list of STs present from an isolate list"""
    ST_List = []
    for entries in input_isolate_list:
        Files = glob.glob('/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/' + entries + '/MLST/*.mlst')
        File = Files[0]
        MLST = MLST_ST(File)
        ST_List.append(MLST)
    return ST_List

def STs_Present(input_ST_list):
    """Makes a string of the STs present in a ST cluster list entry"""
    STs = []
    for entries in input_ST_list:
        STs = STs + entries[1]
    STs = Repeat_Remover(STs)
    Final = []
    for entries in STs:
        if str.isdigit(entries) == True:
            #print("is digit:"+entries)
            Final.append(int(entries))
        else:
            #print("Not digit:"+entries)
            Final.append(entries)
    Final.sort()
    Out = []
    for entries in Final:
        Out.append(str(entries))
    #Out_String = ','.join(Out)
    Out_String = '_'.join(Out)
    return Out_String

def OA_Samples(input_file, outbreak_name):
    """Reads an outbreak sample list and returns the relevant .sample files"""
    Sample_List = []
    f = open(input_file, 'r')
    String1 = f.readline()
    while String1 != '':
        if String1[-1] == '\n':
            Sample_List.append(String1[0:-1])
            String1 = f.readline()
        else:
            Sample_List.append(String1)
            String1 = f.readline()
    f.close()
    #print("SL:",Sample_List)
    Species = Species_Cluster_Maker(Sample_List)
##    for entries in Species:
##        Name= entries[0][1]
##        f = open(outbreak_name + '__' + Name + '.samples', 'w')
##        for isos in entries:
##            f.write(isos[0] + '\n')
##        f.close()
    for entries in Species:
        Name= entries[0][1]
        Isolate_List = []
        for isos in entries:
            Isolate_List.append(isos[0])
        STs = ST_Cluster_Maker(Isolate_List)
        #Clusters = Close_STs(STs)
        Clusters = Same_STs(STs)
        for clusters in Clusters:
            if len(clusters) > 1:
                ST_String = STs_Present(clusters)
                f = open(outbreak_name + '__' + Name + '__' + ST_String + '.samples', 'w')
                clusters = Repeat_Remover(clusters)
                for isos in clusters:
                    f.write(isos[0] + '\n')
                f.close()
        if len(Clusters) > 1:
            f = open(outbreak_name + '__' + Name + '.samples', 'w')
            entries = Repeat_Remover(entries)
            for isos in entries:
                f.write(isos[0] + '\n')
            f.close()

#Input = sys.argv[1]
#Outbreak_Name = sys.argv[2]

args = parseArgs()
OA_Samples(args.input, args.output)
