#!/usr/bin/python
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("File", help="Description")
parser.add_argument("Virus_Coverage", help="Median virus coverage")
parser.add_argument("Virus_Name", help="Virus enter name of virus as per reference fasta (e.g. accession)")
parser.add_argument("--Min_Coverage", help="Minimum counts to include in calculation of Shannon Entropy")
parser.add_argument("-Quiet", action='store_true', help="Minimum counts")
args = parser.parse_args()
if args.Min_Coverage:
    Min_Coverage = int(args.Min_Coverage)
else:
    Min_Coverage = 0
if args.Quiet:
    Quiet = str(args.Quiet)
else:
    Quiet = False

File = str(args.File)
Virus_Coverage = int(args.Virus_Coverage)

Dicts = {}
with open(File, 'r') as File1:
    Data = File1.readline()
    while Data:
        Name = Data[13:-1]
        Dicts[Name] = File1.readline().split("\t")[:-1]
        Data = File1.readline()
        Data = File1.readline()

DictKeys = {}
n = 1
for Gene in Dicts:
    Data = Dicts[Gene]
    if args.Virus_Name in Gene:
        Total_Reads = Virus_Coverage
    else:
        Total_Reads = 0
        print("Unknown Gene")
    Sums = []
    RecTotal = 0
    for i in Data:
        data = i.split("_")
        Freq = int(data[-1])
        if Freq > Min_Coverage:
            RecTotal += Freq
            Sums.append(Freq)
        else:
            pass
    Entropy = 0
    for i in Sums:
        Fraction = i / float(RecTotal)
        Entropy -= math.log(Fraction, 2) * Fraction
    if not Quiet:
        print("Target Gene is: ", Gene)
        print("Entropy not considering WT: ", Entropy)
    else:
        print(Entropy)

    Entropy = 0
    for i in Sums:
        Fraction = i / float(RecTotal + Total_Reads)
        Entropy -= math.log(Fraction, 2) * Fraction
    Fraction = Total_Reads / float(RecTotal + Total_Reads)
    Entropy -= math.log(Fraction, 2) * Fraction
    if not Quiet:
        print("Entropy considering WT: ", Entropy)
