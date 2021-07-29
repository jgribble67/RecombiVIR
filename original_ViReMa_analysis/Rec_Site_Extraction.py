#!/bin/python3
##Last Modifed Feb19 by ALR
from subprocess import check_output
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Input", help="Input BED file with unclustered PASs e.g. hg19_PACs.bed. Must be sorted by count e.g. '$ sort -k4 -rn In.bed > In.sorted.bed'")
parser.add_argument("Genome", help="Genome_Path fasta")
parser.add_argument("Output", help="Unmasked Output BED file for clustered annotated PASs")
parser.add_argument("--Window", help="Nuc Window, default = 10")
args = parser.parse_args()

InFile = str(args.Input)
Genome = str(args.Genome)
if args.Window:
    Window = int(args.Window)
else:
    Window = 10

###################
def Rev_Comp(Seq):
        Seq = Seq.upper()
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        letters = list(Seq)
        letters = [basecomplement[base] for base in letters]
        return ''.join(letters)[::-1]

Output = open(str(args.Output), 'w')
with open(InFile, 'r') as In:
        line = In.readline().rstrip()
        while line:
                Data = line.split('\t')
                FromCoord = int(Data[1])
                ToCoord = int(Data[2])
                Strand = Data[5]
                Fromcmd= Data[0] + ":" + str(FromCoord - Window) + "-" + str(FromCoord + Window)
                FromSeq = check_output(['samtools', 'faidx', '-n', '1000', Genome, Fromcmd], universal_newlines=True).split()[1]
                if FromSeq:
                        FromSeq = FromSeq.upper()
                        if Strand == "-":
                                FromSeq = Rev_Comp(FromSeq)
                        else:
                                pass
                else:
                        print("Failed locus in index: ", Fromcmd)
                Tocmd= Data[0] + ":" + str(ToCoord - Window) + "-" + str(ToCoord + Window)
                ToSeq = check_output(['samtools', 'faidx', '-n', '1000', Genome, Tocmd], universal_newlines=True).split()[1]
                if ToSeq:
                        ToSeq = ToSeq.upper()
                        if Strand == "-":
                                ToSeq = Rev_Comp(ToSeq)
                        else:
                                pass
                else:
                        print("Failed locus in index: ", ToSeq)
                Output.write(line + '\t' + FromSeq + '\t' + ToSeq + '\n')
                line = In.readline().rstrip()
Output.close()


