#ViReMa array transformations
import numpy as np
import pandas as pd

Files = '''/Users/jennifergribble/Dropbox/Perlman_MERS_SARS2/Vantage_RNAseq/raw_junction_files/0-01B_virema_Virus_Recombination_Results.txt'''.split()

Dict = {}

for i in Files:
    with open(i, 'r') as File1:
        Data = File1.readline()
        while Data:
            Name = Data[:-1]
            Dict[Name] = File1.readline().split("\t")[:-1]
            Data = File1.readline()
            Data = File1.readline()
            Data = File1.readline()
            print(Data)