#uhomologyplot
import numpy as np
import pandas as pd
import os
import fnmatch

report = pd.DataFrame()

Dict = {}
N = 20
wd = "/Users/jennifergribble/Dropbox/P250_recombination/Genewiz_samples/P250_cell/Shannon_Entropy/"
od = "/Users/jennifergribble/Dropbox/P250_recombination/Genewiz_samples/P250_cell/"
exp = "P250_cell"

for file in os.listdir(wd):
    if fnmatch.fnmatch(file, "*_Virus_Recombination_Results.txt"):
        sample_name = str(file.split("_")[0])
        Dict[file] = np.array([0]*N)
        with open(wd + file, 'r') as In:
            Data = In.readline()
            while Data:
                if 'RevStrand' in Data:
                    Data = In.readline()
                    Data = In.readline()
                    Data = In.readline()
                else:
                    Data = In.readline()
                    Data = Data.split()
                    for j in Data:
                        Fuzz = int(len(j.split('_')[1][1:]))
                        # Count = int(j.split('_')[-1])
                        Count = 1
                        Dict[file][Fuzz] += Count
                    Data = In.readline()
                    Data = In.readline()
                    # print(i)
                    # print(Dict[i])
                    Dict[file] = Dict[file]/np.sum(Dict[file])
                    # print(Dict[i])
                    report[sample_name] = Dict[file].tolist()
report.to_csv(od + exp + "_homology.txt", sep="\t")

# def MakeTheoreticalDistribution(N):
#     Dist = [0] * (N)
#     Dist[0] = 1.0
#     for i in range(1,N):
#         Prob = 0.25**i
#         Dist[i] = Prob
#     for i in range(N)[::-1]:
#         Dist[i] -= sum(Dist[i+1:])
#     return Dist

#Dist = np.array(MakeTheoreticalDistribution(N))
#plt.bar(np.arange(N), Dist, width=0.2, color='b')
#x = 0
#for i in Files:
#    if 'MERS' in i:
#        x += 0.25
#        plt.bar(np.arange(N) + x, Dict[i], width=0.2, color='orange')
#plt.bar(np.arange(N), Dist, width=0.2, color='b')
#x = 0
#for i in Files:
#    if 'SARS' in i:
#        x += 0.25
#        plt.bar(np.arange(N) + x, Dict[i], width=0.2, color='red')
#plt.bar(np.arange(N), Dist, width=0.2, color='b')
#plt.bar(np.arange(N), Dist, width=0.2, color='b')
#x = 0
#for i in Files:
#    if 'WT' in i:
#        x += 0.25
#        plt.bar(np.arange(N) + x, Dict[i], width=0.2, color='green')
#plt.bar(np.arange(N), Dist, width=0.2, color='b')
#x = 0
#for i in Files:
#    if 'XN' in i:
#        x += 0.25
#        plt.bar(np.arange(N) + x, Dict[i], width=0.2, color='cyan')

