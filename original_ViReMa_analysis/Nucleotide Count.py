import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


##DI = open("Virus_Recombination_Results_0uM.txt", "r")
Seq=open("JX869059.fasta","r")#MERS
#Seq=open("MT020881.fasta","r")#SARS
#Seq=open("SHC014_CoV.fasta","r")#SHC
Seq.readline()
frag=[]
for line in Seq:
    f = line
    frag.append(f[:-1])
Genome = "".join(frag)
Genome = Genome.upper()
Seq.close

#Genome = Genome[:-75]
As = (Genome.count('A')/len(Genome)) * 100
Gs = (Genome.count('G')/len(Genome)) * 100
Ts = (Genome.count('T')/len(Genome)) * 100
Cs = (Genome.count('C')/len(Genome)) * 100

#Window
N = 40

class NTs(object):
    def __init__(self, Name):
        self.DcountA = np.array([0]*N)
        self.DcountG = np.array([0]*N)
        self.DcountC = np.array([0]*N)
        self.DcountT = np.array([0]*N)
        self.AcountA = np.array([0]*N)
        self.AcountG = np.array([0]*N)
        self.AcountC = np.array([0]*N)
        self.AcountT = np.array([0]*N)
        self.total = 0

    def Finish(self):
        self.DFreqA = (self.DcountA/self.total)*100
        self.DFreqG = (self.DcountG/self.total)*100
        self.DFreqC = (self.DcountC/self.total)*100
        self.DFreqT = (self.DcountT/self.total)*100
        self.AFreqA = (self.AcountA/self.total)*100
        self.AFreqG = (self.AcountG/self.total)*100
        self.AFreqC = (self.AcountC/self.total)*100
        self.AFreqT = (self.AcountT/self.total)*100

Data = {}
#Files = '''WT-1_DVGs.txt
# WT16hmA_DVGs.txt
# WT16hmB_DVGs.txt
# WT16hmC_DVGs.txt
# WT-2_DVGs.txt
# WT-3_DVGs.txt
# XN-1_DVGs.txt
# XN-2_DVGs.txt
# XN24hmA_DVGs.txt
# XN24hmB_DVGs.txt
# XN24hmC_DVGs.txt
# XN-3_DVGs.txt
# '''
Files = '''
MERS_LCP1_forward_junctions.txt
MERS_LDP1_forward_junctions.txt
MERS_LEP1_forward_junctions.txt
SARSCoV2_A_forward_junctions.txt
SARSCoV2_B_forward_junctions.txt
SARSCoV2_C_forward_junctions.txt
'''.split()



for File in Files:
    Data[File] = NTs(File)
    with open(File, "r") as DI:
        data = DI.readline()
        #data = DI.readline()
        data = DI.readline()
        while data:
            data = data.split()
            counts = 1
            #counts = int(data[3])
            donorsite = int(data[1])
            acceptorsite = int(data[2])
            Gap = acceptorsite - donorsite - 1
               #if math.fabs(Gap) > 10 and math.fabs(Gap) < 500:
                  #if counts > 3:
            donorsequence = Genome[int(donorsite) - 20:int(donorsite +20)]
            acceptorsequence = Genome[int(acceptorsite) - 21:int(acceptorsite +19)]
            if len(acceptorsequence) >= N and len(donorsequence) >= N:# and donorsite > 100:
                Data[File].total = Data[File].total + counts
                for i in range(N):
                    if acceptorsequence[i] == 'A':
                        Data[File].AcountA[i] = Data[File].AcountA[i] + counts
                    if acceptorsequence[i] == 'G':
                        Data[File].AcountG[i] = Data[File].AcountG[i] + counts
                    if acceptorsequence[i] == 'C':
                        Data[File].AcountC[i] = Data[File].AcountC[i] + counts
                    if acceptorsequence[i] == 'T':
                        Data[File].AcountT[i] = Data[File].AcountT[i] + counts
                for i in range(N):
                    if donorsequence[i] == 'A':
                        Data[File].DcountA[i] = Data[File].DcountA[i] + counts
                    if donorsequence[i] == 'G':
                        Data[File].DcountG[i] = Data[File].DcountG[i] + counts
                    if donorsequence[i] == 'C':
                        Data[File].DcountC[i] = Data[File].DcountC[i] + counts
                    if donorsequence[i] == 'T':
                        Data[File].DcountT[i] = Data[File].DcountT[i] + counts
            data = DI.readline()
        Data[File].Finish()

ticks = np.concatenate((np.arange(1,20,2), np.arange(20,39,2)))         
fusionsites = np.arange(-19,20,2)
ind = np.arange(N)
indc = np.arange(20)
width = 1
fig = plt.figure()

ax1 = fig.add_subplot(2,2,1)
ax1.annotate('A', (2,90))
plt.axvline(x = 19.5, linewidth = 3 )
plt.axhspan(0,As, facecolor = '0.5', alpha = 0.3 )

ax2 = fig.add_subplot(2,2,2)
ax2.annotate('T', (2,90))
plt.axvline(x = 19.5, linewidth = 3 )
plt.axhspan(0,Ts, facecolor = '0.5', alpha = 0.3 )

ax3 = fig.add_subplot(2,2,3)
plt.axvline(x = 19.5, linewidth = 3 )
plt.axhspan(0,Cs, facecolor = '0.5', alpha = 0.3 )
ax3.annotate('C', (2,90))

ax4 = fig.add_subplot(2,2,4)
ax4.annotate('G', (2,90))
plt.axvline(x = 19.5, linewidth = 3 )
plt.axhspan(0,Gs, facecolor = '0.5', alpha = 0.3 )

WtFilesDVGs = '''
 WT-1_DVGs.txt
 WT-2_DVGs.txt
 WT-3_DVGs.txt
'''.split()

XNFilesDVGs = '''
 XN-1_DVGs.txt
 XN-2_DVGs.txt
 XN-3_DVGs.txt
'''.split()

WtTCFilesDVGs = '''
 WT16hmA_DVGs.txt
 WT16hmB_DVGs.txt
 WT16hmC_DVGs.txt
'''.split()

XNTCFilesDVGs = '''
 XN24hmA_DVGs.txt
 XN24hmB_DVGs.txt
 XN24hmC_DVGs.txt
'''.split()

SARS2 = '''
SARSCoV2_A_forward_junctions.txt
SARSCoV2_B_forward_junctions.txt
SARSCoV2_C_forward_junctions.txt
'''.split()

MERS = '''
MERS_LCP1_forward_junctions.txt
MERS_LDP1_forward_junctions.txt
MERS_LEP1_forward_junctions.txt
'''.split()

SHC = '''
SHC_P10_1-1_Virus_Recombination_Results.bed
SHC_P10_1-2_Virus_Recombination_Results.bed
SHC_P10_1-3_Virus_Recombination_Results.bed
'''.split()

for i in MERS:
    rects1 = ax1.plot(ind, Data[i].DFreqA, width, color='r',)
    rects2 = ax2.plot(ind, Data[i].DFreqT, width, color='r',)
    rects3 = ax3.plot(ind, Data[i].DFreqC, width, color='r',)
    rects4 = ax4.plot(ind, Data[i].DFreqG, width, color='r',)
#for i in SARS2:
#    rects1 = ax1.plot(ind, Data[i].DFreqA, width, color='b',)
#    rects2 = ax2.plot(ind, Data[i].DFreqT, width, color='b',)
#    rects3 = ax3.plot(ind, Data[i].DFreqC, width, color='b',)
#    rects4 = ax4.plot(ind, Data[i].DFreqG, width, color='b',)

#ax.annotate('random', (52, 20))
#ax.annotate('non-random', (52, 27))
#ax1.annotate(total, (52, 90))
#ax.annotate('N = ', (47, 90))
ax1.set_xticks((ticks))
ax1.set_xticklabels(fusionsites)
ax2.set_xticks((ticks))
ax2.set_xticklabels(fusionsites)
ax3.set_xticks((ticks))
ax3.set_xticklabels(fusionsites)
ax4.set_xticks((ticks))
ax4.set_xticklabels(fusionsites)
ax1.set_yticks((0,25,50,75,100))
ax2.set_yticks((0,25,50,75,100))
ax3.set_yticks((0,25,50,75,100))
ax4.set_yticks((0,25,50,75,100))
ax1.set_ylabel('Percentage Base Identity')
ax3.set_ylabel('Percentage Base Identity')
ax1.set_title('Base Position Relative to Donor Site')
ax2.set_title('Base Position Relative to Donor Site')
plt.show()


fig1 = plt.figure()
ax1 = fig1.add_subplot(2,2,1)
ax1.annotate('A', (2,90))
plt.axvline(x = 19.5, linewidth = 3 )
plt.axhspan(0,As, facecolor = '0.5', alpha = 0.3 )

ax2 = fig1.add_subplot(2,2,2)
ax2.annotate('T', (2,90))
plt.axvline(x = 19.5, linewidth = 3 )
plt.axhspan(0,Ts, facecolor = '0.5', alpha = 0.3 )

ax3 = fig1.add_subplot(2,2,3)
ax3.annotate('C', (2,90))
plt.axvline(x = 19.5, linewidth = 3 )
plt.axhspan(0,Cs, facecolor = '0.5', alpha = 0.3 )

ax4 = fig1.add_subplot(2,2,4)
ax4.annotate('G', (2,90))
plt.axvline(x = 19.5, linewidth = 3 )
plt.axhspan(0,Gs, facecolor = '0.5', alpha = 0.3 )

for i in MERS:
    rects1 = ax1.plot(ind, Data[i].AFreqA, width, color='r',)
    rects2 = ax2.plot(ind, Data[i].AFreqT, width, color='r',)
    rects3 = ax3.plot(ind, Data[i].AFreqC, width, color='r',)
    rects4 = ax4.plot(ind, Data[i].AFreqG, width, color='r',)
#for i in SARS2:
#    rects1 = ax1.plot(ind, Data[i].AFreqA, width, color='b',)
#    rects2 = ax2.plot(ind, Data[i].AFreqT, width, color='b',)
#    rects3 = ax3.plot(ind, Data[i].AFreqC, width, color='b',)
#    rects4 = ax4.plot(ind, Data[i].AFreqG, width, color='b',)
#ax.annotate('random', (52, 20))
#ax.annotate('non-random', (52, 27))
#ax1.annotate(total, (52, 90))
#ax.annotate('N = ', (47, 90))
ax1.set_xticks((ticks))
ax1.set_xticklabels(fusionsites)
ax2.set_xticks((ticks))
ax2.set_xticklabels(fusionsites)
ax3.set_xticks((ticks))
ax3.set_xticklabels(fusionsites)
ax4.set_xticks((ticks))
ax4.set_xticklabels(fusionsites)
ax1.set_yticks((0,25,50,75,100))
ax2.set_yticks((0,25,50,75,100))
ax3.set_yticks((0,25,50,75,100))
ax4.set_yticks((0,25,50,75,100))
ax1.set_ylabel('Percentage Base Identity')
ax3.set_ylabel('Percentage Base Identity')
ax1.set_title('Base Position Relative to Acceptor Site')
ax2.set_title('Base Position Relative to Acceptor Site')
plt.show()




            
    
