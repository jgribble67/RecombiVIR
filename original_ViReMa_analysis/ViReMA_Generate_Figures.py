import ConfigViReMa as cfg
import matplotlib.pyplot as plt
import numpy as np
import math
import argparse
from pylab import *
from gooey import Gooey, GooeyParser

##      -------------------------------------------------------------------------------------------
##      Take arguments from command line, and send them to the config file for cross-module access
##      -------------------------------------------------------------------------------------------

@Gooey(program_name='GUI for ViReMa: A Virus Recombination Mapper',
       program_description = 'Make figure for ViReMa outputs',
      #return_to_config = True'
       )
def MainArgs():
        parser = GooeyParser()
        #Output_group = parser.add_argument_group("Output Handling", "Specify where output files are put and how they are named")

        ##Required Options        
        parser.add_argument("Input_Data", widget='FileChooser', help= "Enter name of sorted Recombination Data File (e.g. Virus_Recombination_Results): ")
        parser.add_argument("Output", help= "Name of Output  files")
        parser.add_argument("GeneLength", help= "Enter Gene nucleotide length")
        parser.add_argument("SaveShow", choices=['Show', 'Save'], help="Select data format.")       

        ##Optional
        parser.add_argument("--Output_Dir", widget='DirChooser', help= "Enter a directory name that all compiled output files will be saved in.")
        parser.add_argument("--MaxLimit", help= "Enter Gene nucleotide length")        
        parser.add_argument("--MovingAverage", help= "Enter MovingAverage for donor-acceptor plot")
        parser.add_argument("--Library", help= "Enter Library to Plot")
        parser.add_argument("--AveCoverage", help= "Enter average coverage over gene. \nN.B. Will only make sense if one library is chosen with --Library.")

        
        args = parser.parse_args()
        if args.SaveShow == 'Save':
            cfg.Save = 'SAVE'
        else:
            cfg.Save = 'SHOW'
        cfg.File1 = str(args.Input_Data)
        cfg.File3 = str(args.Output)
        cfg.GeneLength = int(args.GeneLength)
        if args.AveCoverage:
            cfg.AveCoverage = float(args.AveCoverage)
        else:
            cfg.AveCoverage = 0
        if args.MaxLimit:
            cfg.MaxLimit = int(args.MaxLimit)
        else:
            cfg.MaxLimit = 1000000000
        if args.MovingAverage:
            cfg.MVave = int(args.MovingAverage)
        else:
            cfg.MVave = 1
        if args.Library:
            cfg.Library = str(args.Library)
        else:
            cfg.Library = ''
        if args.Output_Dir:
            cfg.Output_Dir = str(args.Output_Dir) + '/'
        else:
            cfg.Output_Dir = ''

def MovingAverage(Sums, AverageRange):
        Averages = []
        for i in range(AverageRange/2):
                    Averages.append(0)
        n = 0
        for i in range(len(Sums) - AverageRange):
                Local = sum(Sums[n:AverageRange + n])
                Averages.append(Local)
                n+=1
        return Averages
           
def ExtractData(Data, Name):  
    Array = np.array([0] * (cfg.GeneLength*cfg.GeneLength))
    Array.shape = (cfg.GeneLength,cfg.GeneLength)
    AcceptorSums = []
    DonorSums = []
    for i in range(cfg.GeneLength):
        AcceptorSums.append(0.1)
        DonorSums.append(0.1)
    Gaps = [0] * (cfg.GeneLength + 1) * 2
    DelCover = [0]*cfg.GeneLength
    InsCover = [0]*cfg.GeneLength
        
    for line in Data:
            line = line.split("_")
            Donorsite = int(line[0])
            Acceptorsite = int(line[2])
            Count = int(line[4])
            Gap = Acceptorsite - Donorsite + 1
            Gaps[Gap + cfg.GeneLength] += Count
            if Count > 1:
                Array[Donorsite,Acceptorsite] = Count
                DonorSums[Donorsite] += Count
                AcceptorSums[Acceptorsite] += Count
                if math.fabs(Donorsite - Acceptorsite) < cfg.MaxLimit:
                        if Donorsite < Acceptorsite:
                                while Donorsite < Acceptorsite:
                                        DelCover[Donorsite] += Count
                                        Donorsite +=1
                        elif Acceptorsite < Donorsite:
                                while Acceptorsite < Donorsite:
                                        InsCover[Acceptorsite] += Count
                                        Acceptorsite +=1

    ##Make Figures
    #Heatmap(Array, Name)    
#    PlotGaps(Name)
    #DonorAcceptor(DonorSums, AcceptorSums, Name)
    Conservation(DelCover, InsCover, Name)

def Heatmap(Array, Name):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cmap = cm.jet
    cmap.set_under('w')
    ax.imshow(Array, cmap = cmap, interpolation = 'none', vmin = 1, vmax = 1180, origin = 'lower')
    #cbar = fig.colorbar(cax)
    if cfg.Save == "SHOW":
        plt.show()
    else:
        plt.savefig(cfg.Output_Dir + cfg.File3 + "_" + str(Name) + "_HeatMap.pdf", bbox_inches = 0)

#def PlotGaps(Name):
#    fig = plt.figure()
#    N = np.arange((cfg.GeneLength + 1) * 2)
#    ax = fig.add_subplot(111)
#    plt.semilogy()
#    ax.bar(N + 0.6,Gaps, color = 'black')
#    ##plt.axvspan(-5.6,5.6, facecolor = '0.5', alpha = 0.3)
#    ##ax.set_xbound(-40,40)
#    ##ax.set_ybound(0.1, 100000)
#    plt.show()

def DonorAcceptor(DonorSums, AcceptorSums, Name):
    DonorAve = MovingAverage(DonorSums, cfg.MVave)
    AcceptorAve = MovingAverage(AcceptorSums, cfg.MVave)
    width = 0.5
    fig = plt.figure()
    ax1 = fig.add_subplot(2,1,1)
    ##plt.semilogy()
    ax1.bar(N, DonorSums, width)
    ax1.plot(DonorAve)
    ax1.set_ybound(lower = 1, upper = 1000)
    ax2 = fig.add_subplot(2,1,2)
    ##plt.semilogy()
    ax2.bar(N, AcceptorSums, width, color = 'r')
    ax2.plot(AcceptorAve)
    ax2.set_ybound(lower = 1, upper = 1000)
    if cfg.Save == "SHOW":
        plt.show()
    else:
        plt.savefig(cfg.Output_Dir + cfg.File3 + "_" + str(Name) + "_Donor-Acceptor_Histogram.pdf", bbox_inches = 0)

def Conservation(DelCover, InsCover, Name):
    plt.subplot(211)
    if cfg.AveCoverage:
        DelCover = np.array(DelCover)/cfg.AveCoverage
    else:
        pass
    plt.plot(N, DelCover, 'b')
    plt.xlim(0,cfg.GeneLength)
    #plt.ylim(0, yaxish)
    plt.fill(DelCover, 'b')
    plt.grid(True)
    plt.subplot(212)
    plt.plot(N, InsCover, 'r')
    plt.xlim(0,cfg.GeneLength)
    plt.ylim(0)
    plt.fill(InsCover, 'r')
    plt.grid(True)
    if cfg.Save == "SHOW":
        plt.show()
    else:
        plt.savefig(cfg.Output_Dir + cfg.File3 + "_" + str(Name) + "Conservation_histogram.pdf", bbox_inches = 0)
    
    #plt.rcParams["figure.figsize"] = 5,2
    x = N
    y = np.array(DelCover)
    z = np.array(InsCover)
    fig = plt.figure()
    ax = fig.add_subplot(111)

#    fig, (ax2,ax) = plt.subplots(nrows=2)    
    extent = [x[0]-(x[1]-x[0])/2., x[-1]+(x[1]-x[0])/2.,0,1]
    cax = ax.imshow(y[np.newaxis,:], cmap="plasma", aspect="auto", extent=extent, vmin = 0.0, vmax = 1.0)
    cbar = fig.colorbar(cax)
    ax.set_yticks([])
    ax.set_xlim(extent[0], extent[1])
    
#    cax2 = ax2.imshow(z[np.newaxis,:], cmap="plasma", aspect="auto", extent=extent)
#    #cbar = fig.colorbar(cax2)
#    ax2.set_yticks([])
#    ax2.set_xlim(extent[0], extent[1])

    if cfg.Save == "SHOW":
        plt.show()
    else:
        plt.savefig(cfg.Output_Dir + cfg.File3 + "_" + str(Name) + "_ConservationBar.pdf", bbox_inches = 0)

if __name__ == '__main__':
    parser = GooeyParser()
    MainArgs()
    Dicts = {}
    N = np.arange(cfg.GeneLength)

    with open(cfg.File1,'r') as In:
            Name = In.readline()[13:-1]
            while Name:
                    Dicts[Name] = In.readline().split("\t")[:-1]
                    Finish =  In.readline()
                    Name = (In.readline())[13:-1]
    print "Available Recombination Libraries are:"
    for k in Dicts:
            print "Dict name:", k
         
    if cfg.Library:
        ExtractData(Dicts[str(cfg.Library)], str(cfg.Library))
    else:
        for k in Dicts:
            ExtractData(Dicts[k], k)
    
    

    
    