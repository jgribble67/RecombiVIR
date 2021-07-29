import argparse

parser = argparse.ArgumentParser()
parser.add_argument("Inputs", help="Meta data file.")
args = parser.parse_args()

Files = []
with open(str(args.Inputs), 'r') as In:
        header = In.readline()
        data = In.readline()
        while data:
            data = data.split()
            Files.append(data[0])
            data = In.readline()
        
Events = {}
for File in Files:
    with open(File + "/Virus_Recombination_Results.txt","r") as In:
    	Lib = In.readline()
    	while Lib:
    		Lib = Lib.split()[1]
    		if Lib in Events:
    			pass
    		else:
    			Events[Lib] = {}
    		Data = In.readline().split()
    		for i in Data:
    				i = i.split("_")
    				Event = '_'.join(i[:3])
    				if Event in Events[Lib]:
    					Events[Lib][Event][Files.index(File)] = i[4]
    				else:
    					Events[Lib][Event] = ['0'] * len(Files)
    					Events[Lib][Event][Files.index(File)] = i[4]
    		In.readline()
    		Lib = In.readline()

Output = open('Rec_Counts.txt','w')

Output.write('\t' + '\t'.join([i.split('/')[-1] for i in Files]) + '\n')

for Lib in Events:
    for i in Events[Lib]:
        Output.write(Lib + '_' + i + '\t' + '\t'.join([j for j in Events[Lib][i]]) + '\n')
Output.close()
















