import os, subprocess, time
from os.path import isfile, join
import shutil # to move files from af folder to another
import math

def GetListFile(PathFile, FileExtension):
    return [os.path.splitext(f)[0] for f in os.listdir(PathFile) if isfile(join(PathFile, f)) and os.path.splitext(f)[1] == '.' + FileExtension]


def openfile(Path):
    fileIn = open(Path, "r")
    lines = fileIn.readlines()
    fileIn.close()
    return lines


import re
def parse_Ct_Structure(Path):
    lines = openfile(Path)
    # Get the initial position of the read and it Series ||||| with mutations
    # replace one space in case of multi spaces re.sub( '\s+', ' ', mystring ).strip()
    #print Path
    #print [int(re.sub( '\s+', ' ', elem ).strip().split(' ')[4]) for elem in lines[1:]]
    return [int(re.sub( '\s+', ' ', elem ).strip().split(' ')[4]) for elem in lines[1:]]

def parse_SHAPE(Path):
    lines = openfile(Path)
    lista=[]
    # Get the initial position of the read and it Series ||||| with mutations
    # replace one space in case of multi spaces re.sub( '\s+', ' ', mystring ).strip()

    Intermediate=[re.sub( '\s+', '\t', elem ).strip().split('\t') for elem in lines]
    for elem in  Intermediate:
        if len(elem)>2 and float(elem[2])!=0:
            lista.append(float(elem[2])) 
        else:
            lista.append(-10)
    return lista
def Pairing_status(structure):# a line in ct format with a number for the partner if paired and 0 if not
    status=[]
    for value in structure:
        #print value,'lol'
        if value== '(' or value== ')':
            status.append('P')
	if value =='.':
	    status.append('Un')      
        if value=='x':
            status.append('PK')
    return status

def plot2D(x,y,titre):
    import matplotlib.pyplot as plt
    import numpy as np

    plt.plot(x,y,'r.')
    plt.xlabel('shape')
    plt.ylabel('unpaired probability')
    plt.title(titre)
    # fitting functions
    plt.show()
def distribution_plots(Y,titre):
    # plot the x distribution

    import matplotlib.pyplot as plt
    import numpy as np
    '''
    d = {x: Y.count(x) for x in Y}
    print d
    plt.plot(d.keys(),d.values(),'.')
    '''
    plt.hist(Y, bins=10, color='green')
    plt.ylabel('Frequencies')
    plt.title(titre)
    # fitting functions
    #plt.show()
    plt.savefig(titre)

def GetbasepairsProb(path_Fasta, fil, FileExtensionFasta,Path_dot_plot):

    listProbaPairs=[]
    SingleProba=[]
    FastaPath=os.path.join(path_Fasta, fil + '.' + FileExtensionFasta)
    rna =openfile(FastaPath)[1]
    #print rna
    os.system("RNAfold -p -d2 --noLP <" + FastaPath+ ">output.txt")
    PSPath=os.path.join( Path_dot_plot,fil+"_dp.ps")
    shutil.move(fil+"_dp.ps",PSPath)
    os.remove(fil+"_ss.ps")
    #print fil,'rr'
    bpm = loadDotPlotPS(PSPath)
    dp = DotPlot(rna, bpm)

    for i in range(len(rna)):
        for j in range(i, len(rna)):
            if dp.getBPProb(i, j) > 0:# get only non null probabilities
                listProbaPairs.append((i,j,dp.getBPProb(i, j)))

    SingleProba=dp.getUnpairedProbs()
    return listProbaPairs, SingleProba
#!!!!!!!!!!!!!!!!!!!!!!!!!!!! loadDotPlotPS(path)
def loadDotPlotPS(path):
    res = {}
    outTmp = open(path)
    for l in outTmp:
        data = l[:-1].split()
        if len(data) == 4 and data[3]=="ubox":
            i = int(data[0])-1
            j = int(data[1])-1
            p = math.pow(float(data[2]),2.)
            res[i,j] = p
    outTmp.close()
    return res

def parse_rdat(Path):
	lines = openfile(Path)
	RNA=[]
	seq=dict()
	struct=dict()
	reactivity=dict()
	for line in lines:
		#print line.split('\t')[0]
		if line.split('\t')[0]== 'NAME':
			RNA.append(line.split('\t')[1][:-1])
		        Id= line.split('\t')[1][:-1]
		        #print 'id',Id
		if line.split('\t')[0]=='SEQUENCE':
			seq[Id]= line.split('\t')[1][:-2]
		if line.split('\t')[0]=='STRUCTURE':
			struct[Id]=line.split('\t')[1][:-2]
		if line.split('\t')[0]=='REACTIVITY:1':
                        #print line.split('\t')[1:-1]
			reactivity[Id]=line.split('\t')[1:-1]
	return RNA,seq, struct,reactivity


def create_fasta_shapeFiles(RNA,seq, struct,reactivity,path_fasta,path_SHAPE):
    for Id in RNA:
	    Outfasta=open(os.path.join(path_fasta, Id+'.fa'),'w')
	    OutShape=open(os.path.join(path_SHAPE, Id+'.shape'),'w')
	    Outfasta.write("%s \n" % (">"+ Id))
	    Outfasta.write("%s \n" % (seq[Id]))
	    Outfasta.write("%s " % (struct[Id]))
	    Outfasta.close()
            #print Id, len(reactivity[Id]),len(seq[Id]) 
	    for i, val in enumerate(reactivity[Id][:-1]):
                #print reactivity[Id]
                if i <len(seq[Id])and val!=" ": 
                        print Id,i, seq[Id][i],"FF",val
			OutShape.write("%i \t %s \t %f \n"%(i+1,seq[Id][i],float(val)))
                        #print "done"
class DotPlot:
    """Class for holding/producing base-pair probability matrices"""
    def __init__(self, rna , bpm = None):
        self.rna = rna[:]
        if bpm is None:
            # we will avoid this case to be sure that rnafold from the min works well
            self.bpm = self.runRNAFold()
        else:
            self.bpm = bpm
    def getSeq(self):
        return self.rna

    def getBPProb(self,i,j):
        if (i,j) in self.bpm:
            return self.bpm[i,j]
        else:
            return 0.

    def getUnpairedProbs(self):
        res = [1. for i in self.rna]
        for i,j in self.bpm:
            res[i] -= self.bpm[i,j]
            res[j] -= self.bpm[i,j]
        return res

def Parsefile(Path):
    fileIn = open(Path, "r")
    lines = fileIn.readlines()
    fileIn.close()
    return lines


def parseReactivityfile(fileinput):
   Reactvities=[]
   lines=Parsefile(fileinput)
   for it in range(len(lines)):
                if (lines[it].split("\t")[2][:-1]):
                	Reactvities.append(lines[it].split("\t")[2][:-1])
		else:
			Reactvities.append(-10)
   return Reactvities



if __name__ == '__main__':
    ####################"" To parametrize ##########################
    FileExtensionshape ='shape'
    react='NMIA'
    path_SHAPE='SHAPE_files_NMIA'
    path_Fasta = 'fasta_files'
    FileExtensionFasta = 'fa'


    Shape={}
    States={}

    BP={}
    UNP={}
    Quaternary = {}
    tertiary = {}
    Unpaired2 = {}
    lisTIP = []
    lisTUn = []
    lisHE = []
    lisHES = []
    lisTIPinterne=[]
    lisTIPexterne=[]
    SPTIP = []
    SPTUn = []
    SPHE = []
    SPHES = []
    SPTIPinterne = []
    SPTIPexterne = []

    lisIP = []  # shape values for paired Positions
    lisUn = []  # shape values for unpaired Positions
    SPIP = []  # probability of being unpaired  from McCaskill for IP category
    SPUn = []  # probability of being unpaired  from McCaskill for Un category
    RNA=[]

    
    struct=dict()
    reactivity=dict()

       
    
	
    list3=[]
  
    for filz in GetListFile(path_Fasta, FileExtensionFasta):
    #print reactivity
  
        
        rna = Parsefile(os.path.join(path_Fasta, filz + '.' + FileExtensionFasta))[1]
        structure=Parsefile(os.path.join(path_Fasta, filz + '.' + FileExtensionFasta))[2]
        States[filz] =  Pairing_status(structure)
        reactivity[filz]=parseReactivityfile(os.path.join(path_SHAPE, filz + react+'Shape.txt'))

        # Get the end-Helix positions
        print "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",filz, len(States[filz])
        print States[filz][-1]
        for i,elem in enumerate(States[filz]):
            
            if elem=='Un' or elem=='PK':
                list3.append(elem)
            if elem=='P':
                #print i, elem,i-1,i+1, States[filz][i+1]
                if i in range(1,len(rna)-3) and (States[filz][i-1]=='Un'  or States[filz][i-1]=='PK'  or States[filz][i+1]=='Un' or States[filz][i-1]=='PK' ): # wh should add PK constraint because we are looking for stacked substructures that does not take into account thhe tertiary or pseudoknots extension!!
                    list3.append('HE')
                else:
                    list3.append(elem)
    cum=[]
    for filz in GetListFile(path_Fasta, FileExtensionFasta):
        if reactivity[filz]==[]:
		print "warning!! Empty reactivity",rna
	cum=cum+reactivity[filz]
    Shape=cum
    #print len(Shape)
    lIP = []  # shape values for paired Positions
    lUn = []
    lHE =[]
    # for structural_study 
    #print [elem[:-1] for elem in Shape]
    for nucl,shape in zip(list3,Shape):
		if shape!=-10  and nucl=='P':
                        print "la vaaleur", shape
		        lIP.append( float(shape))
		        #SPIP.append(UnpPb )
		if shape!=-10  and nucl=='Un':
		        lUn.append( float(shape))
		        #SPUn.append(UnpPb )
		if shape!=-10  and nucl=='HE':
		        lHE.append( float(shape))
    import numpy as np
    labels=["Stacked nucleotides","Unpaired nucleotides","Helix-end"]
    lists= [lIP,lUn, lHE]
    for (data, title) in [(lIP, "P" + str(react)), (lUn, "U"+ str(react)),(lHE,"HE"+ str(react))]:
        	print  title ,'\n' , data
