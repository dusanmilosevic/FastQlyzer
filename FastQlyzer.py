"""
Usage: FastQlyzer.py [-h] [--qt=quality_threshold] FASTQ1 FASTQ2

Options:
    -h --help
    --qt Quality threshold between bad and good quality. [Default: 9]
"""
from docopt import docopt
import sys

arguments = docopt(__doc__)

docs=[]

filename1=arguments['FASTQ1']
filename2=arguments['FASTQ2']

docs.append(filename1)
docs.append(filename2)

file=open(docs[0],"r")
line=file.readline()
poorQualityBorder='*'

if arguments['--qt'] == None:
    poorQualityBorder = 9
elif arguments['--qt']>=0 and arguments['--qt']<=93:
    poorQualityBorder = int(arguments['--qt'])

poorQualityBorderValue=poorQualityBorder
numOfSeqflagedAsPoor=0
numOfGCchar=0;
totalChars=0;


#On line bellow we have first sequence
line=file.readline()

#Reading the number of Bases per line
numOfBasesPerRead= len(line)

file.close()


#Now that we know number of bases per line we can start reading
#.fastq file again

qualDict={0:0,1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0,13:0,14:0,15:0,16:0,17:0,18:0,19:0,20:0,21:0,22:0,23:0,24:0,25:0,26:0,27:0,28:0,29:0,30:0,31:0,32:0,33:0,34:0,35:0,36:0,37:0,38:0,39:0,40:0,41:0,42:0,43:0,44:0,45:0,46:0,47:0,48:0,49:0,50:0,51:0,52:0,53:0,54:0,55:0,56:0,57:0,58:0,59:0,60:0,61:0,62:0,63:0,64:0,65:0,66:0,67:0,68:0,69:0,70:0,71:0,72:0,73:0,74:0,75:0,76:0,77:0,78:0,79:0,80:0,81:0,82:0,83:0,84:0,85:0,86:0,87:0,88:0,89:0,90:0,91:0,92:0,93:0}
#Matrix we use to calculate median quality value per position
weightMatrix = [{0:0,1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:0,12:0,13:0,14:0,15:0,16:0,17:0,18:0,19:0,20:0,21:0,22:0,23:0,24:0,25:0,26:0,27:0,28:0,29:0,30:0,31:0,32:0,33:0,34:0,35:0,36:0,37:0,38:0,39:0,40:0,41:0,42:0,43:0,44:0,45:0,46:0,47:0,48:0,49:0,50:0,51:0,52:0,53:0,54:0,55:0,56:0,57:0,58:0,59:0,60:0,61:0,62:0,63:0,64:0,65:0,66:0,67:0,68:0,69:0,70:0,71:0,72:0,73:0,74:0,75:0,76:0,77:0,78:0,79:0,80:0,81:0,82:0,83:0,84:0,85:0,86:0,87:0,88:0,89:0,90:0,91:0,92:0,93:0} for k in range(numOfBasesPerRead)]


#Initializing lists that are needed
numOfAperRead=[ float(i*0) for i in range(numOfBasesPerRead) ]
numOfTperRead=[ float(i*0) for i in range(numOfBasesPerRead) ]
numOfCperRead=[ float(i*0) for i in range(numOfBasesPerRead) ]
numOfGperRead=[ float(i*0) for i in range(numOfBasesPerRead) ]
sumOfQualRead=[ float(i*0) for i in range(numOfBasesPerRead) ]


numOfSequences=0
numOfQualities=0

GCcontentPerPositionInRead=[ float(i*0) for i in range(numOfBasesPerRead) ]
GCdensity=[round(i*0.1,1) for i in range(1001)]
GChowMany=[i*0 for i in range(len(GCdensity))]

#Density represent quality- xlabel
#howMany represents howmanyQualities we have for each average value
densityOfQualities=[]
howMany=[]
i=0.0


while i<95:

	densityOfQualities.append(round(i,0))
	howMany.append(0)
	i+=1


#My program works as a slot machine with four states

order=0

for doc in docs:
	file=open(doc,"r")

	for line in file:

		if order==0:
			#Here we just read the sequence id
			#We don't need to analyze this so we just skip this
			order=1

		elif order==1:

			#READING SEQUENCE
			numOfSequences+=1

			numOfGCperSequence=0
			#per base sequence content
			for index in range(numOfBasesPerRead):

				if line[index]=='A':
					totalChars+=1
					numOfAperRead[index]+=1

				elif line[index]=='T':
					totalChars+=1
					numOfTperRead[index]+=1

				elif line[index]=='C':
					totalChars+=1
					numOfGCchar+=1
					numOfCperRead[index]+=1
					GCcontentPerPositionInRead[index]+=1
					numOfGCperSequence+=1

				elif line[index]=='G':
					totalChars+=1
					numOfGCchar+=1
					numOfGperRead[index]+=1
					GCcontentPerPositionInRead[index]+=1
					numOfGCperSequence+=1

			average=numOfGCperSequence/float(numOfBasesPerRead)
			average*=100.0
			indeksic=GCdensity.index(round(average,1))
			GChowMany[indeksic]+=1

			order=2

		elif order==2:

			#Here we read "+" char so we just skip this
			order=3

		elif order==3:

			#READING QUALITIES
			#Qualities in lower-higher order
			#!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~

			numOfQualities+=1
			seqQualitySum=0.0;

			for index in range(numOfBasesPerRead):
				if line[index]!='\n':
					quality=ord(line[index])
					quality-=33
					seqQualitySum+=quality
					sumOfQualRead[index]+=quality
					if quality in weightMatrix[index].keys():
						weightMatrix[index][quality]+=1

			seqQualitySum/=float(numOfBasesPerRead);
			seqQualitySum=round(seqQualitySum,0);
			middleValue=seqQualitySum
			if middleValue<poorQualityBorderValue:
				numOfSeqflagedAsPoor+=1
			ind= densityOfQualities.index(middleValue)
			howMany[ind]+=1
			order=0
	file.close()





for i in range(numOfBasesPerRead):

	numOfAperRead[i]/=float(numOfSequences)
	numOfAperRead[i]*=100

	numOfTperRead[i]/=float(numOfSequences)
	numOfTperRead[i]*=100

	numOfCperRead[i]/=float(numOfSequences)
	numOfCperRead[i]*=100

	numOfGperRead[i]/=float(numOfSequences)
	numOfGperRead[i]*=100

	GCcontentPerPositionInRead[i]/=float(numOfSequences)
	GCcontentPerPositionInRead[i]*=100

	sumOfQualRead[i]/=float(numOfQualities)

###########################################################################################
##########################Mean quality calculation#########################################
def median(lst):
    lst = sorted(lst)
    if len(lst) < 1:
            return None
    if len(lst) %2 == 1:
            return lst[((len(lst)+1)/2)-1]
    else:
            return float(sum(lst[(len(lst)/2)-1:(len(lst)/2)+1]))/2.0

medians=[0.0 for i in range(numOfBasesPerRead)]
leftIndex=0
rightIndex=0
medianIndex=0
if numOfQualities%2==1:
	medianIndex=numOfQualities/2

	for i in range(len(weightMatrix)):
		sum=0
		kljuc=0
		indicator=False
		for key in weightMatrix[i].keys():
			kljuc=key
			sum+=weightMatrix[i][key]
			if sum>=medianIndex:
				indicator=True
				break
		if indicator:
			medians[i]=kljuc

else:
	rightIndex=numOfQualities/2
	leftIndex=rightIndex-1
	for i in range(len(weightMatrix)):
		sum=0
		kljuc1=0
		kljuc2=0
		indicator=False
		for key in weightMatrix[i].keys():
			kljuc1=key
			sum+=weightMatrix[i][key]
			if sum>=leftIndex:
				indicator1=True
				break
		sum=0
		for key in weightMatrix[i].keys():
			kljuc2=key
			sum+=weightMatrix[i][key]
			if sum>=rightIndex:
				indicator2=True
				break
		if indicator1 and indicator2:
			medians[i]=(kljuc1+kljuc2)/2.0


#############################################################################################
##############################################################################################

import numpy as np
import matplotlib.pyplot as plt

#Library for legend
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from numpy import arange

x=[i for i in range(numOfBasesPerRead)]
#########################################################################################
#########################################################################################

fig0 = plt.figure(0)

plt.title('Basic Statistics')
plt.axis([0, 15, 0, 15])
axes=plt.gca()
axes.xaxis.set_major_formatter(plt.NullFormatter())
axes.yaxis.set_major_formatter(plt.NullFormatter())
plt.text(1, 13, 'Filename1:', ha='left', wrap=True,color='blue')
plt.text(5, 13, filename1, ha='left', wrap=True)
plt.text(1, 11, 'Filename2:', ha='left', wrap=True,color='blue')
plt.text(5, 11, filename2, ha='left', wrap=True)
plt.text(1, 9, 'Poor Quality Border:', ha='left', wrap=True,color='blue')
plt.text(7, 9, str(poorQualityBorder), ha='left', wrap=True)
plt.text(1, 7, 'Total Sequences', ha='left', wrap=True,color='blue')
plt.text(6, 7, str(numOfSequences), ha='left', wrap=True)
g = "Seq. Flaged as"\
	" poor Quality:"
plt.text(1, 5, g, ha='left', wrap=True,color='blue')
plt.text(8, 5, str(numOfSeqflagedAsPoor), ha='left', wrap=True)
plt.text(1, 3, 'Sequence Length:', ha='left', wrap=True, color='blue')
plt.text(6, 3, str(numOfBasesPerRead), ha='left', wrap=True)
plt.text(1, 1, '%GC:', ha='left', wrap=True,color='blue')
plt.text(3, 1, str(round(numOfGCchar/float(totalChars)*100,1)), ha='left', wrap=True)
#########################################################################################
#########################################################################################

#AVERAGE QUALIY PER POSITION
fig1=plt.figure(1)
axes=plt.gca()
axes.set_ylim([0,60])
axes.set_xlim([0,numOfBasesPerRead-2])
plt.xlabel('Position in read(bp)')
plt.ylabel('Quality')
plt.title('Quality scores across all bases')
plt.grid(True)
blue_patch = mpatches.Patch(color='blue', label='Average(Mean) Quality')
red_patch = mpatches.Patch(color='red', label='Median Quality')
plt.legend(handles=[blue_patch,red_patch])
#plt.xticks(arange(numOfBasesPerRead-2, step=5))
plt.plot(x,sumOfQualRead,'b',x,medians,'r')

#########################################################################################
#########################################################################################

#Per sequence quality score
#Quality score distribution over all sequences
fig2=plt.figure(2)
axes=plt.gca()
axes.set_xlim([0,50])
plt.xlabel('Mean Sequence Quality(Phred Score)')
#plt.ylabel('Proportion[%]')
plt.title('Quality score distribution over all sequences')
plt.grid(True)
plt.plot(densityOfQualities,howMany,'r')

#########################################################################################
#########################################################################################

#Per base sequence content
#Sequence content across all bases
fig3=plt.figure(3)

blue_patch = mpatches.Patch(color='blue', label='%A')
red_patch = mpatches.Patch(color='red', label='%G')
green_patch = mpatches.Patch(color='green', label='%T')
yellow_patch = mpatches.Patch(color='yellow', label='%C')

plt.xlabel('Position in read(bp)')
plt.ylabel('Proportion[%]')
plt.title('Sequence content across all bases')

plt.legend(handles=[blue_patch,red_patch,green_patch,yellow_patch])
plt.grid(True)

axes=plt.gca()
axes.set_xlim([0,numOfBasesPerRead-2])
axes.set_ylim([0,100])

plt.plot(x,numOfAperRead,'b',x,numOfTperRead,'g',x,numOfCperRead,'y',x,numOfGperRead,'r',lw=2)

################################################################################################
################################################################################################

#Per base GC content
#GC content across all bases
fig4=plt.figure(4)
plt.xlabel('Position in read(bp)')
plt.ylabel('Proportion[%]')
plt.title('GC content across all bases')
plt.grid(True)
axes=plt.gca()
axes.set_xlim([0,numOfBasesPerRead-2])
axes.set_ylim([0,100])
#axes.xaxis.set_major_formatter(plt.NullFormatter())
plt.plot(x,GCcontentPerPositionInRead,'g')

################################################################################################
################################################################################################

#Per Sequence GC Content
#GC distribution over all sequences
fig5=plt.figure(5)
plt.xlabel('Mean GC content (%)')
#plt.ylabel('Proportion[%]')
plt.title('GC distribution over all sequences')
plt.grid(True)
red_patch = mpatches.Patch(color='red', label='GC count per read')
plt.legend(handles=[red_patch])
plt.plot(GCdensity,GChowMany,'r')


#Testing time!
#plt.show()

#Saving graphs into a single pdf file
import matplotlib.backends.backend_pdf
pdf = matplotlib.backends.backend_pdf.PdfPages("report.pdf")

pdf.savefig(fig0)
pdf.savefig(fig1)
pdf.savefig(fig2)
pdf.savefig(fig3)
pdf.savefig(fig4)
pdf.savefig(fig5)

pdf.close()
