import fileinput
import os
from multiprocessing import Pool
import pandas as pd
mapData=list()
#ancestries = ('African', 'European', 'Amerindian')
ancestry="African"
numbers = range(1, 23)
chromosomes = [number for number in numbers]
famFile = 'masked/masked'+ ancestry+ 'ChrAllAboveMinimum.fam'
famData = pd.read_table(famFile,index_col=1,header=None)
minSNPs = 7
maxSNPs = 20
maxRate = 0.3

def makeHaplotypeBlock(chromosome):
  haplotypes = []
  leftSNP = 0
  rightSNP = 0
  snps =  []
  nrows=len(mapData[chromosome-1].index)
  numbers = range(0,nrows)
  rowindices = [number for number in numbers]
  for i in rowindices:
   if rightSNP > (nrows-2):
     break
   rightSNP = rightSNP + 1
   if mapData[chromosome-1].iat[rightSNP, 1] - mapData[chromosome-1].iat[leftSNP, 1] >= maxRate:
    if rightSNP - leftSNP >= minSNPs:
     snps = mapData[chromosome-1].iloc[leftSNP:(rightSNP - 1), 2].tolist()
    else:
     leftSNP = rightSNP
     rightSNP = leftSNP + 1

   # If the haplotype is long enough
   if rightSNP - leftSNP >= maxSNPs:
    snps = mapData[chromosome-1].iloc[leftSNP: (rightSNP - 2), 2].tolist()

  # If SNPs isn't NULL, then find the save the range of the window
   if (snps):
#    haplotype = haplotype + 1
    mylines=[chromosome-1,leftSNP,rightSNP,mapData[chromosome-1].iat[rightSNP - 2,1] - mapData[chromosome-1].iat[leftSNP,1]]
    haplotypes.append(mylines)
    snps=[]
    leftSNP = rightSNP

  ancestrychrfilename=ancestry+"."+str(chromosome)+".tsv"
  pd.DataFrame(haplotypes).to_csv(ancestrychrfilename,sep='     ', header=False, index=False)



for chromosome in chromosomes:
  mapFile = 'masked/rfMixInChr'+ str(chromosome)+'_chr'+str(chromosome)+'.map'
  mapData.append(pd.read_table(mapFile,header=None))

#spool to chromosomes and ancestries to the makeHaplotypeBlock function
agents = 22
pool=Pool(processes=agents)
result = pool.map(makeHaplotypeBlock,chromosomes)

#concat the output files
with open(ancestry+".tsv", "wb") as outfile:
    for chromosome in chromosomes:
        with open(ancestry+"."+str(chromosome)+".tsv", "rb") as infile:
            outfile.write(infile.read())

for chromosome in chromosomes:
  try:
    filepath=ancestry+"."+str(chromosome)+".tsv"
    os.remove(filepath)
  except:
    print("Error while deleting file ", filePath)
