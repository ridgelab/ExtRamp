#! /usr/bin/env python3
'''
This program is intended to extract a ramp of slowly translated codons from the beginning
of a gene sequence (DNA or RNA). 

'''
import statistics
import numpy as np
from scipy.stats import gmean
from math import exp, log, sqrt
import gzip
import csv
import argparse
import sys
import re
from multiprocessing import Pool, freeze_support
 
def makeArgParser():
    parser = argparse.ArgumentParser(description='Extract the individual Ramp sequences from a collection of genes')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) gzip fasta file containing gene sequences')
    parser.add_argument('-a', '--tAI', type=str, help='(input) csv file containing the species tAI values')
    parser.add_argument('-o', '--ramp', type=str, help='(output) fasta file to write ramp sequences to')
    parser.add_argument('-v', '--vals', type=str, help='(output) csv file to write tAI/proportion values for each gene to')
    parser.add_argument('-n', '--noRamp', type=str, help='(output) txt file to write the gene names that contained no ramp sequence')
    parser.add_argument('-t', '--threads', type=int, help='the number of threads you want to run, default = 10')
    parser.add_argument('-w', '--window', type=int, default = 9, help='the ribosome window size in codons, default = 9 codons')
    parser.add_argument('-d', '--stdev', type=float, default = 2.0, help='the number of standard deviations below the mean the cutoff value will be')
    parser.add_argument('-m', '--middle', type=str, default = 'gmean', help='the type of statistic used to measure the middle (consensus) efficiency')
    parser.add_argument('-r', '--rna', action="store_true", help='flag for RNA sequences. Default: DNA')
    return parser.parse_args()
     
def readSeqFile(args):
    """
    Input: System arguments
    Returns an array of tuples (Name of sequence, Sequence) created from the input fasta file. 
    Sequences which are not divisible by three or have non-standard nucleotide characters are removed.
    """
    seqArray = []  
    curSeqName = ''
    curSeq = '' 
    if args.input.endswith('.gz'):
        inFile = gzip.open(args.input,'rt')
    else:
        inFile = open(args.input, 'r')
    badSeq = 0
    for line in inFile:
        if line[0] == '>':
            if curSeq != '':
                if args.rna:
                    curSeq.replace('U','T')
                match = re.search('[^ATCG]',curSeq)
                if match == None and len(curSeq)%3 ==0:
                    seqArray.append((curSeqName,curSeq))  
                else:
                    badSeq += 1
            curSeqName = line.strip('\n')
            curSeq = ''
        else:
            curSeq += line.strip('\n').upper()  
    if args.rna:
        curSeq.replace('U','T')
    match = re.search('[^ATCG]',curSeq)
    if match == None and len(curSeq)%3 ==0:
        seqArray.append((curSeqName,curSeq))  
    else:
        badSeq += 1
    if badSeq > 0:
        sys.stderr.write(str(badSeq)+ ' sequences contained non-standard nucleotide characters or were not divisible by three, so were removed!\n')
    inFile.close()
    return seqArray
 
def calcCodonSpeeds(args,seqArray):
    """
    Input: System arguments and an array of tuples (header, sequence) of genes.
    Returns a dictionary of codons to relative translation Speed. If the tAI's for
    the species are provided, it reads them from the file. If tAI values are not
    provided it calculates the codon proportions using all of the gene sequences and
    then uses those values to determine translation speed. 
    """
    sys.stderr.write('Calculating Codon Speeds...\n')
    codonToSpeed = {}
    totalCodons = 0
    if args.tAI != None:
        codonToSpeed = csvToDict(args.tAI)
    else:
        codonGroups = AAtoCodonsDict()
        for elem in seqArray:
            codonToSpeed, totalCodons = countCodons(elem[1], codonToSpeed, totalCodons) 
        for elem in codonGroups:
            totalCounts = 0
            for item in elem:
                totalCounts += codonToSpeed[item]
            for item in elem:
                codonToSpeed[item] = codonToSpeed[item]/totalCounts
    return codonToSpeed

def csvToDict(filename):
    """
    Input: path to a csv file with tAI values where the first row is a list of codons 
        and the second row is the list of tAI values for those codons. The csv file could
        also have a codon and its tAI value on the same line, each codon separated by a new line.
    Returns a dictionary of codons to their tAI speed from the user provided tAI file."""
    tempDict = {}
    try:
        with open(filename, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter= ' ')
            for row in reader:
                codon,tAI = row[0].split(',')
                if tAI == 0:
                     tAI = 0.0001
                tempDict[str(codon)] = float(tAI)
    except:
        with open(filename, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter= ',')
            count = 0
            codonArray = []
            for row in reader:
                if count%2 == 0:
                    codonArray = row
                else:
                    for i in range(len(row)):
                        if float(row[i]) == 0:
                            tempDict[codonArray[i]] = 0.0001
                        else:   
                            tempDict[codonArray[i]] = float(row[i])
                count += 1
    return tempDict

def AAtoCodonsDict():
    '''
    returns a list of a list of codons that encode each amino acid. 
    Groupings are in the following order (one letter code):
    I,L,V,F,M,C,A,G,P,T,S,Y,W,Q,N,H,E,D,K,R,
    '''
    return [["ATT","ATC","ATA"],["CTT","CTC","CTA","CTG","TTA","TTG"],["GTT","GTC","GTA","GTG"],["TTT","TTC"],["ATG"],["TGT","TGC"],["GCT","GCC","GCA","GCG"],["GGT","GGC","GGA","GGG"],["CCT","CCC","CCA","CCG"],["ACT","ACC","ACA","ACG"],["TCT","TCC","TCA","TCG","AGT","AGC"],["TAT","TAC"],["TGG"],["CAA","CAG"],["AAT","AAC"],["CAT","CAC"],["GAA","GAG"],["GAT","GAC"],["AAA","AAG"],["CGT","CGC","CGA","CGG","AGA","AGG"]]
 
def countCodons(sequence, codonToSpeed, totalCodons):
    """
    Input: A DNA sequence
    Returns a dictionary of Codons to the total number 
    of times the codon was found in all provided gene sequences. These are later used
    to find the proportions of codon usage (estimate speed of translation). 
    Also returns the total number of codons.
    """
    for i in range(0,len(sequence),3):
        if sequence[i:i+3] not in codonToSpeed:
            codonToSpeed[sequence[i:i+3]] = 0
        totalCodons += 1
        codonToSpeed[sequence[i:i+3]] += 1
    return codonToSpeed, totalCodons

def calcSeqSpeed(elem):
    """
    Maps the gene name to an array of the tAI/speed values for the given sequence.
    Input: tuple (header, sequence) of a fasta record
    Output: A dictionary of the translational speed of each sequence
    """
    seqToSpeed = {}
    seq = elem[1]     
    speedArray = []
    for index in range(0,len(seq),3):
        stopCodons = ['TAA','TAG','TGA']

        if not seq[index:index+3] in stopCodons: 
            speedArray.append(codonToSpeed[seq[index:index+3]])
        else:
            break
    seqToSpeed[elem[0]] = speedArray
    return seqToSpeed  

def createSpeedsDict(result):
    """
    Returns a dictionary of all the sequence names to the array of their tAI/speed values.
    """
    seqToSpeed = {}
    for elem in result:
        seqToSpeed.update(elem)
    return seqToSpeed

def findGeoSTD(array, gMean):
    '''
    Calculates geometric standard deviation from array and geometric mean.
    '''
    summation = 0
    for item in array:
        summation += (log(item/gMean))**2
    return exp(sqrt(summation/(len(array)-1)))
 
def findConsensusSpeeds(seqToSpeed):
    """
    Returns an array of the Consensus tAI/speed values calculated from all the given
    genes. It is used to determine the cutoff point for the Ramp Sequences.
    """
    consensusSpeed = []
    for i in range(200):
        tAIsAtPosition = []
        for elem in seqToSpeed:
            if i < len(seqToSpeed[elem]):
                tAIsAtPosition.append(seqToSpeed[elem][i])
        if args.middle == 'mean':
            consensusSpeed.append(statistics.mean(tAIsAtPosition))
        elif args.middle == 'gmean':
            consensusSpeed.append(gmean(tAIsAtPosition))      
        elif args.middle == 'median':
            consensusSpeed.append(statistics.median(tAIsAtPosition))
    return consensusSpeed

def calcRiboSpeed(speeds):
    """
    Returns the average of the tAI/speed values in the ribosomeWindow.
    """
    speed = np.longfloat(1.0)
    for elem in speeds:
        speed *= elem
    return float(speed**(1/(len(speeds))))
 
def writeSpeedsFile(seqArray,consensusSpeed):
    """
    Calls findSpeeds function and writes the sequence speeds to a csv file.
    """
    p = Pool(args.threads)
    speedSeqs = p.map(findSpeeds, seqArray) 
    csvfile = open(args.vals, 'w', newline='')
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(["seq", 'position', 'speed_value'])
    for i in range(len(consensusSpeed)):
        writer.writerow(['consensus', i, consensusSpeed[i]])
    from tqdm import tqdm
    for item in tqdm(speedSeqs):
        for row in item:
            writer.writerow(row)
    csvfile.close()

def findSpeeds(elem):
    """
    Calculates the speeds for each sequence and returns them in a csv file format.
    """
    riboSmoothedSpeed = []
    csvLines = []
    for i in range(len(seqToSpeed[elem[0]])-ribosomeWindowLength):
        riboSmoothedSpeed.append(calcRiboSpeed(seqToSpeed[elem[0]][i:i+ribosomeWindowLength]))  
        csvLines.append([elem[0],i,riboSmoothedSpeed[-1]])
    return csvLines

def outputRampSeqs(seqArray):
    """
    Calls isolateRamp function and prints the ramp sequences to the terminal
    or a fasta file
	Input: tuple of header,sequence
    """
    p = Pool(args.threads)  
    rampSeqs = p.map(isolateRamp, seqArray)
    totalSeqs = len(rampSeqs)
    rampSeqs = qualityCheck(rampSeqs)
    outPut = sys.stdout
    if args.ramp:
        outPut = open(args.ramp, 'w')
    if args.noRamp:
        noRampFile = open(args.noRamp, 'w')
    count = 0
    for line in rampSeqs:
        if not line.startswith('None'):
            count += 1
            outPut.write(line)
        elif args.noRamp:
            noRampFile.write(line[4:])
    outPut.close()
    sys.stderr.write("\n" + str(count) + " Ramp Sequences found out of " + str(totalSeqs) + " total sequences\n")

def isolateRamp(elem):
    """
    Calculates the cut off point for the individual ramp sequences and returns them.
    """
    i = 0 
    while calcRiboSpeed(seqToSpeed[elem[0]][i:i+ribosomeWindowLength]) < cutOffVal and i+ribosomeWindowLength <= len(seqToSpeed[elem[0]]):
        i += 1
    if i == 0:
        return 'None' + elem[0] + '\n'
    return elem[0] +  '\n' + elem[1][:(i+ribosomeWindowLength)*3] + '\n'
    
def qualityCheck(rampSeqs):
    """
    Look at the lengths of the ramp sequences and only use those that are within 2 standard
    deviations of the mean length. The distribution is right skewed so the data is log transformed. 
    return the list of rampSeqs with the extremes removed
    """
    lengths = []
    tempSeqs = []
    for i in range(len(rampSeqs)):
        if not rampSeqs[i].startswith('None'):
            ramp = rampSeqs[i].split('\n')
            lengths.append(log(len(ramp[1])))
            tempSeqs.append(i)
    meanLen = statistics.mean(lengths)
    std = statistics.stdev(lengths, meanLen)
    for index in tempSeqs:
        ramp = rampSeqs[index].split('\n')
        if log(len(ramp[1])) < meanLen - (2 * std) or log(len(ramp[1])) > meanLen + (2 * std):
            newLine = 'None' + ramp[0] + '\n'
            rampSeqs[index] = newLine
    return rampSeqs

if __name__ == '__main__':
    freeze_support
    args = makeArgParser()
    sys.stderr.write('Reading Sequences...\n')
    seqArray = readSeqFile(args)
    p = Pool(args.threads)
    codonToSpeed = calcCodonSpeeds(args,seqArray)
    sys.stderr.write('Calculating Sequence Speeds...\n')
    p = Pool(args.threads)
    seqToSpeed = createSpeedsDict(p.map(calcSeqSpeed,seqArray))
    
    #create a consensus speed to determine average speed as a cutoff value for the ramp sequence and smooth consensus with ribosomeWindowLength
    ribosomeWindowLength = args.window 
    consensusSpeed = findConsensusSpeeds(seqToSpeed)
    middle = gmean(consensusSpeed[1:])
    STDev = findGeoSTD(consensusSpeed[1:], middle) 
    cutOffVal = middle / (STDev**args.stdev)   
    sys.stderr.write("\nConsensus Codon Efficiency using " + args.middle + ": " + str(middle) + "\n")
    sys.stderr.write("Standard Deviation: " + str(STDev) + "\n")
    sys.stderr.write("Maximum Efficiency in Ramp Sequence " + str(cutOffVal) + "\n\n")

    #output speed values in a csv file if indicated
    if args.vals:
        sys.stderr.write('Writing Speeds File...\n')
        writeSpeedsFile(seqArray,consensusSpeed)

    #write Ramp Sequence to a fasta file
    sys.stderr.write('Isolating Ramp Sequences...\n')      
    outputRampSeqs(seqArray)
             
 
