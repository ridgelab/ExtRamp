#! /usr/bin/env python3
'''
This program is intended to extract a ramp of slowly translated codons from the beginning
of a gene sequence (DNA or RNA). 
'''

import statistics
from statistics import mean,median
import numpy as np
from scipy.stats import gmean, hmean
from math import exp, log, sqrt,ceil
import gzip
import csv
import argparse
import sys
import re
from multiprocessing import Pool, freeze_support

def makeArgParser():
    parser = argparse.ArgumentParser(description='Extract the individual Ramp sequences from a collection of genes')
    parser.add_argument('-i', '--input', type=str, required=True, help='(input) Required fasta file containing gene sequences (.gz or .gzip for gzipped)')
    parser.add_argument('-a', '--tAI', type=str, help='(input) csv file containing the species tAI values')
    parser.add_argument('-u', '--rscu', type=str, help='(input) fasta file used to compute relative synonymous codon usage and relative adaptiveness of codons')
    parser.add_argument('-o', '--ramp', type=str, help='(output) fasta file to write ramp sequences to')
    parser.add_argument('-v', '--verbose', action="store_true", help='flag to print progress to standard error')
    parser.add_argument('-l', '--vals', type=str, help='(output) csv file to write tAI/relative adaptiveness values for ribosome-smoothed efficiency at each window')
    parser.add_argument('-p', '--speeds', type=str, help='(output) speeds file to write tAI/relative adaptiveness values for each position in the sequence from the codon after the start codon to the codon before the stop codon. Format: Header newline list of values')
    parser.add_argument('-n', '--noRamp', type=str, help='(output) txt file to write the gene names that contained no ramp sequence')
    parser.add_argument('-z', '--removedSequences', default = None, type=str, help='(output) Write the header lines that are removed (e.g., sequence not long enough or not divisible by 3) to output file')
    parser.add_argument('-x', '--afterRamp', type=str, required=False, help='(output) fasta file containing gene sequences after the identified ramp sequence')
    parser.add_argument('-t', '--threads', type=int, help='the number of threads you want to run, default = all')
    parser.add_argument('-w', '--window', type=int, default = 9, help='the ribosome window size in codons, default = 9 codons')
    parser.add_argument('-s', '--stdev', type=float, default = -1.0, help='the number of standard deviations below the mean the cutoff value to be included as a ramp for gmean and mean. Not used by default')
    parser.add_argument('-d', '--stdevRampLength', type=float, default = -1.0, help='the number of standard deviations used in the quality check step (lengths of the ramp sequences). Not done by default')
    parser.add_argument('-m', '--middle', type=str, default = 'hmean', help='the type of statistic used to measure the middle (consensus) efficiency. Options: hmean, gmean, mean, median')
    parser.add_argument('-r', '--rna', action="store_true", help='flag for RNA sequences. Default: DNA')
    parser.add_argument('-f', '--determine_cutoff', action="store_true", help='flag to determine outlier percentages for mean cutoff based on species FASTA file. Default: local minimum in first 8 percent of gene')
    parser.add_argument('-c', '--cutoff', type=int, default = 8, help='Cutoff for where the local minimum must occur in the gene for a ramp to be calculated. If --determine_cutoff (-f) is used, then this value may change. Is not used if standard deviations are set.')
    parser.add_argument('-e', '--determine_cutoff_percent', type=str, default = "outlier", help='Cutoff for determining percent of gene that is in an outlier region. Used in conjunction with -f. Default is true outliers. Other options include numbers from 0-99, which indicate the region of a box plot. For instance, 75 means the 75th quartile or above.')
    parser.add_argument('-q', '--seqLength', type=float, default = 300, help='Minimum nucleotide sequence length. Default is 100 amino acids * 3 = 300 nucleotides')
    args = parser.parse_args()
    if args.middle == 'hmean' and args.stdev >=0:
        sys.stderr.write("Warning: hmean cannot be used with standard deviations. -s will be ignored.\n")
        args.stdev = -1.0
    if (args.window *3) >args.seqLength:
        args.seqLength = args.window*3 +1
    return args
     
def readSeqFile(args, inputFile):
    """
    Input: System arguments
    Returns an array of tuples (Name of sequence, Sequence) created from the input fasta file. 
    Sequences which are not divisible by three or have non-standard nucleotide characters are removed.
    """
    seqArray = []  
    curSeqName = ''
    curSeq = '' 
    inFile = ""
    stop = {'TAA','TGA','TAG'}
    if inputFile.endswith('.gz') or inputFile.endswith('.gzip'):
        inFile = gzip.open(inputFile,'rt')
    else:
        inFile = open(inputFile, 'r')
    badSeq = 0
    removedSeq = ""
    if args.removedSequences and inputFile==args.input:
        removedSeq = open(args.removedSequences,'w')

    for line in inFile:
        if line[0] == '>':
            if curSeq != '':
                if args.rna:
                    curSeq.replace('U','T')
                match = re.search('[^ATCG]',curSeq)
                if match == None and len(curSeq)%3 ==0 and len(curSeq) >=args.seqLength:
                    startCodon = curSeq[0:3] 
                    curSeq = curSeq[3:] #remove start codon
                    stopCodon = curSeq[-3:] 
                    curSeq = curSeq[0:-3] #remove stop codon
                    seqArray.append((curSeqName,curSeq,startCodon,stopCodon))  
                else:
                    if args.removedSequences and inputFile == args.input:
                        removedSeq.write(curSeqName +"\n")
                    badSeq += 1
            curSeqName = line.strip('\n')
            curSeq = ''
        else:
            curSeq += line.strip('\n').upper()  
    if args.rna:
        curSeq.replace('U','T')
    match = re.search('[^ATCG]',curSeq)
    if match == None and len(curSeq)%3 ==0 and len(curSeq) >=args.seqLength:
        startCodon = curSeq[0:3] 
        curSeq = curSeq[3:] #remove start codon
        stopCodon = curSeq[-3:] 
        curSeq = curSeq[0:-3] #remove stop codon
        seqArray.append((curSeqName,curSeq,startCodon,stopCodon))  
    else:
        badSeq += 1
    if badSeq > 0 and args.verbose and inputFile==args.input:
        sys.stderr.write( '\n' + str(badSeq)+ ' sequences contained non-standard nucleotide characters or were not divisible by three, so were removed!\n\n')
    inFile.close()
    if args.removedSequences and inputFile==args.input:
        removedSeq.close()
    return seqArray
 
def calcCodonSpeeds(seqArray):
    """
    Input: an array of tuples (header, sequence) of genes.
    Returns a dictionary of codons to relative translation Speed. 
    It calculates the codon proportions using all of the gene sequences and
    then uses those values to determine translation speed. Uses metrics from CAI values
    Relative Synonymous Codon Usage (RSCU) and
    Relative adaptiveness of codon (w) are calculated based on equations presented in
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC340524/pdf/nar00247-0410.pdf
    """
    codonToSpeed = {}
    totalCodons = 0
    codonGroups = AAtoCodonsDict()
    for elem in seqArray:
        codonToSpeed, totalCodons = countCodons(elem[1],codonToSpeed,totalCodons) 
    for aa in codonGroups:
        totalCounts = 0
        numCodons = 0
        for codon in aa: 
            numCodons +=1
            if not codon in codonToSpeed:
                codonToSpeed[codon] = 0.0001 #codonToSpeed has count of codon usages
                continue
            totalCounts += codonToSpeed[codon]
        maxRSCU = 0.0
        for codon in aa:
            if numCodons ==0 or totalCounts==0:
                rscu = 0.0001
                codonToSpeed[codon] = rscu  #codonToSpeed has RSCU
                continue
            rscu = codonToSpeed[codon] / ((1.0/numCodons)*totalCounts)   
            if rscu > maxRSCU:
                maxRSCU = rscu
            codonToSpeed[codon] = rscu  #codonToSpeed has RSCU
        for codon in aa:
            if maxRSCU == 0:
                w = 0.0001
                codonToSpeed[codon] = w #codonToSpeed has relative adaptiveness of codons
                continue
            w = codonToSpeed[codon] / maxRSCU #w=relative adaptiveness of codons
            codonToSpeed[codon] = w #codonToSpeed has relative adaptiveness of codons
    return codonToSpeed

def csvToDict(filename):
    """
    Input: path to a csv file with tAI values where the first row is a list of codons 
        and the second row is the list of tAI values for those codons. The csv file could
        also have a codon and its tAI value on the same line, each codon separated by a new line.
    Returns a dictionary of codons to their tAI speed from the user provided tAI file.
    """
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
 
def countCodons(sequence,codonToSpeed,totalCodons):
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
    Output: A dictionary of seq name to a tuple (translational speed of each sequence, cutOffValue for rampSequence)
    """

    seqToSpeed = {}
    seq = elem[1]
    speedArray = []
    for index in range(0,len(seq),3):
        stopCodons = ['TAA','TAG','TGA']
        if not seq[index:index+3] in stopCodons: 
            speedArray.append( codonToSpeed[seq[index:index+3]])
        else:
            speedArray.append(0.0001)
    
    seqToSpeed[elem[0]] = tuple(speedArray)
    return seqToSpeed  

def createSpeedsDict(result):
    """
    Returns a dictionary of all the sequence names to the array of their tAI/speed values.
    """
    seqToSpeed = {}
    for elem in result:
        seqToSpeed.update(elem)
    return seqToSpeed

def findStDev(array,mean):
    '''
    Calculates geometric standard deviation from array and geometric mean.
    '''
    if args.middle != 'gmean': #If the user specified other mean
        return statistics.stdev(array,mean)
    ##for gMean
    summation = 0
    for item in array:
        summation += (log(item/mean))**2
    return exp(sqrt(summation/len(array)))

def calcRiboSpeed(speeds):
    """
    Returns the average of the tAI/speed values in the ribosomeWindow.
    """
    if len(speeds)==0:
        return 0.0001
    elif args.middle == 'hmean':
        return hmean(speeds)
    elif args.middle == 'gmean':
        return gmean(speeds)
    elif args.middle == 'mean':
        return statistics.mean(speeds)
    elif args.middle == 'median':
        return statistics.median(speeds)
 
def writeSpeedsFile(speedSeqs):
    """
    Writes the sequence speeds to a csv file.
    """
    csvfile = open(args.vals, 'w', newline='')
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(["seq", 'position', 'speed_value'])
    from tqdm import tqdm
    for item in tqdm(speedSeqs):
        for row in item:
            writer.writerow(row)
    csvfile.close()

def findSpeeds(elem):
    """
    Calculates the speeds for each sequence and returns them in a csv file format.
    """
    csvLines = []
    for i in range(len(seqToSpeed[elem[0]])-ribosomeWindowLength):
        riboSmoothedSpeed = middleFunc[args.middle](seqToSpeed[elem[0]][i:i+ribosomeWindowLength])
        csvLines.append([elem[0],i+1 ,riboSmoothedSpeed]) #i +1 because start codon was taken away
    return csvLines

def outputRampSeqs(rampSeqs,args):
    """
    Calls isolateRamp function and prints the ramp sequences to the terminal
    or a fasta file
    Input: tuple of header,sequence
    """
    totalSeqs = len(rampSeqs)
    noRampFile = ""
    if args.stdevRampLength >=0:
        rampSeqs = qualityCheck(rampSeqs,args.stdevRampLength)
    outPut = sys.stdout
    if args.ramp:
        outPut = open(args.ramp, 'w')
    if args.afterRamp:
        afterRampFile = open(args.afterRamp,'w')
    if args.noRamp and args.determine_cutoff:
        noRampFile = open(args.noRamp, 'a')
    elif args.noRamp:
        noRampFile = open(args.noRamp, 'w')
    count = 0
    for line in rampSeqs:
        if not args.afterRamp:
            if not line.startswith('None'):
                count += 1
                outPut.write(line)
            elif args.noRamp:
                noRampFile.write(line[4:])
        else:
            if not line[0].startswith('None'):
                count += 1
                outPut.write(line[0])
                afterRampFile.write(line[1])
            elif args.noRamp:
                noRampFile.write(line[0][4:])
    outPut.close()
    if args.verbose:
        sys.stderr.write("\n" + str(count) + " Ramp Sequences found\n")
    if args.noRamp:
        noRampFile.close()

def isolateRamp(record):
    """
    Calculates the cut off point for the individual ramp sequences and returns them.
    """
    mean = middleFunc[args.middle](seqToSpeed[record[0]])
    stdev = findStDev(seqToSpeed[record[0]],mean)
    cutOffVal = mean - (args.stdev*stdev)
    i = 0 
    while calcRiboSpeed(seqToSpeed[record[0]][i:i+ribosomeWindowLength]) < cutOffVal and (i+ribosomeWindowLength) <= len(seqToSpeed[record[0]]):
        i += 1
    if i == 0:
        if not args.afterRamp:
            return 'None' + record[0] + '\n'
        else:
            return tuple(['None' + record[0] + '\n'])

    if not args.afterRamp:
        return record[0] +  '\n' + record[2] + record[1][:(i+ribosomeWindowLength)*3] + '\n'
    else:
        return tuple([record[0] +  '\n' + record[2] + record[1][:(i+ribosomeWindowLength)*3] + '\n',record[0] +  '\n' + record[1][(i+ribosomeWindowLength)*3:] +record[3]+ '\n'])

def isolateRampHmean(record):
    speeds = seqToSpeed[record[0]]
    mean = hmean(speeds)
    windowMeans = []
    for z in range(len(speeds)-ribosomeWindowLength):
        windowMeans.append(middleFunc[args.middle](speeds[z:z+ribosomeWindowLength]))
    bestMean = min(windowMeans)
    pos = [index for index, value in enumerate(windowMeans) if value == bestMean] #ensures that all local minima are included
    perc = float(ceil(100*(pos[0]/float(len(windowMeans))))) #Percentages x 100
    if perc <=percentThatIsRamp:
        i =pos[0]
        for x in range(i,len(windowMeans)):
            if windowMeans[x] >= mean:
                if not args.afterRamp:
                    return record[0] +  '\n' + record[2] + record[1][:(x+ribosomeWindowLength)*3] + '\n'
                else:
                    return tuple([record[0] +  '\n' + record[2] + record[1][:(x+ribosomeWindowLength)*3] + '\n',record[0] +  '\n' + record[1][(x+ribosomeWindowLength)*3:] +record[3]+ '\n'])
    if not args.afterRamp:
        return 'None' + record[0] + '\n'
    else:
        return tuple(['None' + record[0] + '\n'])

def qualityCheck(rampSeqs,numStDev):
    
    """
    Look at the lengths of the ramp sequences and only use those that are within numStDev standard
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
    if len(lengths) == 0:
        if args.verbose:
            sys.stderr.write('NO RAMP SEQUENCES FOUND\n')
        sys.exit()
    if len(lengths) == 1:
        return rampSeqs

    meanLen = statistics.mean(lengths)
    std = statistics.stdev(lengths, meanLen)
    removed = 0
    for index in tempSeqs:
        ramp = rampSeqs[index].split('\n')
        if log(len(ramp[1])) < meanLen - (numStDev * std) or log(len(ramp[1])) > meanLen + (numStDev * std):
            newLine = 'None' + ramp[0] + '\n'
            rampSeqs[index] = newLine
            removed +=1
    if args.verbose:
        sys.stderr.write(str(removed) + " sequences removed by standard deviation check of all ramp sequences (-d option).\n")
    return rampSeqs


def getBottleneck(header):
    allPercents = []
    speeds = seqToSpeed[header]
    windowMeans = []
    for i in range(len(speeds)-args.window):
        windowMeans.append(middleFunc[args.middle](speeds[i:i+args.window]))
    bestMean = min(windowMeans)
    pos = [index for index, value in enumerate(windowMeans) if value == bestMean] #ensures that all local minima are included
    for p in pos:
        perc = float(ceil(100.0*((p+1)/float(len(windowMeans))))) #Percentages x 100
        allPercents.append((perc,header,p))
    return allPercents

def getCutoffValue(counts):
    '''
    Input: a dictionary of header lines -> tuple array of codon efficiency values for each codon
    Output: Locations that have more local minimum codon efficiencies and are outliers, starting from the front of the gene. If no outliers exist, return 0.
    '''
    instances = list(counts.values())
    instances = np.array(instances)
    upper = 0
    if args.determine_cutoff_percent.isdigit():
        upper = np.percentile(instances,int(args.determine_cutoff_percent))
    else:
        q1= np.percentile(instances,25)
        q3= np.percentile(instances,75)
        upper = q3 + (1.5*(q3-q1))
    outliers = set()
    for percent in instances:
        if percent >=upper:
            outliers.add(percent)
    if args.verbose:
        sys.stderr.write("The following percentages exceed the threshold (" +args.determine_cutoff_percent +"):")
        for x in range(1,101):
            if counts[x] in outliers:
                sys.stderr.write(" " +str(x))
        sys.stderr.write("\n")

    for x in range(1,101):
        #if x in outliers:
        if counts[x] in outliers:
            continue
        return (x-1) #-1 because it goes 1 past the last outlier. Divide by 100 to turn percent into decimal
    return 0

if __name__ == '__main__':
    freeze_support()
    args = makeArgParser()
    noRampFile = ""
    ribosomeWindowLength = args.window 
    middleFunc={'hmean':hmean,'gmean':gmean,'mean':mean,'median':median}
    if not args.middle in middleFunc:
        sys.stderr.write("args.middle must be one of the following options:\nhmean\ngmean\nmean\nmedian\n")
        sys.exit()
    if args.verbose:
        sys.stderr.write('Reading Sequences...\n')
    seqArray = readSeqFile(args, args.input)
    if args.verbose:
        sys.stderr.write("Total Sequences: " +str(len(seqArray)) + '\n')
    if len(seqArray) ==0:
        sys.stderr.write("No sequences passed the initial filter. Ramp sequences were unable to be calculated.\n")
        sys.exit()

    codonToSpeed = {}
    if args.tAI:
        codonToSpeed = csvToDict(args.tAI)
    else:
        if args.verbose:
            sys.stderr.write('Calculating Codon Speeds...\n')
        if args.rscu:
            seqArray_rscu = readSeqFile(args,args.rscu)
            codonToSpeed = calcCodonSpeeds(seqArray_rscu)
        else:
            codonToSpeed = calcCodonSpeeds(seqArray)
    if args.verbose:
        sys.stderr.write('Calculating Sequence Speeds...\n')
    p = Pool(args.threads)
    seqToSpeed = createSpeedsDict(p.map(calcSeqSpeed,seqArray)) ##header -> tuple of array of codon efficiency values
    if args.speeds:
        if args.verbose:
            sys.stderr.write("Writing speeds data to file.\n")
        output_speeds = open(args.speeds,'w')
        for header in seqToSpeed:
            output_speeds.write(header + "\n" + str(seqToSpeed[header])[1:-1] + "\n")
        output_speeds.close()
    #create a consensus speed to determine average speed as a cutoff value for the ramp sequence and smooth consensus with ribosomeWindowLength
    percentThatIsRamp =args.cutoff
    if args.determine_cutoff:
        if args.verbose:
            sys.stderr.write('Determining Cutoff Percentage from Ramp Sequences...\n')
        p = Pool(args.threads)
        bottlenecks = p.map(getBottleneck,list(seqToSpeed.keys()))
        counts = {}
        for percents in bottlenecks:
            for p in percents:
                place = int(p[0])
                if not place in counts:
                    counts[place] = 0
                counts[place] +=1
        for perc in range(1,101):
            if not perc in counts:
                counts[perc] = 0
        percentThatIsRamp = getCutoffValue(counts)
        if args.verbose:
            sys.stderr.write('\tThe cutoff percentage is ' +str(percentThatIsRamp) + '%\n')
            sys.stderr.write('Isolating Ramp Sequences...\n')      
        rampSeqs = []
        posOfRamp = {}
        for percents in bottlenecks:
            for p in percents:
                place = int(p[0])
                if place <= percentThatIsRamp:
                    header = p[1]
                    pos = p[2]
                    speeds = seqToSpeed[header]
                    seqMean = middleFunc[args.middle](speeds)
                    for i in range(pos,len(speeds)-args.window):
                        windowMean =middleFunc[args.middle](speeds[i:i+args.window])
                        if windowMean >= seqMean:
                            posOfRamp[header] = i-1
                            break
        if args.noRamp:
            noRampFile = open(args.noRamp,'w')
        for record in seqArray:
            if record[0] in posOfRamp:
                if not args.afterRamp:
                    rampSeqs.append(record[0] +  '\n' + record[2] + record[1][:(posOfRamp[record[0]]+ribosomeWindowLength)*3] + '\n')
                else:
                    rampSeqs.append(tuple([record[0] +  '\n' + record[2] + record[1][:(posOfRamp[record[0]]+ribosomeWindowLength)*3] + '\n',record[0] +  '\n' + record[1][(posOfRamp[record[0]]+ribosomeWindowLength)*3:] +record[3]+ '\n']))
            elif args.noRamp:
                noRampFile.write(record[0] + '\n')
        if args.noRamp:
            noRampFile.close()
        outputRampSeqs(rampSeqs,args)
    #output speed values in a csv file if indicated
    p = Pool(args.threads)
    if args.vals:
        if args.verbose:
            sys.stderr.write('Writing Speeds File...\n')
        speedSeqs = p.map(findSpeeds, seqArray) 
        writeSpeedsFile(speedSeqs)
    #write Ramp Sequence to a fasta file
    if not args.determine_cutoff:
        if args.verbose:
            sys.stderr.write('Isolating Ramp Sequences...\n')      
        if args.stdev >=0:
            rampSeqs = p.map(isolateRamp, seqArray)
            outputRampSeqs(rampSeqs,args)
        else:
            rampSeqs = p.map(isolateRampHmean, seqArray)
            outputRampSeqs(rampSeqs,args)

