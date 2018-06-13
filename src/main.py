from operator import attrgetter
from Bio import SeqIO
from primer import Primer
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# Biopython has a real complicated way of parsing fastas, so this collects just the sequences in a list. Maybe could
# extend to return a dict maintaining sequence IDs.
def generateSequenceList( fastaFile ):
    with open( fastaFile, "r" ) as sequencesFile:
        sequences = SeqIO.parse( sequencesFile, "fasta" )
        sequenceList = [sequence.seq for sequence in sequences]
    return sequenceList

# Iterates through a list of sequences and calculates all kmers. Should extend to collect a range of ks.
def calculateKmers( sequences, k ):
    print( "Calculating {}-mers".format( k ) )
    kmerDict = dict()

    for sequence in sequences:
        for i in range( 0, len( sequence ) - k + 1 ):
            currentMer = str( sequence[i:i+k] )

            if currentMer in kmerDict:
                kmerDict[currentMer] += 1
            else:
                kmerDict[currentMer] = 1

    # Remove entries which are just ambiguous bp.
    del( kmerDict["N"*k] )

    print( "{} {}-mers found".format( len( kmerDict ), k ) )
    return kmerDict

# Converts Kmers to a Primer object which is a nice data structure.
def covertKmerToPrimers( kmers ):
    print( "Converting kmers to primers and calculating characteristics")
    returnList = list()

    primerIDCounter = 1

    for kmer in kmers.keys():
        returnList.append( Primer( kmer, kmers[kmer], primerIDCounter ) )
        primerIDCounter += 1

    return returnList

# Primers with Tms outside the range specified are removed.
def filterPrimersTm( primers, lowTm, highTm ):
    returnList = [ primer for primer in primers if lowTm < primer.Tm < highTm ]
    print( "{} primers removed by Tm requirements".format( len( primers ) - len( returnList ) ) )
    return returnList

def alignPrimers( primers, sequence, mmAllowed ):
    for i in range( len( sequence ) ):
        print( i )

inputFile = "/Users/natem/Documents/Sequences/LassaV_L_Sierra_Leone.fasta"
seqList = generateSequenceList( inputFile )
outputDict = calculateKmers( seqList, 19 )

# Convert collected Kmers to Primers and sort based on occurance.
primerList = covertKmerToPrimers( outputDict )
primerList.sort( key=attrgetter( "occurrence" ) )

primerList = filterPrimersTm( primerList, 63, 67 )

alignPrimers( primerList, seqList[1], 1 )

#testAlignment = pairwise2.align.localms( primerList[1].sequence, seqList[23], 2, -1, -.5, -.1 )