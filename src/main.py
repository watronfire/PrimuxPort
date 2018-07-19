import alignmentIO

IUPACbases = [ 'A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N' ]
basePairs = [ a + b for a in IUPACbases for b in IUPACbases ]
codes = 'AMRWMRWVHDVHDNNMCSYMVHSYBVHNBNRSGKVRDSBKVNDBNWYKTHDWBYKNHDBNMMVHMVHVHNVHNNNRVRDVRDVNDVNDNNWHDWHDWNHDNHDNNVSSBVVNSBBVNNBNHYBYHNHBYBNHNBNDBKKNDDBBKNNDBNVVVNVVNVNNVNNNNHHNHHNHNHNNHNNNDNDDNDDNNDNNDNNNBBBNNNBBBNNNBNNNNNNNNNNNNNNNN'
AmbiqDict = dict( zip( basePairs, codes ) )

IUPACExpandDict = { "A" : ["A"],
                    "T" : ["T"],
                    "G" : ["G"],
                    "C" : ["C"],
                    "M" : ["A", "C"],
                    "R" : ["A", "G"],
                    "W" : ["A", "T"],
                    "S" : ["C", "G"],
                    "Y" : ["T", "C"],
                    "K" : ["T", "G"],
                    "V" : ["A", "C", "G"],
                    "H" : ["A", "T", "C"],
                    "D" : ["A", "T", "G"],
                    "B" : ["T", "C", "G"],
                    "N" : ["A", "T", "C", "G"] }

# Collapses a dictionary of kmers mapped to alignment locations, where collaps
def collapseKmers( alignment ):
    kmerDict = dict()
    for item in alignment:
        for kmer in item.kmers:
            if kmer.alignLocation in kmerDict:
                if kmer.sequence in kmerDict[kmer.alignLocation]:
                    kmerDict[kmer.alignLocation][kmer.sequence].addConservation( kmer.seqLocation, kmer.parent )
                else:
                    kmerDict[kmer.alignLocation][kmer.sequence] = kmer
            else:
                kmerDict[kmer.alignLocation] = { kmer.sequence : kmer }

    return kmerDict

def kmerOutput( kmers ):
    kmerList = list( kmers.keys() )
    kmerList.sort()
    for item in kmerList :
        if len( kmers[item] ) == 1:
            for variant in kmers[item].keys():
                print( "{}, {}, {}".format( item, len( kmers[item] ), len( kmers[item][variant].parent  ) ) )
        else:
            print( "{}, {}".format( item, len( kmers[item] ) ) )

def generateKmerTree( sequence ):

    kmerTree = list()

    for i in range( len( sequence ) ):
        kmerTree.append( [sequence[k:i+k] for k in range( len( sequence ) - i  + 1) ] )
    return kmerTree

def calculateDegeneracy( seqA, seqB ):
    returnDegen = 1
    for i in range( len( seqA ) ):
        baseA = seqA[i]
        baseB = seqB[i]
        if baseA != baseB:
            returnDegen *= len( set( IUPACExpandDict[baseA] ).union( set( IUPACExpandDict[baseB] ) ) )
        elif baseA not in "ATCG":
            returnDegen *= len( IUPACExpandDict[baseA] )
    return returnDegen

def mergeSequece( seqA, seqB ):
    returnSeq = ""
    for a, b in zip( seqA, seqB ):
        returnSeq += AmbiqDict[a+b]
    return returnSeq

# Finds a pair of sequences such that occurance bonus is maximized and degeneracy cost is minimized.
def optimizeCollapse( sequenceList ):

    maxPair = 0,0
    maxTC = 0

    for i in range( len( sequenceList ) ):
        for k in range( i, len( sequenceList ) ):
            if i != k:

                # Occurance bonus is the amount of sequences the kmer pair would bind to divided by the total number of
                # sequences.
                occuranceBonus = ( sequenceList[i][1] + sequenceList[k][1] ) / 66.0

                # Degeneracy is the number of sequences which a degenerate sequence describes. A measure of entropy, essentially
                degeneracy = calculateDegeneracy( sequenceList[i][0], sequenceList[k][0] )

                # I don't like the variable name total cost but I don't have a better description. TODO: replace totalCost.
                totalCost = occuranceBonus / degeneracy

                # Find the maximum total cost of all combinations of kmers.
                if totalCost > maxTC:
                    maxTC = totalCost
                    maxPair = i,k

    return maxPair, maxTC

with open( "/Users/natem/Documents/Code/Python/PriMuxPort/src/temp.txt", "r" ) as seqInput:
    seqList = list()

    for line in seqInput:
        lineSplit = line.strip().split(",")
        seqList.append( [lineSplit[0], int(lineSplit[1])] )

pairing, tc = optimizeCollapse( seqList )
outputSeq = mergeSequece( seqList[pairing[0]][0], seqList[pairing[1]][0])
print( seqList[pairing[0]][0] )
print( seqList[pairing[1]][0] )
print( outputSeq )

#file = "/Users/natem/Documents/Sequences/LassaL_Sierra_Leone_Alignment.fasta"
#lassaAlignment = alignmentIO.parseAlignment( file )
#for thing in lassaAlignment:
#    thing.filterKmers( 63, 67 )
#totalKmers = collapseKmers( lassaAlignment )
#print( totalKmers[6031] )
#kmerOutput( totalKmers  )

