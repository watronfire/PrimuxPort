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

def optimizeCollapse( file ):
    seqList = list()
    with open( file, "r" ) as inputSeqs:
        for line in inputSeqs:
            line = line.strip()
            seqList.append( [line.split(",")[0], int(line.split(",")[1])] )
    print("\t" + "\t\t".join( str(i) for i in range( len( seqList ) ) ) )

    # Prints a occurrance bonus matrix. Essentially,
    for i in range( len(seqList) ):
        outputString = "{}\t".format( i )
        for k in range(len( seqList ) ):
            if i != k:
                occuranceBonus = ( seqList[i][1] + seqList[k][1] ) / 66.0
                outputString += "{:.2f}\t".format( occuranceBonus )
            else:
                outputString += "-\t\t"
        print( outputString )
    print("\n")

    # Prints the Degeneracy Matrix
    print("\t" + "\t".join( str( i ) for i in range( len( seqList ) ) ))
    for i in range( len( seqList ) ):
        outputString = "{}\t".format( i )
        for k in range( len( seqList ) ):
            if i != k:
                degeneracy = 0
                for j in range( len( seqList[i][0] ) ):
                    base1 = seqList[i][0][j]
                    base2 = seqList[k][0][j]
                    if base1 != base2:
                        degeneracy += len( set( IUPACExpandDict[base1] ).union( IUPACExpandDict[base2] ) )
                outputString += "{}\t".format( degeneracy )
            else:
                outputString += "-\t"
        print( outputString )


optimizeCollapse( "/Users/natem/Documents/Code/Python/PriMuxPort/src/temp.txt" )

#file = "/Users/natem/Documents/Sequences/LassaL_Sierra_Leone_Alignment.fasta"
#lassaAlignment = alignmentIO.parseAlignment( file )

#for thing in lassaAlignment:
#    thing.filterKmers( 63, 67 )

#totalKmers = collapseKmers( lassaAlignment )

#print( totalKmers[120] )

# Not quite a Kmer tree but generates all kmers and returns a lame dictionary-like list where the key is equal to k.
# Ranges of kmers can then be extracted using typical means. kmerlist[kn:kn+] where kn+ in exclusive.

