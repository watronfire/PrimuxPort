from Kmer import Kmer

class Alignment( object ) :
    name = ""
    raw = ""
    sequence = ""
    seqID = 0
    kmers = list()
    sequenceMap = dict()

    def __init__( self, name, sequence, ID ) :

        # Basic Assignment
        self.name = name
        self.raw = sequence
        self.sequence = self.raw.replace( "-","" )
        self.seqID = ID

        # Functional Assignments
        self.generateSequenceMap()
        self.generateKmers( 18 )

    # Map of each position in a sequence to its position in the alignment
    def generateSequenceMap( self ):
        self.sequenceMap = dict()
        sequenceIndex = 0
        for alignIndex, letter in enumerate( self.raw ):
            if letter != "-":
                self.sequenceMap[sequenceIndex] = alignIndex
                sequenceIndex += 1

    def generateKmers( self, k ):
        self.kmers = list()
        for i in range( 0, len( self.sequence ) - k + 1 ):
            currentMer = str( self.sequence[i:i+k] )
            if currentMer != "N" * k:
                self.kmers.append( Kmer( currentMer, i, self.sequenceMap[i], self.seqID ) )

    def getPosition( self, index, fromRaw ):
        if fromRaw:
            return self.raw[index]
        else:
            return self.sequence[index]

    def filterKmers( self, lowTm, highTm ):
        filteredKmers = [kmer for kmer in self.kmers if lowTm < kmer.primer.Tm < highTm]
        print("{} kmers of {} removed by Tm requirements for sequence {}".format( len( self.kmers ) - len( filteredKmers ), len( self.kmers ), self.seqID ) )
        self.kmers = filteredKmers

# Basic parser for alignment in a fasta format. Not particularly fast as kmers are determined and melting temperatures are calculated.
def parseAlignment( filePath ):
    count = 1
    returnList = list()
    with open( filePath, "r" ) as alignment:
        for line in alignment:
            alignmentName = line.strip().replace( ">", "" )
            alignmentString = alignment.next().strip()
            returnList.append( Alignment( alignmentName, alignmentString, count ) )
            print("{} sequences parsed from alignment.".format( count ))
            count += 1
    return returnList

