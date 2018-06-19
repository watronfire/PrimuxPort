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
        self.generateKmers( 18 )
        self.generateSequenceMap()

    def generateSequenceMap( self ):
        self.sequenceMap = dict()
        sequenceIndex = 0
        for alignIndex, letter in enumerate( self.raw ):
            if letter != "-":
                sequenceIndex += 1
                self.sequenceMap[sequenceIndex] = alignIndex

    def generateKmers( self, k ):
        for i in range( 0, len( self.sequence ) - k + 1 ):
            currentMer = str( self.sequence[i:i+k] )
            self.kmers.append( Kmer( currentMer, i, i, self.seqID ) )

    def getPosition( self, index, fromRaw ):
        if fromRaw:
            return self.raw[index]
        else:
            return self.sequence[index]

def parseAlignment( filePath ):
    count = 1
    returnList = list()
    with open( filePath, "r" ) as alignment:

        alignmentName = alignment.next().strip().replace( ">", "" )
        alignmentString = alignment.next().strip()
        for line in alignment:
            if ">" in line:
                returnList.append( Alignment( alignmentName, alignmentString, count ) )
                count += 1
                alignmentName = line.strip().replace( ">", "" )
                alignmentString = ""
            else:
                alignmentString += line.strip()
        returnList.append( Alignment( alignmentName, alignmentString, count ) )
    return returnList



