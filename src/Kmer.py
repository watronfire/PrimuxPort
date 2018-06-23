from primer import Primer


class Kmer( object ):
    
    sequence = ""
    alignLocation = 0
    seqLocation = []
    parent = []
    primer = ""

    def __init__( self, sequence, seqLocation, alignLocation, parent, primer=None ) :
        self.parent = [ parent ]
        self.seqLocation = [ seqLocation ]
        self.alignLocation = alignLocation
        self.sequence = sequence

        if primer is None:
            self.primer = Primer( self.sequence )
        else:
            self.primer = primer

    def addConservation( self, seqLocation, parent ):
        self.seqLocation.extend( seqLocation )
        self.parent.extend( parent )

    def __repr__( self ):
        temperature = "{0:.2f}".format( self.primer.Tm )
        return "{}, in {} sequences and at position {} of alignment. Tm: {}".format( self.sequence, len( self.seqLocation ), self.alignLocation, temperature )