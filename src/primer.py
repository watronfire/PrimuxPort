from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

class Primer( object ):
    sequence = ""
    Tm = 0.0

    def __init__( self, sequence ):
        self.sequence = sequence
        self.Tm = mt.Tm_NN( Seq( self.sequence ), nn_table=mt.DNA_NN4, Na=50, Mg=2.0, dNTPs=0.2  )

    def __repr__(self):
        return "{}: {:.2f}C".format( self.sequence, self.Tm )

    #def fastaRepr( self ):
    #    return "> primer-{}\n{}".format( self.ID, self.sequence )