
import primer3

class Primer( object ):
    sequence = ""
    Tm = 0.0

    def __init__( self, sequence ):
        self.sequence = sequence
        try:
            # Calculates melting temperature of the primer. Primer3 is quite fast but BioPython provides more constumization.
            #self.Tm = mt.Tm_NN( Seq( self.sequence ), nn_table=mt.DNA_NN4, Na=50, Mg=2.0, dNTPs=0.2  )
            self.Tm = primer3.calcTm( self.sequence, mv_conc=50, dv_conc=2.0, dntp_conc=0.2 )
        except IndexError:
            print( self.sequence )
            exit( 69 )

    def __repr__(self):
        return "{}: {:.2f}C".format( self.sequence, self.Tm )

    #def fastaRepr( self ):
    #     #    return "> primer-{}\n{}".format( self.ID, self.sequence )