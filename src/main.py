
import regex
import alignmentIO

IUPACbases = [ 'A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N' ]
basePairs = [ a + b for a in IUPACbases for b in IUPACbases ]
codes = 'AMRWMRWVHDVHDNNMCSYMVHSYBVHNBNRSGKVRDSBKVNDBNWYKTHDWBYKNHDBNMMVHMVHVHNVHNNNRVRDVRDVNDVNDNNWHDWHDWNHDNHDNNVSSBVVNSBBVNNBNHYBYHNHBYBNHNBNDBKKNDDBBKNNDBNVVVNVVNVNNVNNNNHHNHHNHNHNNHNNNDNDDNDDNNDNNDNNNBBBNNNBBBNNNBNNNNNNNNNNNNNNNN'
AmbiqDict = dict( zip( basePairs, codes ) )

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
                print( "{}, {}, {}".format( item, len( kmers[item] ), len( kmers[item][variant].parent ) ) )
        else:
            print( "{}, {}".format( item, len( kmers[item] ) ) )


file = "/Users/natem/Documents/Sequences/LassaL_Sierra_Leone_Alignment.fasta"
lassaAlignment = alignmentIO.parseAlignment( file )

for thing in lassaAlignment:
    thing.filterKmers( 63, 67 )

totalKmers = collapseKmers( lassaAlignment )

kmerOutput( totalKmers )
